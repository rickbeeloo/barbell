use crate::barbell::seq::*;
use crate::barbell::reader::*;
use crate::barbell::misc::*;
use pa_bitpacking::search::*;
use std::collections::HashMap;
use rand::Rng;
use rayon::prelude::*;
use spinners::{Spinner, Spinners};
use std::time::Instant;
use colored::Colorize;

pub trait Strategy {
    fn default() -> Self;
    fn annotate(&self, queries: &[QueryGroup], read: &[u8]) -> Vec<Match>;
    fn final_assignment(&self, matches: &[Match]) -> Vec<Match>;
    fn auto_tune_parmas(&mut self, queries: &[QueryGroup]);
}

#[derive(Clone)]
pub struct SimpleStrategy {
   q: f32,                      // Prefix/suffix penalty semi-global aln (0.5)
   max_edit_fraction: f32,      // Fraction of total (unmasked)query length to consider a match valid (0.4)
   min_barcode_edit_diff: i32,  // Minimum edit distance between best and second best barcode match (6)
   min_mask_available: f32,     // Minimum fraction of mask available to have sufficient space for barcode check(0.5)
   filter_overlap: f32,         // Minimum fraction of overlap between matches to consider competing (90%)
}

impl Strategy for SimpleStrategy {
   
   fn default() -> Self {
       Self { 
           q: 0.5, 
           max_edit_fraction: 0.35, 
           min_barcode_edit_diff: 3, 
           min_mask_available: 0.5,
           filter_overlap: 0.9,
       }
   }



   fn annotate(&self, query_groups: &[QueryGroup], read: &[u8]) -> Vec<Match> {
       let mut all_matches = Vec::new();
       
       for query in query_groups.iter() {
           let flank = query.flank.as_ref().unwrap(); // Flank gauranteed for now
           let (passed_mask, locations) = self.locate_flank(flank, read);
           for (passed, (f_start, f_end, e)) in passed_mask.iter().zip(locations.iter()) {
               if *passed {
                  // println!("\t✅ Flank found: {:?} with edits: {}", (f_start, f_end), e);
                   let barcodes = self.find_valid_barcode((*f_start, *f_end), read, &query.queries, query.prefix_char, query.rc);
                   if let Some(barcode) = barcodes {
                       all_matches.push(barcode);
                   } else {
                     //   println!("\t❌ no barcode found");
                       let rel_dist = rel_dist_to_end(*f_start as isize, read.len());
                       all_matches.push(Match::new("Flank".to_string(), MatchType::Flank, *f_start, *f_end, *e, rel_dist));
                   }
               } else {
                   // We still push the Flank as it hints contamination regardless of barcode presence
                  // println!("\t❌ no barcode region found: {:?} with edits: {}", (f_start, f_end), e);
                   let rel_dist = rel_dist_to_end(*f_start as isize, read.len());
                   all_matches.push(Match::new("Flank".to_string(), MatchType::Flank, *f_start, *f_end, *e, rel_dist));
               }
           }
       }
       all_matches
   }


   fn final_assignment(&self, matches: &[Match]) -> Vec<Match> {
       // Sort matches based on start position
       let mut sorted_matches = matches.to_vec();
       sorted_matches.sort_by_key(|m| m.start);

       // Group overlapping matches
       let mut groups: Vec<Vec<Match>> = Vec::new();
       let mut current_group: Vec<Match> = Vec::new();

       for match_item in sorted_matches {
           if current_group.is_empty() {
               current_group.push(match_item);
               continue;
           }

           let last_match = current_group.last().unwrap();
           let overlap_fraction = {
               let overlap_start = match_item.start.max(last_match.start);
               let overlap_end = match_item.end.min(last_match.end);
               if overlap_end <= overlap_start {
                   0.0
               } else {
                   let overlap_length = overlap_end - overlap_start;
                   let min_length = (match_item.end - match_item.start)
                       .min(last_match.end - last_match.start);
                   overlap_length as f32 / min_length as f32
               }
           };

           if overlap_fraction >= self.filter_overlap {
               current_group.push(match_item);
           } else {
               groups.push(current_group);
               current_group = vec![match_item];
           }
       }
       
       // Don't forget to add the last group
       if !current_group.is_empty() {
           groups.push(current_group);
       }

       // For each group, prefer Barcode over Flank
       let final_matches: Vec<Match> = groups.into_iter()
           .map(|group| {
               // println!("Group: {:?}", group);
               // First try to find the Barcode match with lowest edit distance
               group.iter()
                   .filter(|m| matches!(m.match_type, MatchType::Barcode))
                   .min_by_key(|m| m.edits)
                   .or_else(|| {
                       // If no Barcode, take the Flank with lowest edit distance
                       group.iter()
                           .min_by_key(|m| m.edits)
                   })
                   .unwrap()
                   .clone()
           })
           .collect();

     
       final_matches
   }


  

   fn auto_tune_parmas(&mut self, queries: &[QueryGroup])  {
        let start_time = Instant::now();
        // Generate sequence from 0.05 to 0.5 with 0.01 steps
        let param_values: Vec<f32> = (0..46).map(|i| 0.05 + (i as f32 * 0.01)).collect();
    
        let total_tests = 10_000;
        let tests_per_param = total_tests / param_values.len();

        println!("\n{}", "Parameter Tuning".bold().underline());
        println!("  • Range: {} - {}", "0.05".dimmed(), "0.50".dimmed());
        println!("  • Test sequences: {}\n", "10,000".dimmed());

        let mut sp = Spinner::new(Spinners::OrangeBluePulse, "Tuning on random sequences".into());
    
        // Track counts of number of barcodes found for each parameter value
        let param_counts: HashMap<i32, HashMap<usize, usize>> = param_values.par_iter()
            .map(|&param_value| {
                let mut local_counts: HashMap<usize, usize> = HashMap::new();
                let mut rng = rand::thread_rng();
                let bases = b"ACGT";
                let seq_length = 1000;
    
                for _ in 0..tests_per_param {
                    // Generate random sequence
                    let seq: Vec<u8> = (0..seq_length)
                        .map(|_| bases[rng.gen_range(0..4)])
                        .collect();
    
                    // Clone the necessary parts of self
                    let mut self_clone = self.clone();
                    self_clone.max_edit_fraction = param_value;
                    let annotation = self_clone.annotate(queries, &seq);
                    let matches = self_clone.final_assignment(&annotation);
                    let barcode_count = matches.len();
    
                    local_counts
                        .entry(barcode_count)
                        .and_modify(|count| *count += 1)
                        .or_insert(1);
                }
    
                let param_key = (param_value * 100.0) as i32;
                (param_key, local_counts)
            })
            .collect();
    
        sp.stop();
        // Find the most permissive parameter where we still get zero matches
        let best_param = param_values.iter()
            .rev() // Reverse to start from highest (most permissive) value
            .find(|&&param| {
                let param_key = (param * 100.0) as i32;
                let counts = param_counts.get(&param_key).unwrap();
                let false_positives: usize = counts.iter()
                    .filter(|(&k, _)| k > 0)  // Only count sequences with matches
                    .map(|(_, count)| *count)
                    .sum();
                false_positives == 0
            })
            .unwrap();
    
        println!("  • Selected: {}\n", best_param.to_string().bold());
        println!("Done! tuning took {} seconds ", start_time.elapsed().as_secs());
    
        self.max_edit_fraction = *best_param;
    }

}

impl SimpleStrategy {

   
   fn segregate_scores_in_alns(&self, flank: &FlankSeq, res: &SearchResult, threshold: i32) -> (Vec<bool>, Vec<(usize, usize, i32)>) {
       // Filter scores below threshold and collect positions with their edit distances
       let filtered_scores = res.out.iter()
           .enumerate()
           .filter(|&(_, edits)| *edits <= threshold)
           .map(|(pos, &edits)| (pos, edits))
           .collect::<Vec<_>>();

       // Find local minima in filtered scores
       let minima = find_prominent_minima(&filtered_scores, 1.0);

       // Pre-allocate vectors for results
       let mut traced_ranges = Vec::with_capacity(minima.len());
       let mut passed_mask = Vec::with_capacity(minima.len());

       // Process each minimum to get alignment ranges and mask coverage
       for &min_idx in &minima {
           let (_, path_positions) = res.trace(min_idx);
           
           // Extract alignment boundaries from path
           let alignment_range = {
               let start: i32 = path_positions.first().unwrap().0;
               let end = path_positions.last().unwrap().0;
               (start as usize, end as usize)
           };

           // Check mask coverage and store results
           let (_, mask_passed) = flank.mask_covered(&path_positions, self.min_mask_available);
           traced_ranges.push((alignment_range.0, alignment_range.1, res.out[min_idx]));
           passed_mask.push(mask_passed);
       }

       (passed_mask, traced_ranges)
   }


   fn locate_flank(&self, flank: &FlankSeq, read: &[u8]) -> (Vec<bool>, Vec<(usize, usize, i32)>) {
       // Perform semi-global alignment of flank, using q as gap penalty for prefix/suffix region
       let result = search(flank.seq.as_ref(), read, self.q);

       // We have a masked region that matches  anything, on average we expect around 50% edits when actually alignen
       // this region, so we take mask.len() / 2 as expected edits
       let half_mask = flank.mask_len() / 2; // we count only half as expected edits
       let threshold = ((flank.unmasked_len() + half_mask) as f32 * self.max_edit_fraction) as i32;

        // Use the threshold for filtering
       let (passed_mask, locations) = self.segregate_scores_in_alns(flank, &result, threshold);

       (passed_mask, locations)
   }

   fn select_best_match(&self, barcodes: &mut [Match]) -> Option<Match> {
       // Early return if empty
       if barcodes.is_empty() {
           return None;
       }

       // Sort matches by edit distance
       barcodes.sort_by_key(|m| m.edits);

       // Return single match if it's the only one
       if barcodes.len() == 1 {
           return Some(barcodes[0].clone());
       }

       // Check if best match is significantly better than second best
       let (best, second_best) = (&barcodes[0], &barcodes[1]);
       match second_best.edits - best.edits >= self.min_barcode_edit_diff {
           true => Some(best.clone()),
           false => None,
       }
   }

   fn find_valid_barcode(&self, flank: (usize, usize), read: &[u8], queries: &[Query], prefix_char: u8, rc: bool) -> Option<Match> {
       let mut barcodes = Vec::new();
       let read_target_slice = &read[flank.0..flank.1];

       // Try to match each query sequence against the flanked region
       for query in queries {
           let result = search(query.seq.as_ref(), read_target_slice, self.q);
           let threshold = (query.seq.len() as f32 * self.max_edit_fraction) as i32; // Should we precalculate this?

           // Check if the best alignment score is within our threshold
           if let Some(&min_score) = result.out.iter().min() {
               if min_score <= threshold {
                   let orientation = if rc { "rc" } else { "fw" };
                   let label = format!("{}#{}_{}", query.id, prefix_char as char, orientation);
                   let rel_dist = rel_dist_to_end(flank.0 as isize, read.len());
                   
                   barcodes.push(Match::new(
                       label,
                       MatchType::Barcode,
                       flank.0,
                       flank.1,
                       min_score,
                       rel_dist
                   ));
               }
           }
       }

       self.select_best_match(&mut barcodes)
   }


}   
