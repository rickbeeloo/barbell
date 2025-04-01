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
use crate::barbell::plot::*;

use crate::barbell::pattern_assign::*;

pub trait Strategy {
    fn default() -> Self;
    fn annotate(&self, queries: &[QueryGroup], read: &[u8]) -> Vec<Match>;
    fn final_assignment(&self, matches: &[Match]) -> Vec<Match>;
    fn auto_tune_parmas(&mut self, queries: &[QueryGroup]);
    fn generate_plot(&self, query_groups: &[QueryGroup], read: &[u8]);
}

#[derive(Clone)]
pub struct SimpleStrategy {
    q: f32,                      // Semi-global alignment prefix/suffix penalty
    max_edit_fraction: f32,      // Maximum allowed edit distance as fraction of query length
    min_barcode_edit_diff: i32,  // Required edit distance between best and second-best barcode
    min_mask_available: f32,     // Minimum required mask fraction for barcode validation
    filter_overlap: f32,         // Minimum overlap fraction to consider matches competing
}

impl Strategy for SimpleStrategy {
    fn default() -> Self {
        Self { 
            q: 0.4, 
            max_edit_fraction: 0.35, 
            min_barcode_edit_diff: 6, 
            min_mask_available: 0.5,
            filter_overlap: 0.9,
        }
    }

    fn generate_plot(&self, query_groups: &[QueryGroup], read: &[u8])  {
        let mut all_edits = Vec::new();
        for (query_idx, query) in query_groups.iter().enumerate() {
            let flank = query.flank.as_ref().unwrap(); // Flank gauranteed for now
            let result = search(flank.seq.as_ref(), read, self.q);
            all_edits.push((query_idx, result.out.clone()));
        }
        plot_edit_distances(all_edits).expect("Failed to plot edit distnaces!");
    }

    fn annotate(&self, query_groups: &[QueryGroup], read: &[u8]) -> Vec<Match> {
        let mut all_matches = Vec::new();
        
        for query in query_groups {
            let flank = query.flank.as_ref().expect("Flank sequence required");
            let (passed_mask, locations) = self.locate_flank(flank, read);

            for (passed, &(f_start, f_end, edit_dist)) in passed_mask.iter().zip(locations.iter()) {
                let rel_dist = rel_dist_to_end(f_start as isize, read.len());
                
                if *passed {
                    if let Some(barcode) = self.find_valid_barcode((f_start, f_end), read, query) {
                        all_matches.push(barcode);
                    } else {
                        all_matches.push(self.create_flank_match(f_start, f_end, edit_dist, rel_dist, query));
                    }
                } else {
                    all_matches.push(self.create_flank_match(f_start, f_end, edit_dist, rel_dist, query));
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
        // Should we look at the distance here too?
        let final_matches: Vec<Match> = groups.into_iter()
            .map(|group| {
           //     println!("Group: {:?}", group);
               // First try to find the Barcode match with lowest edit distance
               group.iter()
                   .filter(|m| matches!(m.match_str.match_type, MatchType::Fbarcode | MatchType::Rbarcode))
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
         // Generate sequence from 0.1 to 0.5 with 0.01 steps
         let param_values: Vec<f32> = (0..41).map(|i| 0.1 + (i as f32 * 0.01)).collect();
     
         let total_tests = 10_000;
         let tests_per_param = total_tests / param_values.len();

         println!("\n{}", "Parameter Tuning".bold().underline());
         println!("  • Range: {} - {}", "0.10".dimmed(), "0.50".dimmed());
         println!("  • Test sequences: {}\n", "10,000".dimmed());

         let mut sp = Spinner::new(Spinners::OrangeBluePulse, "Tuning on random sequences".into());
     
         // Track counts of number of barcodes found for each parameter value
         let param_counts: HashMap<i32, HashMap<usize, usize>> = param_values.par_iter()
             .map(|&param_value| {
                 let mut local_counts: HashMap<usize, usize> = HashMap::new();
                 let mut rng = rand::thread_rng();
                 let bases = b"ACGT";
                 let seq_length = 10_000;
     
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
         // self.max_edit_fraction = 0.28; // REMOVE
     }
}

impl SimpleStrategy {
    fn create_flank_match(&self, start: usize, end: usize, edits: i32, rel_dist: isize, query: &QueryGroup) -> Match {
        let match_str = EncodedMatchStr::new(MatchType::Flank, query.orientation.clone(), None);
        Match::new(match_str, start, end, edits, rel_dist)
    }

    fn find_valid_barcode(&self, flank: (usize, usize), read: &[u8], query_group: &QueryGroup) -> Option<Match> {
        let read_slice = &read[flank.0..flank.1];
        let mut barcodes: Vec<Match> = query_group.queries.iter()
            .filter_map(|query| {
                let result = search(query.seq.as_ref(), read_slice, self.q);
                let threshold = (query.seq.len() as f32 * self.max_edit_fraction) as i32;
                
                result.out.iter().min().and_then(|&min_score| {
                    if min_score <= threshold {
                        let rel_dist = rel_dist_to_end(flank.0 as isize, read.len());
                        let match_str = EncodedMatchStr::new(
                            query_group.match_type.clone(),
                            query_group.orientation.clone(),
                            Some(query.id.clone())
                        );
                        Some(Match::new(match_str, flank.0, flank.1, min_score, rel_dist))
                    } else {
                        None
                    }
                })
            })
            .collect();

        self.select_best_match(&mut barcodes)
    }

    fn segregate_scores_in_alns(&self, flank: &FlankSeq, res: &SearchResult, threshold: i32) -> (Vec<bool>, Vec<(usize, usize, i32)>) {
        // println!("Threshold: {}", threshold);

        // Filter scores below threshold and collect positions with their edit distances
        let filtered_scores = res.out.iter()
            .enumerate()
            .filter(|&(_, edits)| *edits <= threshold)
            .map(|(pos, &edits)| (pos, edits))
            .collect::<Vec<_>>();

        // Find local minima in filtered scores
        let minima = find_prominent_minima(&filtered_scores, 1.0);
        // println!("Minima: {:?}", minima);

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
        //let half_mask = flank.mask_len() / 2; // we count only half as expected edits

        // The edits we would normally expect considering unmasked sequence, len * max frac
        let threshold = (flank.unmasked_len()) as f32 * self.max_edit_fraction;
        // The mask will always match and thereby inherently recude the match rate by ~ half the mask 
        //threshold += half_mask as f32;

         // Use the threshold for filtering
        let (passed_mask, locations) = self.segregate_scores_in_alns(flank, &result, threshold as i32);

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
}   
