use crate::barbell::seq::*;
use crate::barbell::reader::*;
use crate::barbell::misc::*;
use pa_bitpacking::search::*;


fn rel_dist_to_end(pos: isize, read_len: usize) -> isize {
    // If already negative, for a match starting before the read we can return 1
    if pos < 0 {
        return 1;
    }

    if pos <= (read_len / 2) as isize {
        if pos == 0 {
            1  // Left end
        } else {
            pos  // Distance from left end (already isize)
        }
    } else if pos == read_len as isize {
        -1  // Right end
    } else {
        -(read_len as isize - pos)  // Distance from right end
    }
}

pub trait Strategy {
    fn default() -> Self;
    fn annotate(&self, flanks: &[QueryGroup], read: &[u8]) -> Vec<Match>;
    fn final_assignment(&self, matches: &[Match]) -> Vec<Match>;
}

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
       let mut all_edits = Vec::new();

       for (i, query) in query_groups.iter().enumerate() {
           // Always flank?
           let flank = query.flank.as_ref().unwrap();
           let (passed_mask, locations, edits) = self.locate_flank(flank, read);
           all_edits.push((i, edits));
           for (passed, (f_start, f_end, e)) in passed_mask.iter().zip(locations.iter()) {
               if *passed {
                  // println!("\t✅ Flank found: {:?} with edits: {}", (f_start, f_end), e);
                   let barcodes = self.find_valid_barcode((*f_start, *f_end), read, &query.queries, query.prefix_char, query.rc);
                   // for b in barcodes.iter() {
                   //     println!("\tBarcode: {:?}", b);
                   // }
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

       // for m in all_matches.iter() {
       //     println!("Match: {:?}", m);
       // }

     self.plot_edit_distances(all_edits);

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
}

impl SimpleStrategy {


   fn segregate_scores_in_alns(&self, flank: &FlankSeq, res: &SearchResult, threshold: i32) -> (Vec<bool>, Vec<(usize, usize, i32)>) {
     // The first thing we can do is discard all scores ABOVE our edit distance threshold
     // for clarify now, collect as tuple of (pos, edits)
     // println!("Using threshold: {}", threshold);
     let filtered_scores: Vec<(usize, i32)> = res.out.iter()
         .enumerate()
         .filter(|&(_, edits)| *edits <= threshold)
         .map(|(pos, &edits)| (pos, edits))
         .collect();
      // In the filtered scores, we now have to find local minima (misc.rs)
      let minima = find_prominent_minima(&filtered_scores, 1.0);
     // println!("Found minima: {:?}", minima);
      // Now for each minimum we do a traceback, then again remove overlapping (but tracebacked alignments)
       let mut traced_ranges: Vec<(usize, usize, i32)> = Vec::with_capacity(minima.len());
       let mut passed_mask = Vec::with_capacity(minima.len());
       for min_idx in minima {
           let (_, path_positions) = res.trace(min_idx);
           let left_side = path_positions.first().unwrap();
           let right_side = path_positions.last().unwrap();
           let (t_from, t_to) = (left_side.0, right_side.0);
           // println!("\t {} - {} with edits: {}", t_from, t_to, res.out[min_idx]);
           let (mask_coverage, passed) = flank.mask_covered(&path_positions[..], self.min_mask_available);
           traced_ranges.push((t_from as usize, t_to as usize, res.out[min_idx]));
           passed_mask.push(passed);
       }
       // Only merge if contained within another
       //let overlap_removed = remove_overlap(traced_ranges);
       //println!("Overlap removed: {:?}", overlap_removed);
     //  println!("Traced ranges: {:?}", traced_ranges);
       (passed_mask, traced_ranges)
   }


   fn locate_flank(&self, flank: &FlankSeq, read: &[u8]) -> (Vec<bool>, Vec<(usize, usize, i32)>, Vec<i32>) {
       // Perform semi-global alignment of flank, using q as gap penalty for prefix/suffix region
       let result = search(flank.seq.as_ref(), read, self.q);
       let half_mask = flank.mask_len() / 2; // we count only half as expected edits
       let threshold = ((flank.unmasked_len() + half_mask) as f32 * self.max_edit_fraction) as i32;
       println!("Using threshold: {} from unmasked len: {}", threshold, flank.unmasked_len());
       let (passed_mask, locations) = self.segregate_scores_in_alns(flank, &result, threshold);
       let edits = result.out;
       (passed_mask, locations, edits)
   }

   fn find_valid_barcode(&self, flank:(usize, usize), read: &[u8], queries: &[Query], prefix_char: u8, rc: bool) -> Option<Match> {
       let mut barcodes = Vec::new();
       for (query_idx, query) in queries.iter().enumerate() {
           let result = search(query.seq.as_ref(), read[flank.0..flank.1].as_ref(), self.q);
           let threshold = (query.seq.len() as f32 * self.max_edit_fraction) as i32;
           let min_score = result.out.iter().min().unwrap();
           if *min_score <= threshold {
               // println!("\t✅ Barcode passed: {:?} with edits: {}", (flank.0, flank.1), *min_score);
               let label = format!("{}#{}_{}", query.id, prefix_char as char, if rc { "rc" } else { "fw" });
               let rel_dist = rel_dist_to_end(flank.0 as isize, read.len());
               barcodes.push(Match::new(label, MatchType::Barcode, flank.0, flank.1, *min_score, rel_dist));
           }
       }
       // Sort barcode matches by edit distance (low to high)
       barcodes.sort_by_key(|m| m.edits);


       // If we have at least two matches, we can check if the second best is at least min_barcode_edit_diff better than the best
       if barcodes.len() >= 2 {
           let best_edit = barcodes[0].edits;
           let second_best_edit = barcodes[1].edits;
          // println!("Best edit: {}, Second best edit: {}", best_edit, second_best_edit);
           if second_best_edit - best_edit >= self.min_barcode_edit_diff {
               Some(barcodes[0].clone())
           } else {
               None
           }
       } else if barcodes.len() == 1 {
           Some(barcodes[0].clone())
       } else {
           None
       }
       
       
   }



   fn plot_edit_distances(&self, all_edits: Vec<(usize, Vec<i32>)>) -> Result<(), Box<dyn std::error::Error>> {
       use plotly::{Plot, Scatter};
       
       let mut plot = Plot::new();
       
       // Create a scatter plot for each flank's edits
       for (flank_idx, edits) in all_edits.iter() {
           let x: Vec<f64> = (0..edits.len()).map(|i| i as f64).collect();
           let y: Vec<f64> = edits.iter().map(|&e| e as f64).collect();
           
           let trace = Scatter::new(x, y)
               .name(format!("Flank {}", flank_idx))
               .mode(plotly::common::Mode::Lines);
           
           plot.add_trace(trace);
       }

       // Calculate and plot difference between best and second-best scores
       if all_edits.len() >= 2 {
           let seq_len = all_edits[0].1.len();
           let mut diff_y: Vec<f64> = Vec::with_capacity(seq_len);
           let mut intersections_x: Vec<f64> = Vec::new();
           let mut intersections_y: Vec<f64> = Vec::new();
           
           // For each position
           for pos in 0..seq_len {
               // Get all edit distances at this position
               let mut scores: Vec<i32> = all_edits.iter()
                   .map(|(_, edits)| edits[pos])
                   .collect();
               scores.sort();
               
               // Calculate difference between best and second-best
               let diff = scores[1] - scores[0];
               diff_y.push(diff as f64);

               // If the best two scores are equal, this is an intersection point
               if scores[0] == scores[1] {
                   intersections_x.push(pos as f64);
                   intersections_y.push(scores[0] as f64);
               }
           }

         //  let x: Vec<f64> = (0..seq_len).map(|i| i as f64).collect();
        //    let diff_trace = Scatter::new(x, diff_y)
        //        .name("Difference (Second Best - Best)")
        //        .mode(plotly::common::Mode::Lines);
           
           // Add intersection points
           let intersection_trace = Scatter::new(intersections_x, intersections_y)
               .name("Intersections")
               .mode(plotly::common::Mode::Markers)
               .marker(plotly::common::Marker::new().color("pink").size(5));

           // plot.add_trace(diff_trace);
           plot.add_trace(intersection_trace);
       }
       
       // Customize the layout
       plot.set_layout(plotly::Layout::new()
           .title(plotly::common::Title::from("Edit Distances Across Read Positions"))
           .x_axis(plotly::layout::Axis::new().title(plotly::common::Title::from("Position")))
           .y_axis(plotly::layout::Axis::new().title(plotly::common::Title::from("Edit Distance"))));
       
       plot.show();
       
       Ok(())
   }
}   
