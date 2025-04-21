use crate::types::MatchType;
use crate::annotate::flank::{FlankGroup, FlankSeq};
use rand::Rng;
use colored::Colorize;
use pa_bitpacking::search::*;
use std::time::Instant;
use crate::types::*;
use crate::annotate::tune::*;
use crate::simdedits::simd::{TransposedQueries, simd_search};

#[derive(Clone)]
pub struct BarMan {
    pub queries: Vec<FlankGroup>,
    pub q: f32, // default set to 0.4
    pub min_mask_available: f32,
    pub max_edits: u8,
    pub filter_overlap: f32, // filter overlapping matches
    pub fp_target: f64, // target false positive rate
    pub edit_distance_thresholds: Vec<i32>,
    pub n_tune_runs: usize,
}


impl BarMan {

    // const LOG_PROB_TUNE_RUNS: usize = 1_000_000;

    pub fn new(queries: Vec<FlankGroup>, q: f32, min_mask_available: f32, 
        filter_overlap: f32, fp_target: f64, max_edits: u8, n_tune_runs: usize) -> Self {
        Self { queries, q, min_mask_available, filter_overlap, 
            fp_target, max_edits, edit_distance_thresholds: Vec::new(), n_tune_runs }
    }

    fn create_barcode_match(&self, f_start: usize, f_end: usize, edits: u8, rel_dist: isize, query: &FlankGroup, mask_index: usize) -> Match {
        let label = query.flank_seq.mask_ids[mask_index].clone();
        Match::new(
            EncodedMatchStr::new(
                query.match_type.clone(), 
                query.orientation.clone(), 
                Some(label)
            ), 
            f_start, f_end, 
            Some(edits as i32), 
            rel_dist
        )
    }

    fn create_flank_match(&self, f_start: usize, f_end: usize, edit_dist: i32, rel_dist: isize, query: &FlankGroup) -> Match {
        Match::new(
            EncodedMatchStr::new(
                MatchType::Flank, 
                query.orientation.clone(), 
                Some("Flank".to_string()) // or maybe just None
            ), 
            f_start, f_end,
            Some(edit_dist),
            rel_dist
        )
    }

    pub fn annotate(&self, read: &[u8]) -> Vec<Match> {
        let mut all_matches = Vec::new();
        
        for (q_i, query) in self.queries.iter().enumerate() {
            let flank = &query.flank_seq;
            let query_cut_off = self.edit_distance_thresholds[q_i];
            let (passed_mask, locations, barcode_ranges) = self.locate_flank(flank, read, query_cut_off);

            for ((passed, &(f_start, f_end, edit_dist)), bar_range) in passed_mask.iter().zip(locations.iter()).zip(barcode_ranges.iter()) {
                let rel_dist = rel_dist_to_end(f_start as isize, read.len());
                if *passed {
                    // We only look  for the barcode in the mask region 
                    let (bar_start, bar_end) = bar_range.unwrap();
                    let barcode_slice = &read[bar_start.saturating_sub(5)..(bar_end + 5).min(read.len())];
                   // panic!("barcode_slice: {:?}", String::from_utf8_lossy(barcode_slice));
                    let mask_query_slices = flank.mask_queries.iter().map(|q| q.as_ref()).collect::<Vec<_>>();
                    if let Some((barcode_idx, edits)) = self.check_barcode(barcode_slice, &mask_query_slices) {
                        all_matches.push(self.create_barcode_match(f_start, f_end, edits, rel_dist, query, barcode_idx));
                    } else {
                        all_matches.push(self.create_flank_match(f_start, f_end, edit_dist, rel_dist, query));
                    }
                }
                else {
                    all_matches.push(self.create_flank_match(f_start, f_end, edit_dist, rel_dist, query));
                }
            }
        }
        self.final_assignment(&all_matches)
    }

 

    pub fn auto_tune_parmas(&mut self) {
        let start_time = Instant::now();
        println!("\n{}", "Parameter Tuning: Edit Distance Fraction".bold().underline());
        println!("  • Range: {} - {}", "0.10".dimmed(), "0.50".dimmed());
        println!("  • Test sequences: {}\n", "10,000".dimmed());

        // Tune edit distance thresholds
        self.edit_distance_thresholds = tune_edit_distance(&self.queries, self.n_tune_runs, self.fp_target);

        // Pretty print each threshold
        println!("{}", "Edit distance thresholds:".green().bold());
        for (i, threshold) in self.edit_distance_thresholds.iter().enumerate() {
            println!("  • Query {:>2}: {}", i, threshold.to_string().bold());
        }

        if self.max_edits == 0 {
            self.max_edits = tune_max_edits(&self.queries, self.n_tune_runs, self.fp_target);
        }


        // Print posterior threshold
        println!(
            "{} {}",
            format!("Max edits (FP: {}):", self.fp_target).green(),
            self.max_edits.to_string().bold()
        );

        // Done message
        println!(
            "{} {} seconds",
            "Done! Tuning took".bold(),
            start_time.elapsed().as_secs()
        );

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

            // Check if the current match overlaps with ANY match in the current group
            let overlaps_with_group = current_group.iter().any(|group_match| {
                let overlap_start = match_item.start.max(group_match.start);
                let overlap_end = match_item.end.min(group_match.end);
                
                if overlap_end <= overlap_start {
                    false
                } else {
                    let overlap_length = overlap_end - overlap_start;
                    let min_length = (match_item.end - match_item.start)
                        .min(group_match.end - group_match.start);
                    
                    // Calculate the actual fraction and then compare
                    let fraction = overlap_length as f32 / min_length as f32;
                    fraction >= self.filter_overlap
                }
            });

            if overlaps_with_group {
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
        // For barcode we select max log prob, for flanks min edits
        let final_matches: Vec<Match> = groups.into_iter()
        .map(|group| {
            group.iter()
                .filter(|m| matches!(m.label.match_type, MatchType::Fbarcode | MatchType::Rbarcode))
                .min_by(|a, b| {
                    a.edit_dist.unwrap().partial_cmp(&b.edit_dist.unwrap()).unwrap()
                })
                .or_else(|| {
                    group.iter()
                        .min_by_key(|m| m.edit_dist.unwrap())
                })
                .unwrap()
                .clone()
        })
        .collect();

        final_matches
    }


    fn split_scores_in_alns(&self, flank: &FlankSeq, result: &SearchResult, threshold: i32 ) -> (Vec<bool>, Vec<(usize, usize, i32)>, Vec<Option<(usize, usize)>>) {
        // Filter scores that pass threshold
        let passed_scores = result.out.iter()
            .enumerate()
            .filter(|&(_, edits)| *edits <= threshold)
            .map(|(pos, &edits)| (pos, edits))
            .collect::<Vec<_>>();

        // Find local minima
        let minima = find_prominent_minima(&passed_scores, 1.0);
       // println!("Minima found: {:?}", minima);
        
        // Pre-allocate result vectors
        let mut traced_ranges = Vec::with_capacity(minima.len());
        let mut passed_mask = Vec::with_capacity(minima.len());
        let mut barcode_ranges = Vec::with_capacity(minima.len());

        // Process each minimum
        for &min_idx in &minima {
            let (_, path_positions) = result.trace(min_idx);
            
            // Get alignment boundaries
            let start = path_positions.first().unwrap().0 as usize;
            let end = path_positions.last().unwrap().0 as usize - 1; // Cause aligns up to x, so ends at x-1
            
            // Check mask coverage
            let (_, mask_passed, r_range) = flank.mask_covered(&path_positions, self.min_mask_available);
            
            traced_ranges.push((start, end, result.out[min_idx])); // Covered read area, with number of edits 
            passed_mask.push(mask_passed); // Whether at least self.min_mask_available is covered
            barcode_ranges.push(r_range); // Range of barcode covered, in the read
        }

        (passed_mask, traced_ranges, barcode_ranges)
    }

     
    // return the index of the best match and the posterior probability
    fn check_barcode(&self, barcode_slice: &[u8], mask_fillers: &[&[u8]]) -> Option<(usize, u8)> {
        // We use the mutation rates to check which barcode is more probable and what it's posterior probability is
        // assuming equal priors 

        // Group mask filters in batches of 32, and encode in TransposedQueries
        let mut mask_transposed_queries = Vec::with_capacity(mask_fillers.len() / 32);
        for group in mask_fillers.chunks(32) {
            mask_transposed_queries.push(TransposedQueries::new(group.to_vec()));
        }

        // Vector to get all scores in, retains orignal order
        let edits = mask_transposed_queries.iter().map(
            |tq| 
            simd_search(tq, barcode_slice)
        ).flatten().collect::<Vec<_>>();

        let mut lowest_edits = u8::MAX;
        let mut lowest_idx = usize::MAX;
        let mut second_lowest_edits = u8::MAX;

        for (idx, edit) in edits.iter().enumerate() {
            if edit < &lowest_edits {
                second_lowest_edits = lowest_edits;
                lowest_edits = *edit;
                lowest_idx = idx;
            } else if edit < &second_lowest_edits {
                second_lowest_edits = *edit;
            }
        }

        if lowest_edits < self.max_edits {
            Some((lowest_idx, lowest_edits))
        } else {
            None
        }
    }
    

    fn locate_flank(&self, flank: &FlankSeq, read: &[u8], query_cut_off: i32) -> (Vec<bool>, Vec<(usize, usize, i32)>, Vec<Option<(usize, usize)>>) {
        // Perform semi-global alignment of flank, using q as gap penalty for prefix/suffix region
        let result = search(flank.seq.as_ref(), read, self.q);
        self.split_scores_in_alns(flank, &result, query_cut_off)
    }
}



pub fn random_dna(length: usize) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let bases: Vec<u8> = vec![b'A', b'T', b'C', b'G'];
    (0..length).map(|_| bases[rng.gen_range(0..bases.len())]).collect()
}

fn is_local_minimum(points: &[(usize, i32)], i: usize, window_size: usize) -> bool {
    let n = points.len();

    // println!("Checking if {} is a local minimum", i);
    
    // Look at left window
    let left_start = i.saturating_sub(window_size);
    for j in left_start..i {
        if points[j].0 + window_size >= points[i].0 && points[j].1 <= points[i].1 {
            // println!("Left window: {} is not a local minimum", i);
            return false;
        }
    }

    // Look at right window
    let right_end = (i + window_size + 1).min(n);
    for j in (i + 1)..right_end {
        if points[j].0 <= points[i].0 + window_size && points[j].1 < points[i].1 {
            // println!("Right window: {} is not a local minimum", i);
            return false;
        }
    }

    true
}

fn calculate_prominence(points: &[(usize, i32)], i: usize, window_size: usize) -> f32 {
    let n = points.len();

    // assert!(points.len() > window_size, "Window size is larger than the number of points");
    
    // Calculate left prominence
    let mut left_max = points[i].1;
    let mut j = i;
    while j > 0 && points[j - 1].0 >= points[i].0.saturating_sub(window_size) {
        j -= 1;
        left_max = left_max.max(points[j].1);
    }
    let left_prominence = (left_max - points[i].1) as f32;

    // Calculate right prominence
    let mut right_max = points[i].1;
    j = i;
    while j < n - 1 && points[j + 1].0 <= (points[i].0 + window_size).min(n - 1) {
        j += 1;
        right_max = right_max.max(points[j].1);
    }
    let right_prominence = (right_max - points[i].1) as f32;

    // Return the larger prominence value
    left_prominence.max(right_prominence)
}

fn find_prominent_minima(points: &[(usize, i32)], min_prominence: f32) -> Vec<usize> {
    let mut minima = Vec::new();
    let n = points.len();
    let window_size = 10;

    if n < 2 {
        return minima;
    }

    for i in 0..n {
        if is_local_minimum(points, i, window_size) {
            // println!("Local minimum at position: {}", points[i].0);
            let prominence = calculate_prominence(points, i, window_size);
            // println!("Prominence: {}", prominence);
            if prominence >= min_prominence {
                minima.push(points[i].0);
            }
        }
    }

    minima
}

pub fn rel_dist_to_end(pos: isize, read_len: usize) -> isize {
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


mod test {

    use super::*;

    #[test]
    fn test_rel_dist_to_end() {
        assert_eq!(rel_dist_to_end(0, 10), 1); // Technically 0, but left side would also be 0, but -0? so using 1 and -1 instead
        assert_eq!(rel_dist_to_end(5, 10), 5);
        assert_eq!(rel_dist_to_end(10, 10), -1);
    }

    #[test]
    fn test_find_prominent_minima() {
        // Generate 20 points with some clear minima
        let points = vec![
            (0, 10), (1, 8), (2, 6), (3, 4), (4, 2),  // Descending
            (5, 1), (6, 3), (7, 5), (8, 7), (9, 9),   // Valley at 5, then ascending
            (10, 10), (11, 8), (12, 6), (13, 4), (14, 2), // Descending again
            (15, 1), (16, 3), (17, 5), (18, 7), (19, 9)   // Valley at 15, then ascending
        ];
        
        let minima = find_prominent_minima(&points, 1.0);
        // You would expect 5 and 15, here however we compare to - window size to + window size
        // since window size is 10 position 15 is being compared to 5 which is also a minimum
        // hence 15 is not cosnidered a minimum. 5 however is.
        assert_eq!(minima, vec![5]); 
    }

    #[test]
    fn test_find_prominent_minima_on_left_edge() {
        let points = vec![
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5),  
            (5, 6), (6, 7), (7, 8), (8, 9), (9, 10),
            (10, 11), (11, 12), (12, 13), (13, 14), 
            (14, 15) 
        ];
        let minima = find_prominent_minima(&points, 1.0);
        assert_eq!(minima, vec![0]);
    }

    #[test]
    fn test_multiple_minima() {
    let points = vec![
            (0, 1), (1, 1), (2, 1), (3, 4), (4, 5),  
            (5, 6), (6, 1), (7, 8), (8, 9), (9, 10),
            (10, 11), (11, 12), (12, 13), (13, 14), 
            (14, 15) 
        ];
        let minima = find_prominent_minima(&points, 1.0);
        println!("Minima: {:?}", minima);
        assert_eq!(minima, vec![0]);
    }



    #[test]
    fn test_find_prominent_minima_on_right_edge() {
        let points = vec![
            (0, 15), (1, 14), (2, 13), (3, 12), (4, 11),  
            (5, 10), (6, 9), (7, 8), (8, 7), (9, 6),
            (10, 5), (11, 4), (12, 3), (13, 2), 
            (14, 1) 
        ];
        let minima = find_prominent_minima(&points, 1.0);
        assert_eq!(minima, vec![14]);
    }

    #[test]
    fn test_find_multiple_prominent_minima_sep_window_size() {
        let points = vec![
            (0, 1), // <- minimum 
            (1, 2), (2, 3), (3, 4), (4, 5),  
            (5, 6), (6, 7), (7, 8), (8, 9), (9, 10),
            (10, 11), (11, 12), (12, 13), (13, 14), (14, 15), 
            (15, 2),// <- minimum
            (16, 3), (17, 4), (18, 5),
            (19, 6), (20, 7), (21, 8), (22, 9), (23, 10),
            (24, 11), (25, 12), (26, 13), (27, 14), 
            (29, 15)
        ];
        let minima = find_prominent_minima(&points, 1.0);
        assert_eq!(minima, vec![0, 15]);
    }


    // #[test]
    // fn test_simple_barcode_search() {
    //     let flank = FlankSeq::init_from_sequences(&[
    //         b"GGGGAAACCCC".to_vec(), 
    //         b"GGGGTTTCCCC".to_vec()],
    //         vec!["seq1".to_string(), "seq2".to_string()]
    //     ).unwrap();

    //     let read = b"ATCTATCAGCATCGACGACTGGGGAAACCCCATCTACTAGCATCGACTAGCTAGATAGATAGTAGATAGATGAGA";
    //     //                                          GGGGAAACCCC
    //     //                                         20         30
    //     let bm = BarMan::new(vec![FlankGroup {
    //         flank_seq: flank,
    //         match_type: MatchType::Fbarcode,
    //         orientation: Orientation::Forward,
    //     }], 0.4, 0.5, -2.0,
    //      0.9, 20);
    //     let matches = bm.annotate(read);
    //     assert_eq!(matches.len(), 1);
    //     let expected_match = Match::new( 
    //         EncodedMatchStr { 
    //             match_type: MatchType::Fbarcode, 
    //             orientation: Orientation::Forward, 
    //             label: Some("seq1".to_string()) 
    //         },
    //         20, 
    //         30, 
    //         None, 
    //         20
    //     );
    //     println!("Matches: {:?}", matches);
    //     assert_eq!(matches[0], expected_match);
    // }


//     #[test]
//     fn test_simple_barcode_search_double_match() {
//         let flank = FlankSeq::init_from_sequences(&[
//             b"GGGGAAACCCC".to_vec(), 
//             b"GGGGTTTCCCC".to_vec()],
//             vec!["seq1".to_string(), "seq2".to_string()]
//         ).unwrap();

//         let read = b"ATCTATCAGCATCGACGACTGGGGAAACCCCATCTACTAGCATCGACTAGCTAGATGGGGTTTCCCCAGATAGTAGATAGATGAGA";
//         //                                          GGGGAAACCCC                         GGGGTTTCCCC
//         //                                         20         30                        56         66
//         let bm = BarMan::new(vec![FlankGroup {
//             flank_seq: flank,
//             match_type: MatchType::Fbarcode,
//             orientation: Orientation::Forward,
//         }], 0.4, 0.5, -2.0,
//          0.9, 20);
//         let matches = bm.annotate(read);
//         let expected_matches = vec![
//             Match::new(EncodedMatchStr { 
//                 match_type: MatchType::Fbarcode, 
//                 orientation: Orientation::Forward, 
//                 label: Some("seq1".to_string()) }, 
//                 20, 
//                 30, 
//                 Some(-0.2177120785045065), 
//                 20 ), 
//             Match::new(EncodedMatchStr { 
//                 match_type: MatchType::Fbarcode, 
//                 orientation: Orientation::Forward, 
//                 label: Some("seq2".to_string()) }, 
//                 56, 
//                 66, 
//                 Some(-0.2177120785045065), 
//                 None, 
//                 -30 )
//         ];
//         assert_eq!(matches, expected_matches);
//     }

//     #[test]
//     fn test_simple_barcode_search_double_close() {
//         let flank = FlankSeq::init_from_sequences(&[
//             b"GGGGAAACCCC".to_vec(), 
//             b"GGGGTTTCCCC".to_vec()],
//             vec!["seq1".to_string(), "seq2".to_string()]
//         ).unwrap();

//         let read = b"ATCTATCAGCATCGACGACTGGGGAAACCCCGGGGAAACCCCTAGATAGATGAGA";
//         //                                          GGGGAAACCCCGGGGAAACCCC                         
//         //                                         20         30         41
//         let bm = BarMan::new(vec![FlankGroup {
//             flank_seq: flank,
//             match_type: MatchType::Fbarcode,
//             orientation: Orientation::Forward,
//         }], 0.4, ErrorRates::default(), 0.5, 
//         -2.0, 0.9, 0.0);
//         let matches = bm.annotate(read);
//         let expected_matches = vec![
//             Match::new(EncodedMatchStr { 
//                 match_type: MatchType::Fbarcode, 
//                 orientation: Orientation::Forward,
//                 label: Some("seq1".to_string()) }, 
//                 20, 
//                 30, 
//                 Some(-0.2177120785045065), 
//                 None, 
//                 20 ), 
//             Match::new(EncodedMatchStr { 
//                 match_type: MatchType::Fbarcode, 
//                 orientation: Orientation::Forward, 
//                 label: Some("seq1".to_string()) 
//             }, 
//                 31, 
//                 41, 
//                 Some(-0.2177120785045065), 
//                 None, 
//                 -24 )
//         ];
//         assert_eq!(matches, expected_matches);
//     }

//     #[test]
//     fn test_match_and_just_flank() {
//         let flank = FlankSeq::init_from_sequences(&[
//             b"GGGGAAACCCCCCCCCCCCCC".to_vec(), 
//             b"GGGGTTTCCCCCCCCCCCCCC".to_vec()],
//             vec!["seq1".to_string(), "seq2".to_string()]
//         ).unwrap();

//         println!("Flank mask region: {:?}", flank.mask_region);

//         let read = b"ATCTATCAGCATCGACGACTGGGGAAACCCCCCCCTAGCTGGGGGGGGGGGCCCCCCCCCCCCCCC";
//         //                                          GGGGAAACCCCCCCC (partial)      CCCCCCCCCCCCCCC (just flank)      
//         //                                                                    GGGGAAACCCCCCCCCCCCCC       
//         let bm = BarMan::new(vec![FlankGroup {
//             flank_seq: flank,
//             match_type: MatchType::Fbarcode,
//             orientation: Orientation::Forward,
//         }], 0.4, ErrorRates::default(), 0.4,
//           -2.0, 0.9, 0.0); // -2 is very strict, but we  have short tests here
//                                                                                    // Only allowed matches with > 99% posterior prob
//         let matches = bm.annotate(read);

//         let expected_matches = [Match { read: None, label: 
//             EncodedMatchStr { 
//                 match_type: MatchType::Fbarcode, 
//                 orientation: Orientation::Forward, 
//                 label: Some("seq1".to_string()) 
//             }, 
//                 start: 20, 
//                 end: 38, 
//                 log_prob: Some(-0.2177120785045065), 
//                 edit_dist: None, 
//                 read_len: None, 
//                 rel_dist_to_end: 20, 
//                 record_set_idx: None, 
//                 record_idx: None, 
//                 cuts: None }, 
            
//             Match { read: None, label: 
//                 EncodedMatchStr { 
//                     match_type: MatchType::Flank, 
//                     orientation: Orientation::Forward, 
//                     label: Some("Flank".to_string()) }, 
//                 start: 44, 
//                 end: 64, 
//                 log_prob: None, 
//                 edit_dist: Some(0), 
//                 read_len: None, 
//                 rel_dist_to_end: -22, 
//                 record_set_idx: None, 
//                 record_idx: None, 
//                 cuts: None 
//             }
//         ];
//         assert_eq!(matches, expected_matches);
//     }

  
// }

// #[cfg(test)]
// mod overlap_tests {
//     use super::*;

//     fn create_test_matches() -> Vec<Match> {
//         vec![
//             // Case 1: Two matches with 90% overlap (second match extends 10% beyond first)
//             Match::new(
//                 EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("barcode1".to_string())),
//                 100, 200, Some(-5.0), None, 100
//             ),
//             Match::new(
//                 EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("barcode2".to_string())),
//                 110, 210, Some(-4.0), None, 110
//             ),
            
//             // Case 2: One match completely contained within another (A is within B)
//             Match::new(
//                 EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("barcode3".to_string())),
//                 300, 400, Some(-6.0), None, 300
//             ),
//             Match::new(
//                 EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("barcode4".to_string())),
//                 320, 380, Some(-3.0), None, 320
//             ),
            
//             // Case 3: Overruling - one match extends beyond the other on both sides
//             Match::new(
//                 EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("barcode5".to_string())),
//                 500, 600, Some(-7.0), None, 500
//             ),
//             Match::new(
//                 EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("barcode6".to_string())),
//                 480, 620, Some(-2.0), None, 480
//             ),
            
//             // Non-overlapping match to ensure separation works
//             Match::new(
//                 EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("barcode7".to_string())),
//                 700, 800, Some(-1.0), None, 700
//             ),
//         ]
//     }
    
//     #[test]
//     fn test_partial_overlap_90_percent() {
//         // Create a BarMan with filter_overlap = 0.9 (90%)
//         let bar_man = BarMan::new(
//             vec![], 0.4, ErrorRates::default(), 
//             0.5, 0.0, 0.9, 0.0
//         );
        
//         let test_matches = create_test_matches();
        
//         // Filter to just the first test case
//         let case1_matches = test_matches.iter()
//             .filter(|m| m.start == 100 || m.start == 110)
//             .cloned()
//             .collect::<Vec<_>>();
            
//         let final_matches = bar_man.final_assignment(&case1_matches);
        
//         // Should be grouped and we pick the higher log prob (-4.0)
//         assert_eq!(final_matches.len(), 1);
//         assert_eq!(final_matches[0].label.label, Some("barcode2".to_string()));
//         assert_eq!(final_matches[0].log_prob, Some(-4.0));
//     }
    
//     #[test]
//     fn test_contained_match() {
//         // Create a BarMan with filter_overlap = 0.9 (90%)
//         let bar_man = BarMan::new(
//             vec![], 0.4, ErrorRates::default(), 
//             0.5, 0.0, 0.9, 0.0
//         );
        
//         let test_matches = create_test_matches();
        
//         // Filter to just the second test case
//         let case2_matches = test_matches.iter()
//             .filter(|m| m.start == 300 || m.start == 320)
//             .cloned()
//             .collect::<Vec<_>>();
            
//         let final_matches = bar_man.final_assignment(&case2_matches);
        
//         // Should be grouped and we pick the higher log prob (-3.0)
//         assert_eq!(final_matches.len(), 1);
//         assert_eq!(final_matches[0].label.label, Some("barcode4".to_string()));
//         assert_eq!(final_matches[0].log_prob, Some(-3.0));
//     }
    
//     #[test]
//     fn test_overruling_match() {
//         // Create a BarMan with filter_overlap = 0.9 (90%)
//         let bar_man = BarMan::new(
//             vec![], 0.4, ErrorRates::default(), 
//             0.5, 0.0, 0.9, 0.0
//         );
        
//         let test_matches = create_test_matches();
        
//         // Filter to just the third test case
//         let case3_matches = test_matches.iter()
//             .filter(|m| m.start == 500 || m.start == 480)
//             .cloned()
//             .collect::<Vec<_>>();
            
//         let final_matches = bar_man.final_assignment(&case3_matches);
        
//         // Should be grouped and we pick the higher log prob (-2.0)
//         assert_eq!(final_matches.len(), 1);
//         assert_eq!(final_matches[0].label.label, Some("barcode6".to_string()));
//         assert_eq!(final_matches[0].log_prob, Some(-2.0));
//     }
    
//     #[test]
//     fn test_all_overlap_cases_together() {
//         // Create a BarMan with filter_overlap = 0.9 (90%)
//         let bar_man = BarMan::new(
//             vec![], 0.4, ErrorRates::default(), 
//             0.5, 0.0, 0.9, 0.0
//         );
        
//         let test_matches = create_test_matches();
//         let final_matches = bar_man.final_assignment(&test_matches);
        
//         // Should have 4 matches total (3 groups + 1 singleton)
//         assert_eq!(final_matches.len(), 4);
        
//         // Check that each group kept the highest log prob match
//         let mut has_barcode2 = false;
//         let mut has_barcode4 = false;
//         let mut has_barcode6 = false;
//         let mut has_barcode7 = false;
        
//         for m in &final_matches {
//             match m.label.label.as_deref() {
//                 Some("barcode2") => has_barcode2 = true,
//                 Some("barcode4") => has_barcode4 = true,
//                 Some("barcode6") => has_barcode6 = true,
//                 Some("barcode7") => has_barcode7 = true,
//                 _ => {}
//             }
//         }
        
//         assert!(has_barcode2, "Missing barcode2 which should be selected from group 1");
//         assert!(has_barcode4, "Missing barcode4 which should be selected from group 2");
//         assert!(has_barcode6, "Missing barcode6 which should be selected from group 3");
//         assert!(has_barcode7, "Missing barcode7 which should be a singleton");
//     }
    
//     #[test]
//     fn test_identical_matches() {
//         // Create two identical matches that should be grouped
//         let matches = vec![
//             Match::new(
//                 EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("barcode1".to_string())),
//                 100, 200, Some(-5.0), None, 100
//             ),
//             Match::new(
//                 EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("barcode2".to_string())),
//                 100, 200, Some(-4.0), None, 100
//             ),
//         ];
        
//         let bar_man = BarMan::new(
//             vec![], 0.4, ErrorRates::default(), 
//             0.5, 0.0, 0.9, 0.0
//         );
        
//         let final_matches = bar_man.final_assignment(&matches);
        
//         // Should be 1 match with the better log prob
//         assert_eq!(final_matches.len(), 1);
//         assert_eq!(final_matches[0].label.label, Some("barcode2".to_string()));
//     }
// }
}