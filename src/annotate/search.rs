use crate::annotate::flank::{FlankGroup, FlankSeq};
use crate::annotate::tune::*;
use crate::types::MatchType;
use crate::types::*;
use colored::Colorize;
use rand::Rng;
use sassy::profiles::*;
use sassy::search::*;
use std::cell::RefCell;
use std::sync::{Arc, RwLock};
use std::time::Instant;

#[derive(Clone)]
pub struct BarMan {
    pub queries: Vec<FlankGroup>,
    pub min_mask_available: f32,
    pub max_edits: i32,
    pub filter_overlap: f32, // filter overlapping matches
    pub fp_target: f64,      // target false positive rate
    pub edit_distance_thresholds: Vec<i32>,
    pub n_tune_runs: usize,
    pub q: f32, // Store q parameter for creating new searchers
}

impl BarMan {
    // const LOG_PROB_TUNE_RUNS: usize = 1_000_000;

    pub fn new(
        queries: Vec<FlankGroup>,
        q: f32,
        min_mask_available: f32,
        filter_overlap: f32,
        fp_target: f64,
        max_edits: i32,
        n_tune_runs: usize,
    ) -> Self {
        Self {
            queries,
            min_mask_available,
            filter_overlap,
            fp_target,
            max_edits,
            edit_distance_thresholds: Vec::new(),
            n_tune_runs,
            q,
        }
    }

    pub fn create_searcher(&self) -> Searcher<Iupac> {
        Searcher::new_fwd_with_overhang(self.q)
    }

    fn create_barcode_match(
        &self,
        f_start: usize,
        f_end: usize,
        edits: u8,
        rel_dist: isize,
        query: &FlankGroup,
        mask_index: usize,
    ) -> BarbellMatch {
        let label = query.flank_seq.mask_ids[mask_index].clone();
        BarbellMatch::new(
            EncodedMatchStr::new(
                query.match_type.clone(),
                query.orientation.clone(),
                Some(label),
            ),
            f_start,
            f_end,
            Some(edits as i32),
            rel_dist,
        )
    }

    fn create_flank_match(
        &self,
        f_start: usize,
        f_end: usize,
        edit_dist: i32,
        rel_dist: isize,
        query: &FlankGroup,
    ) -> BarbellMatch {
        BarbellMatch::new(
            EncodedMatchStr::new(
                MatchType::Flank,
                query.orientation.clone(),
                Some("Flank".to_string()), // or maybe just None?
            ),
            f_start,
            f_end,
            Some(edit_dist),
            rel_dist,
        )
    }

    pub fn annotate(&self, read: &[u8]) -> Vec<BarbellMatch> {
        let mut searcher = self.create_searcher();
        let mut all_matches = Vec::new();
        let static_read = StaticText::new(read);

        for (q_i, query) in self.queries.iter().enumerate() {
            let flank = &query.flank_seq;
            let query_cut_off = self.edit_distance_thresholds[q_i];
            let (passed_mask, locations, barcode_ranges) =
                self.locate_flank(&mut searcher, flank, &static_read, query_cut_off);

            for ((passed, &(f_start, f_end, edit_dist)), bar_range) in passed_mask
                .iter()
                .zip(locations.iter())
                .zip(barcode_ranges.iter())
            {
                let rel_dist = rel_dist_to_end(f_start as isize, read.len());
                if *passed {
                    // We only look  for the barcode in the mask region with some wiggle room in case of gaps
                    // at the flank boundaries, +/- 5 bases
                    let (bar_start, bar_end) = bar_range.unwrap();
                    let barcode_slice =
                        &read[bar_start.saturating_sub(5)..(bar_end + 5).min(read.len())];

                    // If there are barcode matches, we save all of them as such
                    let barcode_matches = self.check_barcodes(&mut searcher, barcode_slice, &flank);
                    if !barcode_matches.is_empty() {
                        for (barcode_idx, m) in barcode_matches {
                            all_matches.push(self.create_barcode_match(
                                f_start + m.start.1 as usize,
                                f_start + m.end.1 as usize,
                                m.cost as u8,
                                rel_dist,
                                query,
                                barcode_idx,
                            ));
                        }
                    // Without barcode matches it was still a "valid" flanking region, and therefore saved
                    } else {
                        all_matches.push(
                            self.create_flank_match(f_start, f_end, edit_dist, rel_dist, query),
                        );
                    }
                // If we did not have enough space to search a barcode we still had a flank
                } else {
                    all_matches
                        .push(self.create_flank_match(f_start, f_end, edit_dist, rel_dist, query));
                }
            }
        }
        self.final_assignment(&all_matches)
    }

    pub fn auto_tune_parmas(&mut self) {
        let mut searcher = self.create_searcher();
        let start_time = Instant::now();
        println!(
            "\n{}",
            "Parameter Tuning: Edit Distance Fraction"
                .bold()
                .underline()
        );
        println!("  • Range: {} - {}", "0.10".dimmed(), "0.50".dimmed());
        println!("  • Test sequences: {}\n", "10,000".dimmed());

        // Tune edit distance thresholds
        self.edit_distance_thresholds = tune_edit_distance(
            &mut searcher,
            &self.queries,
            self.n_tune_runs,
            self.fp_target,
        );

        // Pretty print each threshold
        println!("{}", "Edit distance thresholds:".green().bold());
        for (i, threshold) in self.edit_distance_thresholds.iter().enumerate() {
            println!("  • Query {:>2}: {}", i, threshold.to_string().bold());
        }

        if self.max_edits == 0 {
            self.max_edits = tune_max_edits(
                &mut searcher,
                &self.queries,
                self.n_tune_runs,
                self.fp_target,
            );
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

    fn final_assignment(&self, matches: &[BarbellMatch]) -> Vec<BarbellMatch> {
        // Sort matches based on start position
        let mut sorted_matches = matches.to_vec();
        sorted_matches.sort_by_key(|m| m.start);

        // Group overlapping matches
        let mut groups: Vec<Vec<BarbellMatch>> = Vec::new();
        let mut current_group: Vec<BarbellMatch> = Vec::new();

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
        // For barcode and flank we select lowest edit distance
        let final_matches: Vec<BarbellMatch> = groups
            .into_iter()
            .map(|group| {
                group
                    .iter()
                    .filter(|m| {
                        matches!(
                            m.label.match_type,
                            MatchType::Fbarcode | MatchType::Rbarcode
                        )
                    })
                    .min_by(|a, b| {
                        a.edit_dist
                            .unwrap()
                            .partial_cmp(&b.edit_dist.unwrap())
                            .unwrap()
                    })
                    .or_else(|| group.iter().min_by_key(|m| m.edit_dist.unwrap()))
                    .unwrap()
                    .clone()
            })
            .collect();

        final_matches
    }

    // return the index of the best match and the posterior probability
    fn check_barcodes(
        &self,
        searcher: &mut Searcher<Iupac>,
        barcode_slice: &[u8],
        mask_fillers: &FlankSeq,
    ) -> Vec<(usize, Match)> {
        let mut barcode_matches = Vec::new();
        for (q_idx, query) in mask_fillers.mask_queries.iter().enumerate() {
            let matches = searcher.search(query.as_ref(), &barcode_slice, self.max_edits as usize);
            for m in matches {
                barcode_matches.push((q_idx, m.clone()));
            }
        }
        barcode_matches
    }

    fn locate_flank(
        &self,
        searcher: &mut Searcher<Iupac>,
        flank: &FlankSeq,
        read: &StaticText,
        query_cut_off: i32,
    ) -> (
        Vec<bool>,
        Vec<(usize, usize, i32)>,
        Vec<Option<(usize, usize)>>,
    ) {
        // Use sassy search to find local minima, if present within k edits
        let matches = searcher.search(flank.seq.as_ref(), read, query_cut_off as usize);

        let mut traced_ranges = Vec::with_capacity(matches.len());
        let mut passed_mask = Vec::with_capacity(matches.len());
        let mut barcode_ranges = Vec::with_capacity(matches.len());

        // For each of the matches, we extract the mask region, if sufficient is recovered add it as passed
        // if not it still counts as possible flank, just not enough "space" to search for the barcode
        for m in matches {
            let path = m.to_path();
            let (_, mask_passed, r_range) = flank.mask_covered(&path, self.min_mask_available);
            traced_ranges.push((m.start.1 as usize, m.end.1 as usize, m.cost)); // Covered read area, with number of edits 
            passed_mask.push(mask_passed); // Whether at least self.min_mask_available is covered
            barcode_ranges.push(r_range); // Range of barcode covered, in the read
        }

        (passed_mask, traced_ranges, barcode_ranges)
    }
}

pub fn random_dna(length: usize) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let bases: Vec<u8> = vec![b'A', b'T', b'C', b'G'];
    (0..length)
        .map(|_| bases[rng.gen_range(0..bases.len())])
        .collect()
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

pub fn rel_dist_to_end(pos: isize, read_len: usize) -> isize {
    // If already negative, for a match starting before the read we can return 1
    if pos < 0 {
        return 1;
    }

    if pos <= (read_len / 2) as isize {
        if pos == 0 {
            1 // Left end
        } else {
            pos // Distance from left end (already isize)
        }
    } else if pos == read_len as isize {
        -1 // Right end
    } else {
        -(read_len as isize - pos) // Distance from right end
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
