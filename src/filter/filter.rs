use crate::annotate::barcodes::BarcodeType;
use crate::annotate::searcher::BarbellMatch;
use crate::filter::filter_strategy::FilterStrategy;
use crate::filter::pattern::*;
use crate::pattern_from_str;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::time::Duration;

fn create_progress_bar() -> (ProgressBar, ProgressBar, ProgressBar) {
    // Create multiprogress bar
    let multi_progress = MultiProgress::new();
    let total_bar = multi_progress.add(ProgressBar::new_spinner());
    let kept_bar = multi_progress.add(ProgressBar::new_spinner());
    let dropped_bar = multi_progress.add(ProgressBar::new_spinner());

    // Style the progress bars with colors - more compact layout
    total_bar.set_style(
        ProgressStyle::with_template(
            "{spinner:.cyan} {prefix:.bold.white:<8} {msg:.bold.cyan:>6} {elapsed:.dim}",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    kept_bar.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} {prefix:.bold.white:<8} {msg:.bold.green:>6} {elapsed:.dim}",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    dropped_bar.set_style(
        ProgressStyle::with_template(
            "{spinner:.red} {prefix:.bold.white:<8} {msg:.bold.red:>6} {elapsed:.dim}",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );

    total_bar.enable_steady_tick(Duration::from_millis(100));
    kept_bar.enable_steady_tick(Duration::from_millis(120)); // Slightly different timing for visual variety
    dropped_bar.enable_steady_tick(Duration::from_millis(140));

    total_bar.set_prefix("Total:");
    kept_bar.set_prefix("Kept:");
    dropped_bar.set_prefix("Dropped:");

    (total_bar, kept_bar, dropped_bar)
}

pub fn filter(
    annotated_file: &str,
    output_file: &str,
    dropped_out_file: Option<&str>,
    filters: &[Pattern],
    filter_strategy: &FilterStrategy,
) -> Result<(), Box<dyn Error>> {
    // Setup progress bar
    let (total_bar, kept_bar, dropped_bar) = create_progress_bar();

    println!("{}", filter_strategy.display());

    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(annotated_file)?;

    // Serde guided writer, with serialize and deserialize
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output_file)
        .unwrap();

    // In case dropped reads should also be written to another file open
    // dropped writer if provided
    let mut dropped_writer: Option<csv::Writer<File>> = None;
    if let Some(dropped_out_file) = dropped_out_file {
        dropped_writer = Some(
            csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_path(dropped_out_file)
                .unwrap(),
        );
    }

    // Process reads one group at a time
    let mut current_read_id: Option<String> = None;
    let mut current_group: Vec<BarbellMatch> = Vec::new();

    for result in reader.deserialize() {
        let record: BarbellMatch = result?;

        if let Some(read_id) = &current_read_id {
            if read_id != &record.read_id {
                total_bar.inc(1);
                // Process previous group
                let passed = filter_annotations(&mut current_group, filters, filter_strategy);
                if passed {
                    kept_bar.inc(1);
                    for annotation in &current_group {
                        writer.serialize(annotation)?;
                    }
                } else {
                    dropped_bar.inc(1);
                    if let Some(dropped_writer) = &mut dropped_writer {
                        for annotation in &current_group {
                            dropped_writer.serialize(annotation)?;
                        }
                    }
                }
                current_group.clear();
                current_read_id = Some(record.read_id.clone());

                // Update messages with current counts
                total_bar.set_message(total_bar.position().to_string());
                kept_bar.set_message(kept_bar.position().to_string());
                dropped_bar.set_message(dropped_bar.position().to_string());
            }
        } else {
            current_read_id = Some(record.read_id.clone());
        }

        current_group.push(record);
    }

    // Process the last group
    if !current_group.is_empty() {
        total_bar.inc(1);
        let passed = filter_annotations(&mut current_group, filters, filter_strategy);
        if passed {
            kept_bar.inc(1);
            for annotation in &current_group {
                writer.serialize(annotation)?;
            }
        } else {
            dropped_bar.inc(1);
            if let Some(dropped_writer) = &mut dropped_writer {
                for annotation in &current_group {
                    dropped_writer.serialize(annotation)?;
                }
            }
        }

        // Final live update before finishing
        total_bar.set_message(total_bar.position().to_string());
        kept_bar.set_message(kept_bar.position().to_string());
        dropped_bar.set_message(dropped_bar.position().to_string());
    }

    writer.flush()?;

    // Flush dropped writer if it exists
    if let Some(mut dropped_writer) = dropped_writer {
        dropped_writer.flush()?;
    }

    // Finish the progress bars with final counts
    let total_count = total_bar.position();
    let kept_count = kept_bar.position();
    let dropped_count = dropped_bar.position();

    total_bar.finish_with_message(format!("{total_count} reads"));
    kept_bar.finish_with_message(format!("{kept_count} reads"));
    dropped_bar.finish_with_message(format!("{dropped_count} reads"));

    Ok(())
}

pub fn filter_from_pattern_str(
    annotated_file: &str,
    pattern_str: &str,
    output_file: &str,
    dropped_out_file: Option<&str>,
    filter_strategy: &FilterStrategy,
) -> Result<(), Box<dyn Error>> {
    // use pattern macro to convert pattern_str to pattern
    let pattern = pattern_from_str!(pattern_str);
    filter(
        annotated_file,
        output_file,
        dropped_out_file,
        &[pattern],
        &filter_strategy,
    )
}

pub fn filter_from_text_file(
    annotated_file: &str,
    text_file: &str,
    output_file: &str,
    dropped_out_file: Option<&str>,
    filter_strategy: &FilterStrategy,
) -> Result<(), Box<dyn Error>> {
    // read the text file into a vector of strings
    let patterns = std::fs::read_to_string(text_file)?;

    let patterns = patterns
        .split('\n')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect::<Vec<&str>>();

    let patterns = patterns
        .iter()
        .map(|s| pattern_from_str!(s))
        .collect::<Vec<Pattern>>();
    filter(
        annotated_file,
        output_file,
        dropped_out_file,
        &patterns,
        &filter_strategy,
    )
}

#[derive(Debug)]
struct Interval {
    start: usize,
    end: usize,
    end_inclusive: bool, // true for end], false for end)
    cut_idx: usize,
}

fn greedy_solve(all_cuts: Vec<Vec<(usize, Cut)>>) -> Vec<Vec<(usize, Cut)>> {
    let mut intervals = Vec::new();

    for (i, cut) in all_cuts.iter().enumerate() {
        match cut.len() {
            1 => {
                let pos = cut[0].0;
                match cut[0].1.direction {
                    CutDirection::After => {
                        // [pos, ∞) - start is inclusive, end is open (at infinity)
                        intervals.push(Interval {
                            start: pos,
                            end: usize::MAX,
                            end_inclusive: false,
                            cut_idx: i,
                        });
                    }
                    CutDirection::Before => {
                        // (-∞, pos] - start is open (at -infinity), end is inclusive
                        intervals.push(Interval {
                            start: 0,
                            end: pos,
                            end_inclusive: true,
                            cut_idx: i,
                        });
                    }
                }
            }
            2 => {
                let (pos_a, dir_a) = (cut[0].0, &cut[0].1.direction);
                let (pos_b, dir_b) = (cut[1].0, &cut[1].1.direction);

                // Determine the interval bounds
                let (start, end) = (pos_a.min(pos_b), pos_a.max(pos_b));

                // For a closed interval [a, b], we need:
                // - After at the lower position (cuts after a, so includes a)
                // - Before at the upper position (cuts before b, so includes b)
                let end_inclusive = if pos_a > pos_b {
                    matches!(dir_a, CutDirection::Before)
                } else {
                    matches!(dir_b, CutDirection::Before)
                };

                intervals.push(Interval {
                    start,
                    end,
                    end_inclusive,
                    cut_idx: i,
                });
            }
            _ => {
                panic!("Warning: unexpected cut length {}", cut.len());
            }
        }
    }

    // Sort by end coordinate, then by end_inclusive (exclusive ends first to prefer earlier cuts)
    intervals.sort_by_key(|interval| (interval.end, interval.end_inclusive));

    let mut result = Vec::new();
    let mut last_interval: Option<Interval> = None;

    for interval in intervals {
        let overlaps = if let Some(ref last) = last_interval {
            // Check if intervals overlap
            // They overlap only if new_start is strictly before last_end
            // Touching at a single boundary point (start == end) is NOT an overlap
            interval.start < last.end
        } else {
            false
        };

        if !overlaps {
            result.push(all_cuts[interval.cut_idx].clone());
            last_interval = Some(interval);
        }
    }

    result
}

fn no_internal(annotations: &mut [BarbellMatch]) -> bool {
    for annotation in annotations.iter() {
        let end_diff = annotation.read_len.saturating_sub(250);
        if annotation.read_start_flank > 250 && annotation.read_end_flank < end_diff {
            return false;
        }
    }
    true
}

fn check_diff_barcodes(annotations: &mut [BarbellMatch]) -> bool {
    let mut unique = HashSet::new();
    for a in annotations.iter() {
        if a.match_type == BarcodeType::Ftag
            || a.match_type == BarcodeType::Rtag
            || a.match_type == BarcodeType::Fbar
            || a.match_type == BarcodeType::Rbar
        {
            unique.insert(&a.label);
            if unique.len() > 1 {
                return false;
            }
        }
    }
    true
}

fn check_if_clean(annotations: &mut [BarbellMatch]) -> bool {
    let read_len = annotations[0].read_len;

    let mut cov = 0;
    for a in annotations.iter() {
        let start = a.read_start_flank;
        let end = a.read_end_flank;
        cov += end - start;
    }

    let cov_perc = cov as f32 / read_len as f32;

    if cov_perc < 0.5 {
        return true;
    }

    false
}

fn filter_annotations(
    annotations: &mut [BarbellMatch],
    patterns: &[Pattern],
    strategy: &FilterStrategy,
) -> bool {
    // In case exact is necessary, we just a custom filtering
    // where the read *only* passes when the matches exactly match the pattern
    if strategy.exact_enabled() {
        return check_filter_pass_full_pattern(annotations, patterns);
    }

    // In case not exact we first use the prefilters
    let mut passed = true;

    if strategy.terminal_enabled() {
        passed &= no_internal(annotations);
    }

    if strategy.unique_labels_enabled() {
        passed &= check_diff_barcodes(annotations);
    }

    if strategy.clean_enabled() {
        passed &= check_if_clean(annotations);
    }

    if strategy.flexible_enabled() {
        // No extra check here
    }

    if !passed {
        return false;
    }

    // If the prefilter passed, we can now match the patterns
    check_filter_pass_sub_pattern(annotations, patterns)
}

fn check_filter_pass_full_pattern(annotations: &mut [BarbellMatch], patterns: &[Pattern]) -> bool {
    // Track both the maximum number of matches and the cut positions
    let mut max_matches = 0;
    let mut best_cut_positions: Option<Vec<(usize, Cut)>> = None;

    for pattern in patterns {
        // cut position is (match_idx, cut)
        let (is_match, cut_positions) = match_full_pattern(annotations, pattern);
        if is_match {
            let pattern_len = pattern.elements.len();
            if pattern_len > max_matches {
                max_matches = pattern_len;
                best_cut_positions = Some(cut_positions);
            }
        }
    }

    // If we have a match and cut positions, update all annotations in the group
    if max_matches > 0
        && let Some(cut_positions) = best_cut_positions
    {
        for (cut_match_idx, cut) in cut_positions {
            if let Some(existing_cuts) = &mut annotations[cut_match_idx].cuts {
                existing_cuts.push((cut, cut_match_idx));
            } else {
                annotations[cut_match_idx].cuts = Some(vec![(cut, cut_match_idx)]);
            }
        }
    }

    max_matches == annotations.len()
}

fn check_filter_pass_sub_pattern(annotations: &mut [BarbellMatch], patterns: &[Pattern]) -> bool {
    // We just keep track of every sub pattern we can find for all patterns
    let mut all_cuts = vec![];

    for pattern in patterns.iter() {
        // cut position is (match_idx, cut)
        let (is_match, cut_positions) = match_sub_pattern(annotations, pattern);

        if is_match {
            all_cuts.extend(cut_positions);
        }
    }

    if all_cuts.is_empty() {
        return false;
    }

    // After collecting all cuts, we have to collapse them to single coverage
    // for now, use greedy solver
    let mut best_cut_positions = greedy_solve(all_cuts);

    // Since in sub pattern matching we allow multiple matches our cut ids are not meaninful
    // we make sure each group has a unique cut id
    // todo: check input such that each pattern has ONLY at most two cut identifiers (>>, <<)
    // if not, panic
    for (id_accum, cut_group) in best_cut_positions.iter_mut().enumerate() {
        for (_, cut) in cut_group.iter_mut() {
            cut.group_id = id_accum;
        }
    }

    // Per group that we found, we have to update the cuts for all annotations WITHIN
    // that group
    for cut_group in best_cut_positions.iter() {
        for (cut_match_idx, cut) in cut_group {
            if let Some(existing_cuts) = &mut annotations[*cut_match_idx].cuts {
                existing_cuts.push((cut.clone(), *cut_match_idx));
            } else {
                annotations[*cut_match_idx].cuts = Some(vec![(cut.clone(), *cut_match_idx)]);
            }
        }
    }

    !best_cut_positions.is_empty()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotate::barcodes::BarcodeType;
    use crate::annotate::searcher::BarbellMatch;
    use crate::filter::filter_strategy::FilterMode;
    use sassy::Strand;

    #[test]
    fn test_greedy_solve() {
        // Say the user supplied two patterns
        let pattern1 = pattern_from_str!("Ftag[fw, *, @left(0..250), >>]");
        let pattern2 =
            pattern_from_str!("Ftag[fw, *, @left(0..250), >>]__Rtag[<<, fw, *, @right(0..250)]");

        let mut matches = vec![
            BarbellMatch::new(
                0,   // read_start_bar
                100, // read_end_bar
                0,   // read_start_flank
                100, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                250, // read_len
                "test".to_string(),
                0,
                None,
            ),
            BarbellMatch::new(
                100, // read_start_bar
                200, // read_end_bar
                100, // read_start_flank
                200, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Rtag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                250, // read_len
                "test".to_string(),
                0,
                None,
            ),
        ];

        // This mutates the matches with cut information
        let mut matches2 = matches.clone(); // For second test 
        let passed = check_filter_pass_sub_pattern(&mut matches, &[pattern1, pattern2]);
        assert!(passed);
        let match1 = matches[0].clone();
        let match2 = matches[1].clone();
        // If indeed greedy, second pattern Ftag__Rtag should be used
        // meaning first match should have Before label, and second should have After label
        assert_eq!(
            match1.cuts,
            Some(vec![(
                Cut {
                    group_id: 0,
                    direction: CutDirection::After
                },
                0
            )])
        );
        assert_eq!(
            match2.cuts,
            Some(vec![(
                Cut {
                    group_id: 0,
                    direction: CutDirection::Before
                },
                1
            )])
        );

        // Say we did the same using two patterns seperate, it should still pass but being split
        let pattern1 = pattern_from_str!("Ftag[fw, *, @left(0..250), >>]");
        let pattern2 = pattern_from_str!("Rtag[<<, fw, *, @right(0..250)]");
        let pattern3 = pattern_from_str!("Rtag[fw, *, @right(0..250), >>]");

        // This is a complex case, first of all both the Ftag and Rtag
        // effectively cover the same region, just from two different directions
        // we should just pick one (arbitrary). The second complicated thing is
        // that the second Rtag is the same Rtag in this case as the first one
        // just that we have to cut it to the right now
        /*
           Ftag---------Rtag-------
           [>>>>>>>>>>>>]
           [<<<<<<<<<<<<] (label based on Ftag,Rtag)
                        [>>>>>>>>>>>] (label based on Rtag)
                        ^ Note, Rtag should have two cuts in this case one Before, and one After
        */
        let first_expected_cuts = Some(vec![(
            Cut {
                group_id: 0,
                direction: CutDirection::After,
            },
            0,
        )]);

        let second_expected_cuts = Some(vec![
            (
                Cut {
                    group_id: 0,
                    direction: CutDirection::Before,
                },
                1,
            ),
            (
                Cut {
                    group_id: 1, // In this case indicating second pattern cut
                    direction: CutDirection::After,
                },
                1,
            ),
        ]);

        let passed = check_filter_pass_sub_pattern(&mut matches2, &[pattern1, pattern2, pattern3]);
        assert!(passed);
        assert_eq!(matches2.len(), 2);
        assert_eq!(matches2[0].cuts, first_expected_cuts);
        assert_eq!(matches2[1].cuts, second_expected_cuts);
    }

    #[test]
    fn test_sub_pattern_different_barcodes() {
        // Say the user supplied two patterns
        let pattern1 = pattern_from_str!("Ftag[fw, *, @left(0..250), >>]");
        let pattern2 =
            pattern_from_str!("Ftag[fw, *, @left(0..250), >>]__Rtag[<<, fw, *, @right(0..250)]");

        let matches = vec![
            BarbellMatch::new(
                0,   // read_start_bar
                100, // read_end_bar
                0,   // read_start_flank
                100, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                250, // read_len
                "test".to_string(),
                0,
                None,
            ),
            BarbellMatch::new(
                100, // read_start_bar
                200, // read_end_bar
                100, // read_start_flank
                200, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Rtag,
                0, // flank_cost
                0, // barcode_cost
                "YYY".to_string(),
                Strand::Fwd,
                250, // read_len
                "test".to_string(),
                0,
                None,
            ),
        ];

        // This mutates the matches with cut information

        let mut matches_test_1 = matches.clone();
        let passed = filter_annotations(
            &mut matches_test_1,
            &[pattern1.clone(), pattern2.clone()],
            &FilterStrategy::from_modes(&[FilterMode::Exact]),
        );
        // Without any sub patttern strategy all of them should pass
        assert!(passed);

        // If we change to Flexible it should still pass
        let mut matches_test_2 = matches.clone();
        let passed = filter_annotations(
            &mut matches_test_2,
            &[pattern1.clone(), pattern2.clone()],
            &FilterStrategy::from_modes(&[FilterMode::Flexible]),
        );
        assert!(passed);

        // If we change to terminal it should still pass
        let mut matches_test_3 = matches.clone();
        let passed = filter_annotations(
            &mut matches_test_3,
            &[pattern1.clone(), pattern2.clone()],
            &FilterStrategy::from_modes(&[FilterMode::Terminal]),
        );
        assert!(passed);

        // However if we changed the read len to 1K or so
        // the terminal filter should fail
        let mut matches_test_4 = matches.clone();
        matches_test_4.iter_mut().for_each(|m| {
            m.read_len = 1000;
            m.read_start_flank = 251;
        });

        let passed = filter_annotations(
            &mut matches_test_4,
            &[pattern1.clone(), pattern2.clone()],
            &FilterStrategy::from_modes(&[FilterMode::Terminal]),
        );
        assert!(!passed);

        // If we use both filters (complete) it should also fail
        let mut matches_test_5 = matches.clone();
        matches_test_5.iter_mut().for_each(|m| m.read_len = 1000);
        let passed = filter_annotations(
            &mut matches_test_5,
            &[pattern1.clone(), pattern2.clone()],
            &FilterStrategy::from_modes(&[
                FilterMode::Exact,
                FilterMode::Clean,
                FilterMode::Terminal,
            ]),
        );
        assert!(!passed);
    }
}
