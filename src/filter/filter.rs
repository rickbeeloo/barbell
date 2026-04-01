use crate::annotate::searcher::BarbellMatch;
use crate::config::FilterConfig;
use crate::filter::pattern::*;
use crate::pattern_from_str;
use crate::progress::progress::{FILTER_PROGRESS_SPECS, ProgressTracker};
use std::error::Error;
use std::fs::File;
use std::path::Path;

pub fn filter(
    annotated_file: &str,
    output_file: &str,
    dropped_out_file: Option<&str>,
    filters: &[Pattern],
    config: &FilterConfig,
) -> Result<(), Box<dyn Error>> {
    // Setup progress tracking using the shared progress module.
    let progress = if config.verbose {
        let log_dir = Path::new(output_file)
            .parent()
            .unwrap_or_else(|| Path::new("."));
        ProgressTracker::new_with_logging(&FILTER_PROGRESS_SPECS, "filter", log_dir)
    } else {
        ProgressTracker::new(&FILTER_PROGRESS_SPECS)
    };

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
                progress.inc(0);
                // Process previous group
                if check_filter_pass(&mut current_group, filters) {
                    progress.inc(1);
                    for annotation in &current_group {
                        writer.serialize(annotation)?;
                    }
                } else {
                    progress.inc(2);
                    if let Some(dropped_writer) = &mut dropped_writer {
                        for annotation in &current_group {
                            dropped_writer.serialize(annotation)?;
                        }
                    }
                }
                current_group.clear();
                current_read_id = Some(record.read_id.clone());

                // Update messages with current counts
                progress.refresh();
            }
        } else {
            current_read_id = Some(record.read_id.clone());
        }

        current_group.push(record);
    }

    // Process the last group
    if !current_group.is_empty() {
        progress.inc(0);
        if check_filter_pass(&mut current_group, filters) {
            progress.inc(1);
            for annotation in &current_group {
                writer.serialize(annotation)?;
            }
        } else {
            progress.inc(2);
            if let Some(dropped_writer) = &mut dropped_writer {
                for annotation in &current_group {
                    dropped_writer.serialize(annotation)?;
                }
            }
        }

        // Final live update before finishing
        progress.refresh();
    }

    writer.flush()?;

    // Flush dropped writer if it exists
    if let Some(mut dropped_writer) = dropped_writer {
        dropped_writer.flush()?;
    }

    // Finish the progress bars with final counts
    progress.finish("reads");

    Ok(())
}

pub fn filter_from_pattern_str(
    annotated_file: &str,
    pattern_str: &str,
    output_file: &str,
    dropped_out_file: Option<&str>,
    config: &FilterConfig,
) -> Result<(), Box<dyn Error>> {
    // use pattern macro to convert pattern_str to pattern
    let pattern = pattern_from_str!(pattern_str);
    filter(annotated_file, output_file, dropped_out_file, &[pattern], config)
}

pub fn filter_from_text_file(
    annotated_file: &str,
    text_file: &str,
    output_file: &str,
    dropped_out_file: Option<&str>,
    config: &FilterConfig,
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
    filter(annotated_file, output_file, dropped_out_file, &patterns, config)
}

fn check_filter_pass(annotations: &mut [BarbellMatch], patterns: &[Pattern]) -> bool {
    // Track both the maximum number of matches and the cut positions
    let mut max_matches = 0;
    let mut best_cut_positions: Option<Vec<(usize, Cut)>> = None;

    for pattern in patterns {
        // cut position is (match_idx, cut)
        let (is_match, cut_positions) = match_pattern(annotations, pattern);
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
