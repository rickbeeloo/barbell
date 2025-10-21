use crate::annotate::searcher::BarbellMatch;
use crate::filter::pattern::*;
use crate::pattern_from_str;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
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
) -> Result<(), Box<dyn Error>> {
    // Setup progress bar
    let (total_bar, kept_bar, dropped_bar) = create_progress_bar();

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
                if check_filter_pass(&mut current_group, filters) {
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
        if check_filter_pass(&mut current_group, filters) {
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
) -> Result<(), Box<dyn Error>> {
    // use pattern macro to convert pattern_str to pattern
    let pattern = pattern_from_str!(pattern_str);
    filter(annotated_file, output_file, dropped_out_file, &[pattern])
}

pub fn filter_from_text_file(
    annotated_file: &str,
    text_file: &str,
    output_file: &str,
    dropped_out_file: Option<&str>,
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
    filter(annotated_file, output_file, dropped_out_file, &patterns)
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
