use crate::filter::pattern::*;
use crate::pattern_from_str;
use crate::progress::{ProgressTracker, print_header, print_summary_stats};
use crate::search::barcodes::BarcodeType;
use crate::search::searcher::BarbellMatch;
use colored::Colorize;
use indicatif::ProgressStyle;
use std::error::Error;
use std::time::Instant;

pub fn filter(
    annotated_file: &str,
    output_file: &str,
    filters: Vec<Pattern>,
) -> Result<(), Box<dyn Error>> {
    let mut progress = ProgressTracker::new();
    print_header("Pattern Filtering");

    progress.step("Configuration");
    progress.indent();
    progress.substep(&format!("Input: {}", annotated_file));
    progress.substep(&format!("Output: {}", output_file));
    progress.substep(&format!("Patterns: {}", filters.len()));
    progress.dedent();

    // Setup progress bar
    let progress_bar = indicatif::ProgressBar::new_spinner();
    progress_bar.set_style(
        ProgressStyle::with_template("{spinner:.blue} {prefix:<12} {msg:>6} {elapsed_precise}")
            .unwrap()
            .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    progress_bar.set_prefix("Processing:");

    progress.step("Reading annotations");
    progress.indent();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(annotated_file)?;

    // Serde guided writer, with serialize and deserialize
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output_file)
        .unwrap();
    progress.success("Files initialized");
    progress.dedent();

    // Counters
    let mut total_reads = 0;
    let mut kept_reads = 0;

    // Process reads one group at a time
    let mut current_read_id: Option<String> = None;
    let mut current_group: Vec<BarbellMatch> = Vec::new();

    progress.step("Processing read groups");
    progress.indent();
    for result in reader.deserialize() {
        let record: BarbellMatch = result?;

        if let Some(read_id) = &current_read_id {
            if read_id != &record.read_id {
                // Process previous group
                if check_filter_pass(&mut current_group, &filters) {
                    kept_reads += 1;
                    for annotation in &current_group {
                        writer.serialize(annotation)?;
                    }
                }
                current_group.clear();
                current_read_id = Some(record.read_id.clone());
                total_reads += 1;
            }
        } else {
            current_read_id = Some(record.read_id.clone());
            total_reads += 1;
        }

        current_group.push(record);
        progress_bar.set_message(format!("{}", total_reads));
    }

    // Process the last group
    if !current_group.is_empty() && check_filter_pass(&mut current_group, &filters) {
        kept_reads += 1;
        for annotation in &current_group {
            writer.serialize(annotation)?;
        }
    }

    writer.flush()?;
    progress_bar.finish_with_message("Done!");
    progress.dedent();

    // Print summary
    print_summary_stats(
        total_reads,
        kept_reads,
        kept_reads, // For filtering, mapped = kept
        kept_reads, // For filtering, trimmed = kept
        progress.elapsed(),
    );

    progress.success("Filtering completed successfully");
    progress.print_elapsed();

    Ok(())
}

pub fn filter_from_pattern_str(
    annotated_file: &str,
    pattern_str: &str,
    output_file: &str,
) -> Result<(), Box<dyn Error>> {
    // use pattern macro to convert pattern_str to pattern
    let pattern = pattern_from_str!(pattern_str);
    filter(annotated_file, output_file, vec![pattern])
}

pub fn filter_from_text_file(
    annotated_file: &str,
    text_file: &str,
    output_file: &str,
) -> Result<(), Box<dyn Error>> {
    // read the text file into a vector of strings
    let patterns = std::fs::read_to_string(text_file)?;
    let patterns = patterns
        .split("\n")
        .map(|s| s.trim())
        .collect::<Vec<&str>>();
    let patterns = patterns
        .iter()
        .map(|s| pattern_from_str!(s))
        .collect::<Vec<Pattern>>();
    filter(annotated_file, output_file, patterns)
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
    if max_matches > 0 && best_cut_positions.is_some() {
        let cut_positions = best_cut_positions.unwrap();
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
