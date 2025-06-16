use crate::pattern::pattern::*;
use crate::pattern_from_str;
use crate::types::*;
use colored::Colorize;
use indicatif::ProgressStyle;
use std::error::Error;
use std::time::Instant;

pub fn filter(
    annotated_file: &str,
    output_file: &str,
    filters: Vec<Pattern>,
) -> Result<(), Box<dyn Error>> {
    let start_time = Instant::now();

    println!("\n{}", "Filtering".bold().underline());
    println!("  • Input:  {}", annotated_file.bold());
    println!("  • Output: {}", output_file.bold());
    println!("  • Patterns: {}\n", filters.len().to_string().bold());

    // Setup progress bar
    let progress_bar = indicatif::ProgressBar::new_spinner();
    progress_bar.set_style(
        ProgressStyle::with_template("{spinner:.blue} {prefix:<12} {msg:>6} {elapsed_precise}")
            .unwrap()
            .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    progress_bar.set_prefix("Processing:");

    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(annotated_file)?;

    // Serde guided writer, with serialize and deserialize
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output_file)
        .unwrap();

    // Counters
    let mut total_reads = 0;
    let mut kept_reads = 0;

    // Process reads one group at a time
    let mut current_read_id: Option<String> = None;
    let mut current_group: Vec<BarbellMatch> = Vec::new();

    println!("Started deser");
    for result in reader.deserialize() {
        let record: BarbellMatch = result?;

        if let Some(read_id) = &current_read_id {
            if read_id != record.read.as_ref().unwrap() {
                // Process previous group
                if check_filter_pass(&mut current_group, &filters) {
                    kept_reads += 1;
                    for annotation in &current_group {
                        writer.serialize(annotation)?;
                    }
                }
                current_group.clear();
                current_read_id = Some(record.read.clone().unwrap());
                total_reads += 1;
            }
        } else {
            current_read_id = Some(record.read.clone().unwrap());
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

    println!("\n{}", "Summary".bold().underline());
    println!(
        "  • Time: {} seconds",
        start_time.elapsed().as_secs().to_string().bold()
    );
    println!("  • Total reads: {}", total_reads.to_string().bold());
    println!("  • Kept reads: {}", kept_reads.to_string().green().bold());
    println!(
        "  • Removed reads: {}\n",
        (total_reads - kept_reads).to_string().red().bold()
    );

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
        for (i, cut) in cut_positions {
            if let Some(existing_cuts) = &mut annotations[i].cuts {
                existing_cuts.push(cut);
            } else {
                annotations[i].cuts = Some(vec![cut]);
            }
        }
    }

    max_matches == annotations.len()
}
