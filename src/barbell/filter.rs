use crate::barbell::pattern_assign::*;

use crate::pattern_from_str;

use std::collections::HashMap;
use serde::Deserialize;
use std::error::Error;
use std::time::Instant;
use indicatif::{ProgressBar, ProgressStyle};
use colored::Colorize;

#[derive(Debug, Deserialize)]
pub struct AnnotationLine {
    pub read: String,
    pub label: String,
    pub start: usize,
    pub end: usize,
    pub edits: i32,
    #[serde(rename = "dist.to.end")]
    pub dist_to_end: isize,
    #[serde(rename = "read.len")]
    pub read_len: usize,
}

impl AnnotationLine {
    pub fn to_record(&self) -> csv::StringRecord {
        csv::StringRecord::from(vec![
            self.read.as_str(),
            self.label.as_str(),
            &self.start.to_string(),
            &self.end.to_string(),
            &self.edits.to_string(),
            &self.dist_to_end.to_string(),
            &self.read_len.to_string(),
        ])
    }
}

pub fn filter(annotated_file: &str, output_file: &str, filters: Vec<Pattern>) -> Result<(), Box<dyn Error>> {
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
            .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
    );
    progress_bar.set_prefix("Processing:");

    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(annotated_file)?;

    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output_file)?;
    
    // Counters
    let mut total_reads = 0;
    let mut kept_reads = 0;
    
    // Process reads one group at a time
    let mut current_read_id: Option<String> = None;
    let mut current_group: Vec<AnnotationLine> = Vec::new();

    for result in reader.deserialize() {
        let record: AnnotationLine = result?;
        
        if let Some(read_id) = &current_read_id {
            if *read_id != record.read {
                // Process previous group
                if check_filter_pass(&current_group, &filters) {
                    kept_reads += 1;
                    for annotation in &current_group {
                        writer.write_record(&annotation.to_record())?;
                    }
                }
                current_group.clear();
                current_read_id = Some(record.read.clone());
            }
        } else {
            current_read_id = Some(record.read.clone());
        }
        
        current_group.push(record);
        total_reads += 1;
        progress_bar.set_message(format!("{}", total_reads));
    }
    
    // Process the last group
    if !current_group.is_empty() && check_filter_pass(&current_group, &filters) {
        kept_reads += 1;
        for annotation in &current_group {
            writer.write_record(&annotation.to_record())?;
        }
    }
    
    writer.flush()?;
    progress_bar.finish_with_message("Done!");

    println!("\n{}", "Summary".bold().underline());
    println!("  • Time: {} seconds", start_time.elapsed().as_secs().to_string().bold());
    println!("  • Total reads: {}", total_reads.to_string().bold());
    println!("  • Kept reads: {}", kept_reads.to_string().green().bold());
    println!("  • Removed reads: {}\n", (total_reads - kept_reads).to_string().red().bold());

    Ok(())
}

pub fn filter_from_pattern_str(annotated_file: &str, pattern_str: &str, output_file: &str) -> Result<(), Box<dyn Error>> {
    // use pattern macro to convert pattern_str to pattern
    let pattern = pattern_from_str!(pattern_str);
    filter(annotated_file, output_file, vec![pattern])
}

pub fn filter_from_text_file(annotated_file: &str, text_file: &str, output_file: &str) -> Result<(), Box<dyn Error>> {
    // read the text file into a vector of strings
    let patterns = std::fs::read_to_string(text_file)?;
    let patterns = patterns.split("\n").map(|s| s.trim()).collect::<Vec<&str>>();
    let patterns = patterns.iter().map(|s| pattern_from_str!(s)).collect::<Vec<Pattern>>();
    filter(annotated_file, output_file, patterns)
}

fn check_filter_pass(annotations: &[AnnotationLine], patterns: &[Pattern]) -> bool {


    // Check if any of the annotations mamtches any filter
    // first we have to convert the label filed to encoded match string
    let matches: Vec<Match> = annotations.iter().map(
        |annotation| {
            let match_str = EncodedMatchStr::unstringify(&annotation.label);
            Match::new(match_str, annotation.start, annotation.end, annotation.edits, annotation.dist_to_end)
        }
    ).collect();


    // Read length should be the same for all 
    let read_len = annotations[0].read_len;

    // now we have to check if any of the encoded match strings match any of the filters
    let mut max_matches = 0;
    for pattern in patterns {
        if match_pattern(&matches, pattern, read_len) {

            let pattern_len = pattern.elements.len();
            if pattern_len > max_matches {
                max_matches = pattern_len;
            }
        }
    }

    // Max matches should not be zero, 
    // Total number of matches 
    max_matches == annotations.len()
}