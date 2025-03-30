use crate::barbell::pattern_assign::*;

use crate::pattern_from_str;

use std::collections::HashMap;
use serde::Deserialize;
use std::error::Error;
use std::time::Instant;
use indicatif::{ProgressBar, ProgressStyle};
use colored::Colorize;

#[derive(Debug, Deserialize, Clone)]
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
    pub record_set_idx: usize,
    pub record_idx: usize,
    #[serde(default)]
    #[serde(deserialize_with = "deserialize_cuts")]
    pub cuts: Option<Vec<Cut>>,
}

// Update the deserializer
fn deserialize_cuts<'de, D>(deserializer: D) -> Result<Option<Vec<Cut>>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let s: String = String::deserialize(deserializer)?;
    
    if s.is_empty() || s == "-" {
        return Ok(None);
    }
    
    let cuts = if s.contains(',') {
        s.split(',')
            .map(str::trim)
            .flat_map(Cut::from_string)
            .collect::<Vec<_>>()
    } else {
        Cut::from_string(&s).into_iter().collect()
    };
    
    Ok(Some(cuts))
}

impl AnnotationLine {
    pub fn to_record(&self) -> csv::StringRecord {
        let cut_str = match &self.cuts {
            Some(cuts) if !cuts.is_empty() => cuts.iter()
                .map(|cut| cut.to_string())
                .collect::<Vec<_>>()
                .join(","),
            _ => "-".to_string()
        };
            
        let start_str = self.start.to_string();
        let end_str = self.end.to_string();
        let edits_str = self.edits.to_string();
        let dist_str = self.dist_to_end.to_string();
        let len_str = self.read_len.to_string();
        let set_idx_str = self.record_set_idx.to_string();
        let idx_str = self.record_idx.to_string();

        let fields = vec![
            self.read.as_str(),
            self.label.as_str(),
            &start_str,
            &end_str,
            &edits_str,
            &dist_str,
            &len_str,
            &set_idx_str,
            &idx_str,
            &cut_str,
        ];
        
        csv::StringRecord::from(fields)
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

    // Get the header before starting deserialization
    let mut headers = reader.headers()?.clone();
    
    // Always add cut_positions column if it doesn't exist
    if !headers.iter().any(|h| h == "cuts") {
        headers.push_field("cuts");
    }

    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output_file)?;

    // Write the headers to the output file
    writer.write_record(&headers)?;
    
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
                if check_filter_pass(&mut current_group, &filters) {
                    kept_reads += 1;
                    for annotation in &current_group {
                        writer.write_record(&annotation.to_record())?;
                    }
                }
                current_group.clear();
                current_read_id = Some(record.read.clone());
                total_reads += 1;
            }
        } else {
            current_read_id = Some(record.read.clone());
            total_reads += 1;
        }
        
        current_group.push(record);
        progress_bar.set_message(format!("{}", total_reads));
    }
    
    // Process the last group
    if !current_group.is_empty() && check_filter_pass(&mut current_group, &filters) {
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

fn check_filter_pass(annotations: &mut [AnnotationLine], patterns: &[Pattern]) -> bool {
    // Convert annotations to matches
    let matches: Vec<Match> = annotations.iter().map(
        |annotation| {
            let match_str = EncodedMatchStr::unstringify(&annotation.label);
            Match::new(match_str, annotation.start, annotation.end, annotation.edits, annotation.dist_to_end)
        }
    ).collect();

    let read_len = annotations[0].read_len;

    // Track both the maximum number of matches and the cut positions
    let mut max_matches = 0;
    let mut best_cut_positions: Option<Vec<Cut>> = None;

    for pattern in patterns {
        let (is_match, cut_positions) = match_pattern(&matches, pattern, read_len);
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
        for annotation in annotations.iter_mut() {
            annotation.cuts = Some(best_cut_positions.clone().unwrap());
        }
    }

    max_matches == annotations.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use crate::barbell::pattern_assign::{Cut, CutDirection};

    #[test]
    fn test_deserialize_cuts_string_format() {
        // Test single cut with full line
        let json_value = json!({ 
            "read": "dac457ea-df0b-4c27-a58b-a52a411dc1ee",
            "label": "barcode35#fw#Fbar",
            "start": 0,
            "end": 107,
            "edits": 19,
            "dist.to.end": 1,
            "read.len": 2374,
            "record_set_idx": 0,
            "record_idx": 0,
            "cuts": "After(118)"
        });
        let annotation: AnnotationLine = serde_json::from_value(json_value).unwrap();
        println!("Annotation: {:?}", annotation);
        assert_eq!(annotation.cuts.unwrap()[0], Cut::new(118, CutDirection::After));

        // Test multiple cuts with full line
        let json_value = json!({ 
            "read": "dac457ea-df0b-4c27-a58b-a52a411dc1ee",
            "label": "barcode35#fw#Fbar",
            "start": 0,
            "end": 107,
            "edits": 19,
            "dist.to.end": 1,
            "read.len": 2374,
            "record_set_idx": 0,
            "record_idx": 0,
            "cuts": "After(118),Before(50)"
        });
        let annotation: AnnotationLine = serde_json::from_value(json_value).unwrap();
        let cuts = annotation.cuts.unwrap();
        assert_eq!(cuts[0], Cut::new(118, CutDirection::After));
        assert_eq!(cuts[1], Cut::new(50, CutDirection::Before));

        // Test empty cut with full line
        let json_value = json!({ 
            "read": "dac457ea-df0b-4c27-a58b-a52a411dc1ee",
            "label": "barcode35#fw#Fbar",
            "start": 0,
            "end": 107,
            "edits": 19,
            "dist.to.end": 1,
            "read.len": 2374,
            "record_set_idx": 0,
            "record_idx": 0,
            "cuts": "-"
        });
        let annotation: AnnotationLine = serde_json::from_value(json_value).unwrap();
        assert_eq!(annotation.cuts, None);
    }
}