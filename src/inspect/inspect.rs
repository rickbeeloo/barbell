use crate::types::*;
use colored::Colorize;
use std::collections::HashMap;
use std::error::Error;

const BUCKET_SIZE: usize = 250;

pub fn get_group_structure(group: &[BarbellMatch]) -> String {
    if group.is_empty() {
        return String::new();
    }

    // Helper function to bucket positions into ranges
    // For example: 0-50 -> (0 to 50), 51-100 -> (50 to 100), etc.
    fn bucket_position(pos: usize) -> usize {
        // We can adjust this
        (pos / BUCKET_SIZE) * BUCKET_SIZE
    }

    // Hashmap to keep track of all the labels we see to later replace in pattern if applicable
    // just count with labels
    let mut label_map: HashMap<String, usize> = HashMap::new();

    let pattern_elements: Vec<String> = group
        .iter()
        .map(|annotation| {
            // Bucket the positions into ranges
            let start_bucket = bucket_position(annotation.start);
            let end_bucket = bucket_position(annotation.end);

            // Determine if it's closer to left or right end
            let position = if annotation.rel_dist_to_end > 0 {
                format!("@left({} to {})", start_bucket, end_bucket + BUCKET_SIZE)
            } else {
                let right_start =
                    bucket_position(annotation.read_len.unwrap().saturating_sub(annotation.end));
                let right_end = bucket_position(
                    annotation
                        .read_len
                        .unwrap()
                        .saturating_sub(annotation.start),
                );
                format!("@right({} to {})", right_start, right_end + BUCKET_SIZE)
            };

            // Get the label of the annotation
            let label = annotation
                .label
                .label
                .clone()
                .unwrap_or("Flank".to_string());

            // Count the label
            if label != "Flank" {
                *label_map.entry(label.to_string()).or_insert(0) += 1;
            }

            // Get cut direction if present
            let cut = if let Some(cuts) = &annotation.cuts {
                if !cuts.is_empty() {
                    match cuts[0].direction {
                        CutDirection::Before => ", <<",
                        CutDirection::After => ", >>",
                    }
                } else {
                    ""
                }
                .to_string()
            } else {
                "".to_string()
            };

            // Get match type matching pattern defined in filter
            let match_type = match annotation.label.match_type {
                MatchType::Flank => "Flank",
                MatchType::Fbarcode => "Fbarcode",
                MatchType::Rbarcode => "Rbarcode",
            };

            format!(
                "{}[{}, *{}, {}]",
                match_type,
                if annotation.label.orientation == Orientation::Forward {
                    "fw"
                } else {
                    "rc"
                },
                cut,
                position
            )
        })
        .collect();

    // Pattern as string
    let mut pattern_string = pattern_elements.join("__");

    // Now we use incremental wirldcards to replace the labels
    let mut wildcard_count = 0;
    for (label, count) in label_map {
        if count > 1 {
            pattern_string = pattern_string.replace(&label, &format!("?{}", wildcard_count));
            wildcard_count += 1;
        }
    }

    pattern_string
}

pub fn inspect(annotated_file: &str, top_n: usize) -> Result<(), Box<dyn Error>> {
    println!("\n{}", "Inspecting".bold().underline());
    println!("  • Input:  {}", annotated_file.bold());

    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(annotated_file)
        .expect("Failed to open annotated file");

    // Process reads one group at a time
    let mut current_read_id: Option<String> = None;
    let mut current_group: Vec<BarbellMatch> = Vec::new();

    // Keep track of how often we see a pattern
    let mut pattern_count: HashMap<String, usize> = HashMap::new();

    for result in reader.deserialize() {
        let record: BarbellMatch = result?;

        if let Some(read_id) = &current_read_id {
            if *read_id != record.read.clone().unwrap() {
                // Process previous group
                let label = get_group_structure(&current_group);
                *pattern_count.entry(label).or_insert(0) += 1;
                current_group.clear();
                current_read_id = Some(record.read.clone().unwrap());
            }
        } else {
            current_read_id = Some(record.read.clone().unwrap());
        }

        current_group.push(record);
    }

    // Process the last group
    if !current_group.is_empty() {
        let label = get_group_structure(&current_group);
        *pattern_count.entry(label).or_insert(0) += 1;
    }

    // Show top 10 most common patterns
    let mut pattern_count_vec: Vec<(String, usize)> = pattern_count.into_iter().collect();
    pattern_count_vec.sort_by(|a, b| b.1.cmp(&a.1));
    for (pattern, count) in pattern_count_vec.iter().take(top_n) {
        println!("{}\nCount: {}\n", pattern, count.to_string().green());
    }

    Ok(())
}
