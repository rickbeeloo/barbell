use crate::barbell::filter::AnnotationLine;
use seq_io::fastq::{Reader,Record};
use csv;
use std::io::{BufWriter, Write};
use std::fs::File;
use std::collections::HashMap;
use crate::barbell::pattern_assign::*;
use std::path::Path;
use std::time::Instant;
use indicatif;
use indicatif::ProgressStyle;
use colored::Colorize;

pub struct LabelConfig {
    include_label: bool,
    include_orientation: bool,
    include_flank: bool,
}

impl LabelConfig {
    pub fn new(include_label: bool, include_orientation: bool, include_flank: bool) -> Self {
        Self { include_label, include_orientation, include_flank }
    }

    fn create_label(&self, annotations: &[AnnotationLine]) -> String {
        if !self.include_label {
            return "none".to_string();
        }

        let label_parts: Vec<String> = annotations.iter()
            .filter_map(|a| {
                let encoded = EncodedMatchStr::unstringify(&a.label);
                let label = encoded.label.unwrap_or_else(|| "Flank".to_string());
                
                // Skip if it's a flank and we don't want flanks
                if !self.include_flank && label == "Flank" {
                    return None;
                }

                let mut result = label;
                
                if self.include_orientation {
                    let ori = match encoded.orientation {
                        Orientation::Forward => "fw",
                        Orientation::ReverseComplement => "rc",
                    };
                    result = format!("{}_{}", result, ori);
                }
                
                Some(result)
            })
            .collect();

        if label_parts.is_empty() {
            "none".to_string()
        } else {
            label_parts.join("__")
        }
    }
}

#[derive(Debug)]
struct CompleteSlice {
    start: usize,
    end: usize,
    annotations: Vec<AnnotationLine>,  // Store all relevant annotations
}

fn preprocess_cuts(annotations: &[AnnotationLine], seq_len: usize) -> Vec<CompleteSlice> {
    let mut slices: Vec<CompleteSlice> = Vec::new();
    
    // Group cuts by their IDs
    let mut cut_groups: HashMap<usize, Vec<(usize, usize, &Cut, &AnnotationLine)>> = HashMap::new();
    for anno in annotations {
        if let Some(cuts) = &anno.cuts {
            println!("Cuts: {:?}", cuts);
            for cut in cuts {
                cut_groups.entry(cut.group_id)
                    .or_default()
                    .push((anno.start, anno.end, cut, anno));
            }
        }
    }
    // Sort groups by their leftmost position
    let mut sorted_groups: Vec<_> = cut_groups.into_iter().collect();
    sorted_groups.sort_by_key(|(_, group)| group.first().map(|(start, _, _, _)| *start).unwrap_or(usize::MAX));

    // Process each group
    for (i, (_, group)) in sorted_groups.iter().enumerate() {
        if group.len() == 2 {
            // We have two annotations so get start and end based on their cuts 
            println!("Got two annotations");
            let group1 = &group[0];
            let group2 = &group[1];
            
            // Get start position based on first group's cut direction
            let start = match &group1.2.direction {
                CutDirection::Before => group1.0,  // If first is Before, use second's start
                CutDirection::After => group1.1,   // If first is After, use first's start
            };
            
            // Get end position based on second group's cut direction
            let end = match &group2.2.direction {
                CutDirection::Before => group2.0,  // If second is Before, use second's end
                CutDirection::After => group2.1,   // If second is After, use first's end
            };

            println!("!!!!!!!!!!!!! Start: {}, End: {}", start, end);

            // Both annotations are relevant for this slice
            let annotations = vec![group1.3.clone(), group2.3.clone()];

            slices.push(CompleteSlice {
                start,
                end,
                annotations,
            });
        } else if group.len() == 1 {
            println!("Group length is 1");
            let &(start, end, cut, anno) = &group[0];
            
            match cut.direction {
                CutDirection::Before => {
                    // Look left for start position and annotation
                    let (slice_start, left_anno) = if i > 0 {
                        let prev_group = &sorted_groups[i-1].1;
                        let max_end_idx = prev_group.iter()
                            .enumerate()
                            .max_by_key(|(_, (_, end, _, _))| end)
                            .map(|(idx, _)| idx)
                            .unwrap_or(0);
                        (prev_group[max_end_idx].1, Some(prev_group[max_end_idx].3.clone()))
                    } else {
                        (0, None)
                    };
                    
                    let mut annotations = Vec::new();
                    if let Some(left) = left_anno {
                        annotations.push(left);
                    }
                    annotations.push(anno.clone());

                    slices.push(CompleteSlice {
                        start: slice_start,
                        end,
                        annotations,
                    });
                },
                CutDirection::After => {
                    // Look right for end position and annotation
                    let (slice_end, right_anno) = if i < sorted_groups.len() - 1 {
                        let next_group = &sorted_groups[i+1].1;
                        let min_start_idx = next_group.iter()
                            .enumerate()
                            .min_by_key(|(_, (start, _, _, _))| start)
                            .map(|(idx, _)| idx)
                            .unwrap_or(0);
                        (next_group[min_start_idx].0, Some(next_group[min_start_idx].3.clone()))
                    } else {
                        println!("No next group, should slice till end");
                        (seq_len, None)
                    };
                    
                    let mut annotations = Vec::new();
                    annotations.push(anno.clone());
                    if let Some(right) = right_anno {
                        annotations.push(right);
                    }

                    println!("Adding slice with start: {}, end: {}", start, slice_end);

                    slices.push(CompleteSlice {
                        start,
                        end: slice_end,
                        annotations,
                    });
                }
            }
        }
    }

    slices
}

pub fn process_read_and_anno(seq: &[u8], qual: &[u8], annotations: &[AnnotationLine], label_config: &LabelConfig) -> Vec<(Vec<u8>, Vec<u8>, String, String)> {
    let mut results = Vec::new();
    let seq_len = seq.len();

    // Preprocess cuts to get complete slices
    let slices = preprocess_cuts(annotations, seq_len);


    // Group slices by cut group ID
    for slice in &slices {
        //if slice.start < slice.end && slice.end <= seq_len {
            println!(">>>> Start: {}, End: {}", slice.start, slice.end);
            let trimmed_seq = seq[slice.start..slice.end].to_vec();
            let trimmed_qual = qual[slice.start..slice.end].to_vec();

            let label_matches: Vec<AnnotationLine> = slice.annotations.clone();
            
            let group_label = label_config.create_label(&label_matches);
            let read_suffix = format!("_{}", group_label);
            results.push((trimmed_seq, trimmed_qual, group_label, read_suffix));
       // } else {
        //     println!("Skipping invalid positions for group {:?}: start={}, end={}", slice.annotations, slice.start, slice.end);
        // }
    }

    // println!("Results: {:?}", results);
    results
}

/// Extracts the clean read ID from a FASTQ record ID by taking the first part before any whitespace
fn clean_read_id(id: &str) -> &str {
    id.split_whitespace().next().unwrap_or(id)
}

pub fn trim_matches(filtered_match_file: &str, read_fastq_file: &str, output_folder: &str, add_labels: bool, add_orientation: bool, add_flank: bool) {
    
    let start_time = Instant::now();

    println!("\n{}", "Trimming".bold().underline());
    println!("  • Input matches: {}", filtered_match_file.bold());
    println!("  • Input reads:   {}", read_fastq_file.bold());
    println!("  • Output folder: {}", output_folder.bold());
    println!("  • Labels:        {}", if add_labels { "yes".green() } else { "no".red() });
    println!("  • Orientation:   {}", if add_orientation { "yes".green() } else { "no".red() });
    println!("  • Flanks:        {}\n", if add_flank { "yes".green() } else { "no".red() });

    let mut reader = Reader::from_path(read_fastq_file).unwrap();
    let mut writers: HashMap<String, BufWriter<File>> = HashMap::new();

    // If output folder does not exist, create it
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
    }

    // Label formatting config
    let label_config = LabelConfig::new(add_labels, add_orientation, add_flank);
    
    // Setup progress bar
    let progress_bar = indicatif::ProgressBar::new_spinner();
    progress_bar.set_style(
        ProgressStyle::with_template("{spinner:.blue} {prefix:<12} {msg:>6} {elapsed_precise}")
            .unwrap()
            .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
    );
    progress_bar.set_prefix("Processing:");

    let mut total_reads = 0;
    let mut mapped_reads = 0;
    let mut trimmed_reads = 0;
    let mut trimmed_split_reads = 0;
    
    let mut matches_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(filtered_match_file)
        .expect("Failed to open matches file");

    // If output folder does not exist, create it
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
    }

    // Group annotations by read ID while preserving order
    let mut current_read_id: Option<Vec<u8>> = None;
    let mut current_annotations: Vec<AnnotationLine> = Vec::new();

    for result in matches_reader.deserialize() {
        let anno: AnnotationLine = result.expect("Failed to parse annotation line");

        let read_id = anno.read.clone().into_bytes();

        // If we encounter a new read ID
        if current_read_id.as_ref() != Some(&read_id) {
            // Process previous group if it exists
            if let Some(prev_read_id) = current_read_id.take() {

                // Search through FASTQ records until we find matching read
                while let Some(record) = reader.next() {
                    let record = record.expect("Error reading record");
                    total_reads += 1;
                    
 
                    //todo! why we need this, seq io should handle this 
                    let record_id = clean_read_id(record.id().unwrap());

                    
                    if record_id.as_bytes() == prev_read_id.as_slice() {
                        mapped_reads += 1;
                        let results = process_read_and_anno(
                            record.seq(),
                            record.qual(),
                            &current_annotations,
                            &label_config
                        );

                        if !results.is_empty() {
                            trimmed_reads += 1; // Count once for multi split
                        }
                        
                        for (trimmed_seq, trimmed_qual, group, _) in results {
                            trimmed_split_reads += 1;
                            
                            // Get or create writer for this group
                            let writer = writers.entry(group.clone()).or_insert_with(|| {
                                let output_file = format!("{}/{}.trimmed.fastq", output_folder, group);
                                BufWriter::new(File::create(&output_file)
                                    .unwrap_or_else(|_| panic!("Failed to create output file: {}", output_file)))
                            });
                            
                            // Write FASTQ format
                            writeln!(writer, "@{}", String::from_utf8_lossy(clean_read_id(record.id().unwrap()).as_bytes()))
                                .expect("Failed to write header");
                            writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_seq))
                                .expect("Failed to write sequence");
                            writeln!(writer, "+").expect("Failed to write separator");
                            writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_qual))
                                .expect("Failed to write quality");
                        }
                        break;
                    }
                }
            }
            
            // Start new group
            current_read_id = Some(read_id);
            current_annotations.clear();
        }
        
        current_annotations.push(anno);
    }

    // Process the last group
    if let Some(last_read_id) = current_read_id {
        println!("Processing last group");
        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            total_reads += 1;
            
            if record.id().unwrap().as_bytes() == last_read_id.as_slice() {
                mapped_reads += 1;
                let results = process_read_and_anno(
                    record.seq(),
                    record.qual(),
                    &current_annotations,
                    &label_config
                );

                if !results.is_empty() {
                    trimmed_reads += 1;
                }
                
                for (trimmed_seq, trimmed_qual, group, _) in results {
                    trimmed_split_reads += 1;
                    
                    // Get or create writer for this group
                    let writer = writers.entry(group.clone()).or_insert_with(|| {
                        let output_file = format!("{}/{}.trimmed.fastq", output_folder, group);
                        BufWriter::new(File::create(&output_file)
                            .expect(&format!("Failed to create output file: {}", output_file)))
                    });
                    
                    // Write FASTQ format
                    writeln!(writer, "@{}", String::from_utf8_lossy(clean_read_id(record.id().unwrap()).as_bytes()))
                        .expect("Failed to write header");
                    writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_seq))
                        .expect("Failed to write sequence");
                    writeln!(writer, "+").expect("Failed to write separator");
                    writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_qual))
                        .expect("Failed to write quality");
                }
                break;
            }
        }
    }

    // Flush all writers
    for (_, writer) in writers.iter_mut() {
        writer.flush().expect("Failed to flush output");
    }

    // Finish progress bar and print summary
    progress_bar.finish_with_message("Done!");

    // Split is the difference between the trimmed reads and the 
    let extra_splits = trimmed_split_reads - trimmed_reads;

    println!("\n{}", "Summary".bold().underline());
    println!("  • Time: {} seconds", start_time.elapsed().as_secs().to_string().bold());
    println!("  • Total reads: {}", total_reads.to_string().bold());
    println!("  • Mapped reads: {}", mapped_reads.to_string().green().bold());
    println!("  • Trimmed reads: {}", trimmed_reads.to_string().green().bold());
    println!("  • Filter (keep rate): {}%", format!("{:.2}", (mapped_reads as f64 / total_reads as f64) * 100.0).bold());
    println!("  • Trimming rate: {}%", format!("{:.2}", (trimmed_reads as f64 / mapped_reads as f64) * 100.0).bold());
    println!("  • Extra splits: {}", extra_splits.to_string().blue().bold());
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::barbell::filter::AnnotationLine;
//     use crate::barbell::pattern_assign::{Cut, CutDirection, Orientation};

//     #[test]
//     fn test_single_cut() {
//         let seq =  b"CCCCCCCCAAAACCCCCCCCCCCC";
//         let qual = b"________IIII____________";

//         let annotations = vec![
            
//             AnnotationLine {
//                 read: "read1".to_string(),
//                 start: 4,
//                 end: 8,
//                 label: "Fbar#fw#Fbar".to_string(),
//                 cuts: Some(vec![Cut::new(0, CutDirection::After)]),
//                 dist_to_end: 0,
//                 edits: 0,
//                 record_set_idx: 0,
//                 record_idx: 0,
//                 read_len: seq.len(),
//             },

//             AnnotationLine {
//                 read: "read1".to_string(),
//                 start: 12,
//                 end: 16,
//                 label: "Rbar#fw#Rbar".to_string(),
//                 cuts: Some(vec![Cut::new(0, CutDirection::Before)]),
//                 dist_to_end: 0,
//                 edits: 0,
//                 record_set_idx: 0,
//                 record_idx: 0,
//                 read_len: seq.len(),
//             },
//         ];

//         let label_config = LabelConfig::new(true, true, true);
//         let results = process_read_and_anno(seq, qual, &annotations, &label_config);

//         assert_eq!(results.len(), 1);
//         let (trimmed_seq, trimmed_qual, group_label, _) = &results[0];
//         println!("trimmed_seq: {}", String::from_utf8_lossy(trimmed_seq));
//         assert_eq!(trimmed_seq, b"AAAA");
//         assert_eq!(trimmed_qual, b"IIII");
//         assert_eq!(group_label, "Fbar_fw__Rbar_fw");
//     }

//     #[test]
//     fn test_multiple_cuts() {
//         let seq = b"CCCCCCCCAAAAACCCCCCCCCTTT";
//         let qual = b"________IIII____________";

//         let annotations = vec![
//             AnnotationLine {
//                 read: "read2".to_string(),
//                 start: 4,
//                 end: 8,
//                 label: "Fbar#fw#Fbar".to_string(),
//                 cuts: Some(vec![Cut::new(1, CutDirection::After)]),
//                 dist_to_end: 0,
//                 edits: 0,
//                 record_set_idx: 0,
//                 record_idx: 0,
//                 read_len: seq.len(),
//             },
//             AnnotationLine {
//                 read: "read2".to_string(),
//                 start: 12,
//                 end: 16,
//                 label: "Rbar#fw#Rbar".to_string(),
//                 cuts: Some(vec![Cut::new(1, CutDirection::Before)]),
//                 dist_to_end: 0,
//                 edits: 0,
//                 record_set_idx: 0,
//                 record_idx: 0,
//                 read_len: seq.len(),
//             },
//             AnnotationLine {
//                 read: "read2".to_string(),
//                 start: 18,
//                 end: 22,
//                 label: "Fbar#fw#Fbar".to_string(),
//                 cuts: Some(vec![Cut::new(2, CutDirection::After)]),
//                 dist_to_end: 0,
//                 edits: 0,
//                 record_set_idx: 0,
//                 record_idx: 0,
//                 read_len: seq.len(),
//             },
//             AnnotationLine {
//                 read: "read2".to_string(),
//                 start: 24,
//                 end: 28,
//                 label: "Rbar#fw#Rbar".to_string(),
//                 cuts: Some(vec![Cut::new(2, CutDirection::Before)]),
//                 dist_to_end: 0,
//                 edits: 0,
//                 record_set_idx: 0,
//                 record_idx: 0,
//                 read_len: seq.len(),
//             },
//         ];

//         let label_config = LabelConfig::new(true, true, true);
//         let results = process_read_and_anno(seq, qual, &annotations, &label_config);

//         println!("results: {:?}", results);
//         assert_eq!(results.len(), 2);

//         let (trimmed_seq1, trimmed_qual1, group_label1, _) = &results[0];
//         assert_eq!(trimmed_seq1, b"AAAA");
//         assert_eq!(trimmed_qual1, b"IIII");
//         assert_eq!(group_label1, "Fbar_fw__Rbar_fw");

//         let (trimmed_seq2, trimmed_qual2, group_label2, _) = &results[1];
//         assert_eq!(trimmed_seq2, b"ACGTACGT");
//         assert_eq!(trimmed_qual2, b"IIIIIIII");
//         assert_eq!(group_label2, "Fbar_fw__Rbar_fw");
//     }
// }