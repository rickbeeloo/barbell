use crate::barbell::filter::AnnotationLine;
use crate::barbell::merge_sort::merge_sort_files;
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

pub fn process_read_and_anno(seq: &[u8], qual: &[u8], annotations: &[AnnotationLine], label_config: &LabelConfig) -> Option<(Vec<u8>, Vec<u8>, String)> {
    let cuts = annotations.first().unwrap().cuts.clone().unwrap_or_default();
    let group = label_config.create_label(annotations);
    
    let (trimmed_seq, trimmed_qual) = match cuts.len() {
        0 => {
            // No cuts means no trimming needed
            return None;
        },
        1 => {
            let cut = &cuts[0];
            let cut_pos = cut.position;
            
            // Check if cut position is valid
            if cut_pos >= seq.len() {
                return None;
            }
            
            match cut.direction {
                CutDirection::Before => (seq[..cut_pos].to_vec(), qual[..cut_pos].to_vec()),
                CutDirection::After => (seq[cut_pos..].to_vec(), qual[cut_pos..].to_vec()),
            }
        },
        2 => {
            let (cut1, cut2) = (&cuts[0], &cuts[1]);
            
            let (start, end) = match (cut1.direction.clone(), cut2.direction.clone()) {
                (CutDirection::After, CutDirection::Before) => (cut1.position, cut2.position),
                (CutDirection::Before, CutDirection::After) => (cut2.position, cut1.position),
                _ => return None, // Invalid cut directions
            };
            
            // Check if positions are valid
            if start >= end || start >= seq.len() || end > seq.len() {
                return None;
            }
            
            // Check if resulting sequence would be empty
            if end - start < 1 {
                return None;
            }
            
            (seq[start..end].to_vec(), qual[start..end].to_vec())
        },
        _ => return None, // More than 2 cuts not supported
    };

    // Check if we actually trimmed something
    if trimmed_seq.is_empty() {
        return None;
    }

    Some((trimmed_seq, trimmed_qual, group))
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
                    
                    if record.id().unwrap().as_bytes() == prev_read_id.as_slice() {
                        mapped_reads += 1;
                        if let Some((trimmed_seq, trimmed_qual, group)) = process_read_and_anno(
                            record.seq(),
                            record.qual(),
                            &current_annotations,
                            &label_config
                        ) {
                            trimmed_reads += 1;
                            
                            // Get or create writer for this group
                            let writer = writers.entry(group.clone()).or_insert_with(|| {
                                let output_file = format!("{}/{}.trimmed.fastq", output_folder, group);
                                BufWriter::new(File::create(&output_file)
                                    .expect(&format!("Failed to create output file: {}", output_file)))
                            });
                            
                            // Write FASTQ format
                            writeln!(writer, "@{}", String::from_utf8_lossy(record.id().unwrap().as_bytes()))
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
        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            total_reads += 1;
            
            if record.id().unwrap().as_bytes() == last_read_id.as_slice() {
                mapped_reads += 1;
                if let Some((trimmed_seq, trimmed_qual, group)) = process_read_and_anno(
                    record.seq(),
                    record.qual(),
                    &current_annotations,
                    &label_config
                ) {
                    trimmed_reads += 1;
                    
                    // Get or create writer for this group
                    let writer = writers.entry(group.clone()).or_insert_with(|| {
                        let output_file = format!("{}/{}.trimmed.fastq", output_folder, group);
                        BufWriter::new(File::create(&output_file)
                            .expect(&format!("Failed to create output file: {}", output_file)))
                    });
                    
                    // Write FASTQ format
                    writeln!(writer, "@{}", String::from_utf8_lossy(record.id().unwrap().as_bytes()))
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

    println!("\n{}", "Summary".bold().underline());
    println!("  • Time: {} seconds", start_time.elapsed().as_secs().to_string().bold());
    println!("  • Total reads: {}", total_reads.to_string().bold());
    println!("  • Mapped reads: {}", mapped_reads.to_string().green().bold());
    println!("  • Trimmed reads: {}", trimmed_reads.to_string().green().bold());
    println!("  • Filter (keep rate): {}%", format!("{:.2}", (mapped_reads as f64 / total_reads as f64) * 100.0).bold());
    println!("  • Trimming rate: {}%\n", format!("{:.2}", (trimmed_reads as f64 / mapped_reads as f64) * 100.0).bold());
}