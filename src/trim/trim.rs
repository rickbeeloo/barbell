use crate::annotate::searcher::BarbellMatch;
use crate::filter::pattern::{Cut, CutDirection};
use colored::Colorize;
use csv;
use indicatif::{ProgressBar, ProgressStyle};
use sassy::Strand;
use seq_io::fastq::{Reader, Record};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

pub struct LabelConfig {
    include_label: bool,
    include_orientation: bool,
    include_flank: bool,
    sort_labels: bool,
}

impl LabelConfig {
    pub fn new(
        include_label: bool,
        include_orientation: bool,
        include_flank: bool,
        sort_labels: bool,
    ) -> Self {
        Self {
            include_label,
            include_orientation,
            include_flank,
            sort_labels,
        }
    }

    fn create_label(&self, annotations: &[BarbellMatch]) -> String {
        if !self.include_label {
            return "none".to_string();
        }

        let mut label_parts: Vec<String> = annotations
            .iter()
            .filter_map(|m| {
                let label = m.label.clone();

                // Skip if it's a flank and we don't want flanks
                if !self.include_flank && label.contains("flank") {
                    return None;
                }

                let mut result = label;

                if self.include_orientation {
                    let ori = match m.strand {
                        Strand::Fwd => "fw",
                        Strand::Rc => "rc",
                    };
                    result = format!("{result}_{ori}");
                }

                Some(result)
            })
            .collect();

        if label_parts.is_empty() {
            "none".to_string()
        } else if self.sort_labels {
            label_parts.sort();
            label_parts.join("__")
        } else {
            label_parts.join("__")
        }
    }
}

#[derive(Debug)]
struct CompleteSlice {
    start: usize,
    end: usize,
    annotations: Vec<BarbellMatch>,
}

fn preprocess_cuts(annotations: &[BarbellMatch], seq_len: usize) -> Vec<CompleteSlice> {
    let mut slices: Vec<CompleteSlice> = Vec::new();

    // Group cuts by their IDs
    let mut cut_groups: HashMap<usize, Vec<(usize, usize, &Cut, &BarbellMatch)>> = HashMap::new();
    for anno in annotations {
        if let Some(cuts) = &anno.cuts {
            for (cut, _) in cuts {
                cut_groups.entry(cut.group_id).or_default().push((
                    anno.read_start_flank,
                    anno.read_end_flank,
                    cut,
                    anno,
                ));
            }
        }
    }

    // Sort groups by their leftmost position
    let mut sorted_groups: Vec<_> = cut_groups.into_iter().collect();
    sorted_groups.sort_by_key(|(_, group)| {
        group
            .first()
            .map(|(start, _, _, _)| *start)
            .unwrap_or(usize::MAX)
    });

    // Process each group
    for (i, (_, group)) in sorted_groups.iter().enumerate() {
        if group.len() == 2 {
            // We have two annotations so get start and end based on their cuts
            let group1 = &group[0];
            let group2 = &group[1];

            // Get start position based on first group's cut direction
            let start = match &group1.2.direction {
                CutDirection::Before => group1.0,
                CutDirection::After => group1.1,
            };

            // Get end position based on second group's cut direction
            let end = match &group2.2.direction {
                CutDirection::Before => group2.0,
                CutDirection::After => group2.1,
            };

            let annotations = vec![group1.3.clone(), group2.3.clone()];

            slices.push(CompleteSlice {
                start,
                end,
                annotations,
            });
        } else if group.len() == 1 {
            let &(start, end, cut, anno) = &group[0];

            match cut.direction {
                CutDirection::Before => {
                    // Look left for start position and annotation
                    let (slice_start, left_anno) = if i > 0 {
                        let prev_group = &sorted_groups[i - 1].1;
                        let max_end_idx = prev_group
                            .iter()
                            .enumerate()
                            .max_by_key(|(_, (_, end, _, _))| end)
                            .map(|(idx, _)| idx)
                            .unwrap_or(0);
                        (
                            prev_group[max_end_idx].1,
                            Some(prev_group[max_end_idx].3.clone()),
                        )
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
                        end: start,
                        annotations,
                    });
                }
                CutDirection::After => {
                    // Look right for end position and annotation
                    let (slice_end, right_anno) = if i < sorted_groups.len() - 1 {
                        let next_group = &sorted_groups[i + 1].1;
                        let min_start_idx = next_group
                            .iter()
                            .enumerate()
                            .min_by_key(|(_, (start, _, _, _))| start)
                            .map(|(idx, _)| idx)
                            .unwrap_or(0);
                        (
                            next_group[min_start_idx].0,
                            Some(next_group[min_start_idx].3.clone()),
                        )
                    } else {
                        (seq_len, None)
                    };

                    let mut annotations = Vec::new();
                    annotations.push(anno.clone());
                    if let Some(right) = right_anno {
                        annotations.push(right);
                    }

                    slices.push(CompleteSlice {
                        start: end,
                        end: slice_end,
                        annotations,
                    });
                }
            }
        }
    }

    slices
}

pub fn process_read_and_anno(
    seq: &[u8],
    qual: &[u8],
    annotations: &[BarbellMatch],
    label_config: &LabelConfig,
) -> Vec<(Vec<u8>, Vec<u8>, String, String)> {
    let mut results = Vec::new();
    let seq_len = seq.len();

    // Preprocess cuts to get complete slices
    let slices = preprocess_cuts(annotations, seq_len);

    // Group slices by cut group ID
    for slice in &slices {
        if slice.start > slice.end {
            continue;
        }

        if slice.start == slice.end {
            continue;
        }

        let trimmed_seq = seq[slice.start..slice.end].to_vec();
        let trimmed_qual = qual[slice.start..slice.end].to_vec();

        let label_matches: Vec<BarbellMatch> = slice.annotations.clone();

        let group_label = label_config.create_label(&label_matches);
        let read_suffix = format!("_{group_label}");
        results.push((trimmed_seq, trimmed_qual, group_label, read_suffix));
    }

    results
}

/// Extracts the clean read ID from a FASTQ record ID by taking the first part before any whitespace
fn clean_read_id(id: &str) -> &str {
    id.split_whitespace().next().unwrap_or(id)
}

pub fn trim_matches(
    filtered_match_file: &str,
    read_fastq_file: &str,
    output_folder: &str,
    add_labels: bool,
    add_orientation: bool,
    add_flank: bool,
    sort_labels: bool,
    _n_threads: usize,
) {
    // Create output folder if it doesn't exist
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
    }

    // Label formatting config
    let label_config = LabelConfig::new(add_labels, add_orientation, add_flank, sort_labels);

    // Read all annotations and group by read ID
    let mut annotations_by_read: HashMap<String, Vec<BarbellMatch>> = HashMap::new();

    let mut matches_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(filtered_match_file)
        .expect("Failed to open matches file");

    let mut annotation_count = 0;
    for result in matches_reader.deserialize() {
        let anno: BarbellMatch = result.expect("Failed to parse annotation line");
        annotations_by_read
            .entry(anno.read_id.clone())
            .or_default()
            .push(anno);
        annotation_count += 1;
    }

    // Create writers
    let mut writers: HashMap<String, BufWriter<File>> = HashMap::new();

    // Create progress tracking
    let mut total_reads = 0;
    let mut mapped_reads = 0;
    let mut trimmed_reads = 0;
    let mut trimmed_split_reads = 0;

    // Create progress bar
    let progress_bar = ProgressBar::new(annotations_by_read.len() as u64);
    progress_bar.set_style(
        ProgressStyle::with_template("{spinner:.blue} {prefix:<12} [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) {elapsed_precise}")
            .unwrap()
            .progress_chars("#>-")
    );
    progress_bar.set_prefix("Processing");

    // Process reads
    let mut reader = Reader::from_path(read_fastq_file).expect("Failed to open FASTQ file");

    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        let read_id = clean_read_id(record.id().unwrap()).to_string();
        total_reads += 1;

        // Check if this read has annotations
        if let Some(annotations) = annotations_by_read.get(&read_id) {
            mapped_reads += 1;

            let results =
                process_read_and_anno(record.seq(), record.qual(), annotations, &label_config);

            if !results.is_empty() {
                trimmed_reads += 1;
            }

            for (trimmed_seq, trimmed_qual, group, _) in results {
                trimmed_split_reads += 1;

                // Get or create writer for this group
                let writer = writers.entry(group.clone()).or_insert_with(|| {
                    let output_file = format!("{}/{}.trimmed.fastq", output_folder, group);
                    BufWriter::new(File::create(&output_file).unwrap_or_else(|_| {
                        panic!("Failed to create output file: {}", output_file)
                    }))
                });

                // Write FASTQ format
                writeln!(writer, "@{}", read_id).expect("Failed to write header");
                writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_seq))
                    .expect("Failed to write sequence");
                writeln!(writer, "+").expect("Failed to write separator");
                writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_qual))
                    .expect("Failed to write quality");
            }

            progress_bar.inc(1);
        }
    }

    // Flush all writers
    for (group, writer) in writers.iter_mut() {
        writer.flush().expect("Failed to flush output");
    }

    // Finish progress bar and print summary
    progress_bar.finish_with_message("Done!");
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotate::barcodes::BarcodeType;
    use crate::filter::pattern::{Cut, CutDirection};

    #[test]
    fn test_single_cut() {
        let seq = b"CCCCCCCCAAAACCCCCCCCCCCC";
        let qual = b"________IIII____________";

        let annotations = vec![
            BarbellMatch::new(
                4, // read_start_bar
                8, // read_end_bar
                4, // read_start_flank
                8, // read_end_flank
                0, // bar_start
                4, // bar_end
                BarcodeType::Fbar,
                0, // flank_cost
                0, // barcode_cost
                "Fbar".to_string(),
                Strand::Fwd,
                seq.len(),
                "read1".to_string(),
                0, // rel_dist_to_end
                Some(vec![(Cut::new(0, CutDirection::After), 8)]),
            ),
            BarbellMatch::new(
                12, // read_start_bar
                16, // read_end_bar
                12, // read_start_flank
                16, // read_end_flank
                0,  // bar_start
                4,  // bar_end
                BarcodeType::Rbar,
                0, // flank_cost
                0, // barcode_cost
                "Rbar".to_string(),
                Strand::Fwd,
                seq.len(),
                "read1".to_string(),
                0, // rel_dist_to_end
                Some(vec![(Cut::new(0, CutDirection::Before), 12)]),
            ),
        ];

        let label_config = LabelConfig::new(true, true, true, true);
        let results = process_read_and_anno(seq, qual, &annotations, &label_config);

        assert_eq!(results.len(), 1);
        let (trimmed_seq, trimmed_qual, group_label, _) = &results[0];
        println!("trimmed_seq: {}", String::from_utf8_lossy(trimmed_seq));
        assert_eq!(trimmed_seq, b"AAAA");
        assert_eq!(trimmed_qual, b"IIII");
        assert_eq!(group_label, "Fbar_fw__Rbar_fw");
    }
}
