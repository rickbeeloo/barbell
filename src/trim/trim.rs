use crate::annotate::barcodes::BarcodeType;
use crate::annotate::searcher::BarbellMatch;
use crate::filter::pattern::{Cut, CutDirection};
use crate::io::io::open_fastq;
use crate::progress::progress::{ProgressTracker, TRIM_PROGRESS_SPECS};
use anyhow::anyhow;
use csv;
use sassy::Strand;
use seq_io::fastq::Record;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use clap::ValueEnum;

const TOTAL_IDX: usize = 0;
const TRIMMED_IDX: usize = 1;
const TRIMMED_SPLIT_IDX: usize = 2;
const FAILED_IDX: usize = 3;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum LabelSide {
    Left,
    Right,
}

pub struct LabelConfig {
    pub include_label: bool,
    pub include_orientation: bool,
    pub include_flank: bool,
    pub sort_labels: bool,
    pub only_side: Option<LabelSide>,
}

impl LabelConfig {
    pub fn new(
        include_label: bool,
        include_orientation: bool,
        include_flank: bool,
        sort_labels: bool,
        only_side: Option<LabelSide>,
    ) -> Self {
        Self {
            include_label,
            include_orientation,
            include_flank,
            sort_labels,
            only_side,
        }
    }

    pub fn create_label(&self, annotations: &[BarbellMatch]) -> String {
        if !self.include_label {
            return "none".to_string();
        }

        let mut label_parts: Vec<String> = annotations
            .iter()
            .filter_map(|m| {
                let label = m.label.clone();

                // Skip if it's a flank and we don't want flanks
                // this also prevents having a flank come before label when
                // selecting only left or only right
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

        if self.sort_labels && self.only_side.is_some() {
            panic!("Cannot enable only keeping left label and sorting as this makes it ambiguous");
        }

        if label_parts.is_empty() {
            "none".to_string()
        } else if self.sort_labels {
            label_parts.sort();
            label_parts.join("__")
        } else if self.only_side.is_some() {
            let side = self.only_side.unwrap();
            if side == LabelSide::Left {
                label_parts.first().unwrap().clone()
            } else {
                label_parts.last().unwrap().clone()
            }
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
    skip_trim: bool,
    flip: bool,
) -> Vec<(Vec<u8>, Vec<u8>, String, String)> {
    let mut results = Vec::new();
    let seq_len = seq.len();

    // Preprocess cuts to get complete slices
    let slices = preprocess_cuts(annotations, seq_len);

    // Group slices by cut group ID
    for slice in &slices {
        if slice.start >= slice.end {
            continue;
        }

        // For now if trimming is disabled, we just
        // return the full sequence and quality
        let mut trimmed_seq = if skip_trim {
            seq.to_vec()
        } else {
            seq[slice.start..slice.end].to_vec()
        };
        let mut trimmed_qual = if skip_trim {
            qual.to_vec()
        } else {
            qual[slice.start..slice.end].to_vec()
        };

        if flip && should_flip(&slice.annotations) {
            trimmed_seq = reverse_complement(&trimmed_seq);
            trimmed_qual.reverse();
        }

        let label_matches: Vec<BarbellMatch> = slice.annotations.clone();

        let group_label = label_config.create_label(&label_matches);
        let read_suffix = format!("_{group_label}");
        results.push((trimmed_seq, trimmed_qual, group_label, read_suffix));
    }

    results
}

/// Extracts the clean read ID from a FASTQ record ID by taking the first part before any whitespace
#[allow(unused)]
fn clean_read_id(id: &str) -> &str {
    id.split_whitespace().next().unwrap_or(id)
}

fn format_output_file_error(output_file: &str, err: &std::io::Error) -> String {
    let mut msg = format!("Failed to create output file '{output_file}': {err}");
    if err.raw_os_error() == Some(24) {
        msg.push_str("\nTry setting ulimit higher: \"ulimit -n 65000\"");
    }
    msg
}

fn should_flip(annotations: &[BarbellMatch]) -> bool {
    // If we matched an Ftag in rc we flip
    annotations
        .iter()
        .any(|anno| anno.match_type == BarcodeType::Ftag && anno.strand == Strand::Rc)
}

pub fn trim_matches(
    filtered_match_file: &str,
    read_fastq_file: &str,
    output_folder: &str,
    add_labels: bool,
    add_orientation: bool,
    add_flank: bool,
    sort_labels: bool,
    only_side: Option<LabelSide>,
    failed_trimmed_writer: Option<String>, // if provided we write ids of failed trimmed reads to this file, like empty reads
    write_full_header: bool,
    skip_trim: bool,
    // Experimental, if any Ftag = rc, flip it
    flip: bool,
    verbose: bool,
) -> anyhow::Result<()> {
    // Create output folder if it doesn't exist
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
    }

    // Label formatting config
    let label_config = LabelConfig::new(
        add_labels,
        add_orientation,
        add_flank,
        sort_labels,
        only_side,
    );

    if sort_labels && only_side.is_some() {
        return Err(anyhow!(
            "Cannot enable only keeping left/right label and sorting; this is ambiguous"
        ));
    }

    // Read all annotations and group by read ID
    let mut annotations_by_read: HashMap<String, Vec<BarbellMatch>> = HashMap::new();

    // Create progress bars
    let progress = if verbose {
        ProgressTracker::new_with_logging(&TRIM_PROGRESS_SPECS, "trim", output_folder)
    } else {
        ProgressTracker::new(&TRIM_PROGRESS_SPECS)
    };

    let mut matches_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(filtered_match_file)
        .expect("Failed to open matches file");

    for result in matches_reader.deserialize() {
        let anno: BarbellMatch = result.expect("Failed to parse annotation line");
        annotations_by_read
            .entry(anno.read_id.clone())
            .or_default()
            .push(anno);
    }

    // Create writers
    let mut writers: HashMap<String, BufWriter<File>> = HashMap::new();

    // If there is a failed trimmed writer, create it
    let mut failed_trimmed_writer = failed_trimmed_writer.map(|failed_trimmed_writer_path| {
        BufWriter::new(File::create(failed_trimmed_writer_path).unwrap())
    });

    // Process reads
    let mut reader = open_fastq(read_fastq_file);

    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        let (read_id, desc) = record.id_desc().unwrap();
        let read_id = read_id.to_string();
        let desc: &str = desc.unwrap_or_default();
        let full_header = format!("{read_id} {desc}");
        progress.inc(TOTAL_IDX);

        // Check if this read has annotations
        if let Some(annotations) = annotations_by_read.get(&read_id) {
            // mapped_reads += 1;

            let results: Vec<(Vec<u8>, Vec<u8>, String, String)> = process_read_and_anno(
                record.seq(),
                record.qual(),
                annotations,
                &label_config,
                skip_trim,
                flip,
            );

            if !results.is_empty() {
                progress.inc(TRIMMED_IDX);
            } else {
                progress.inc(FAILED_IDX);
                if let Some(ref mut failed_trimmed_writer) = failed_trimmed_writer {
                    writeln!(failed_trimmed_writer, "{read_id}")
                        .expect("Failed to write to failed trimmed writer");
                }
            }

            if results.len() > 1 {
                progress.inc(TRIMMED_SPLIT_IDX);
            }

            for (trimmed_seq, trimmed_qual, group, _) in results {
                // Get or create writer for this group
                if !writers.contains_key(&group) {
                    let output_file = format!("{output_folder}/{group}.trimmed.fastq");
                    let file = File::create(&output_file).map_err(|err| {
                        let msg = format_output_file_error(&output_file, &err);
                        progress.print_error(msg.clone());
                        progress.clear();
                        anyhow!(msg)
                    })?;
                    writers.insert(group.clone(), BufWriter::new(file));
                }
                let writer = writers
                    .get_mut(&group)
                    .expect("writer should exist after insertion");

                // Write FASTQ format
                if write_full_header {
                    writeln!(writer, "@{full_header}").expect("Failed to write header");
                } else {
                    writeln!(writer, "@{read_id}").expect("Failed to write header");
                }
                writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_seq))
                    .expect("Failed to write sequence");
                writeln!(writer, "+").expect("Failed to write separator");
                writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_qual))
                    .expect("Failed to write quality");
            }
        }

        progress.refresh();
    }

    // Flush all writers
    for (_, writer) in writers.iter_mut() {
        writer.flush().expect("Failed to flush output");
    }

    progress.finish("reads");
    Ok(())
}

#[inline(always)]
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&c| RC[c as usize]).collect()
}

const RC: [u8; 256] = {
    let mut rc = [0; 256];
    let mut i = 0;
    while i < 256 {
        rc[i] = i as u8;
        i += 1;
    }
    // Standard bases
    rc[b'A' as usize] = b'T';
    rc[b'C' as usize] = b'G';
    rc[b'T' as usize] = b'A';
    rc[b'G' as usize] = b'C';
    rc[b'a' as usize] = b't';
    rc[b'c' as usize] = b'g';
    rc[b't' as usize] = b'a';
    rc[b'g' as usize] = b'c';
    // IUPAC ambiguity codes
    rc[b'R' as usize] = b'Y'; // A|G -> T|C
    rc[b'Y' as usize] = b'R'; // C|T -> G|A
    rc[b'S' as usize] = b'S'; // G|C -> C|G
    rc[b'W' as usize] = b'W'; // A|T -> T|A
    rc[b'K' as usize] = b'M'; // G|T -> C|A
    rc[b'M' as usize] = b'K'; // A|C -> T|G
    rc[b'B' as usize] = b'V'; // C|G|T -> G|C|A
    rc[b'D' as usize] = b'H'; // A|G|T -> T|C|A
    rc[b'H' as usize] = b'D'; // A|C|T -> T|G|A
    rc[b'V' as usize] = b'B'; // A|C|G -> T|G|C
    rc[b'N' as usize] = b'N'; // A|C|G|T -> T|G|C|A
    rc[b'X' as usize] = b'X';
    // Lowercase versions
    rc[b'r' as usize] = b'y';
    rc[b'y' as usize] = b'r';
    rc[b's' as usize] = b's';
    rc[b'w' as usize] = b'w';
    rc[b'k' as usize] = b'm';
    rc[b'm' as usize] = b'k';
    rc[b'b' as usize] = b'v';
    rc[b'd' as usize] = b'h';
    rc[b'h' as usize] = b'd';
    rc[b'v' as usize] = b'b';
    rc[b'n' as usize] = b'n';
    rc[b'x' as usize] = b'x';
    rc
};

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
                BarcodeType::Ftag,
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
                BarcodeType::Rtag,
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

        let label_config = LabelConfig::new(true, true, true, true, None);
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, false, false);

        assert_eq!(results.len(), 1);
        let (trimmed_seq, trimmed_qual, group_label, _) = &results[0];
        println!("trimmed_seq: {}", String::from_utf8_lossy(trimmed_seq));
        assert_eq!(trimmed_seq, b"AAAA");
        assert_eq!(trimmed_qual, b"IIII");
        assert_eq!(group_label, "Fbar_fw__Rbar_fw");
    }

    #[test]
    fn test_two_cut_groups_produce_two_slices() {
        // seq indices: 0..8 C, 8..20 A, 20..26 C, 26..28 G, 28..30 C
        let seq = b"CCCCCCCCAAAAAAAAAAAACCCCCCGGCC";
        let qual = b"________IIIIIIIIIIII______II__";

        let read_len = seq.len();

        let annotations = vec![
            // Group 1: start at After(end_flank=8), end at Before(start_flank=20) -> slice 8..20
            BarbellMatch::new(
                4, // read_start_bar
                8, // read_end_bar
                4, // read_start_flank
                8, // read_end_flank
                0, // bar_start
                4, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "F1".to_string(),
                Strand::Fwd,
                read_len,
                "read1".to_string(),
                0,
                Some(vec![(Cut::new(1, CutDirection::After), 8)]),
            ),
            BarbellMatch::new(
                20, // read_start_bar
                24, // read_end_bar
                20, // read_start_flank
                24, // read_end_flank
                0,
                4,
                BarcodeType::Rtag,
                0,
                0,
                "R1".to_string(),
                Strand::Fwd,
                read_len,
                "read1".to_string(),
                0,
                Some(vec![(Cut::new(1, CutDirection::Before), 20)]),
            ),
            // Group 2: start at After(end_flank=26), end at Before(start_flank=28) -> slice 26..28
            BarbellMatch::new(
                24,
                26,
                24,
                26,
                0,
                2,
                BarcodeType::Ftag,
                0,
                0,
                "F2".to_string(),
                Strand::Fwd,
                read_len,
                "read1".to_string(),
                0,
                Some(vec![(Cut::new(2, CutDirection::After), 26)]),
            ),
            BarbellMatch::new(
                28,
                30,
                28,
                30,
                0,
                2,
                BarcodeType::Rtag,
                0,
                0,
                "R2".to_string(),
                Strand::Fwd,
                read_len,
                "read1".to_string(),
                0,
                Some(vec![(Cut::new(2, CutDirection::Before), 28)]),
            ),
        ];

        let label_config = LabelConfig::new(true, true, true, true, None);
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, false, false);

        assert_eq!(results.len(), 2);

        let (trimmed_seq1, trimmed_qual1, label1, _) = &results[0];
        assert_eq!(trimmed_seq1, b"AAAAAAAAAAAA");
        assert_eq!(trimmed_qual1, b"IIIIIIIIIIII");
        assert_eq!(label1, "F1_fw__R1_fw");

        let (trimmed_seq2, trimmed_qual2, label2, _) = &results[1];
        assert_eq!(trimmed_seq2, b"GG");
        assert_eq!(trimmed_qual2, b"II");
        assert_eq!(label2, "F2_fw__R2_fw");
    }

    #[test]
    fn test_trim_skipping() {
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
                BarcodeType::Ftag,
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
                BarcodeType::Rtag,
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

        let label_config = LabelConfig::new(true, true, true, true, None);
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, true, false);

        assert_eq!(results.len(), 1);
        let (trimmed_seq, trimmed_qual, group_label, _) = &results[0];
        println!("trimmed_seq: {}", String::from_utf8_lossy(trimmed_seq));
        assert_eq!(trimmed_seq, b"CCCCCCCCAAAACCCCCCCCCCCC");
        assert_eq!(trimmed_qual, b"________IIII____________");
        assert_eq!(group_label, "Fbar_fw__Rbar_fw");
    }

    #[test]
    fn test_flipping() {
        let seq = b"CCCCCCCCAGGCCCCCCCCCCCCC";
        let qual = b"________IIIA____________";

        let mut annotations = vec![
            BarbellMatch::new(
                4, // read_start_bar
                8, // read_end_bar
                4, // read_start_flank
                8, // read_end_flank
                0, // bar_start
                4, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "Fbar".to_string(),
                Strand::Rc, // Note RC match we have to flip
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
                BarcodeType::Rtag,
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

        let label_config = LabelConfig::new(true, true, true, true, None);
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, false, true);

        assert_eq!(results.len(), 1);
        let (trimmed_seq, trimmed_qual, group_label, _) = &results[0];
        println!("trimmed_seq: {}", String::from_utf8_lossy(trimmed_seq));
        assert_eq!(trimmed_seq, b"GCCT");
        assert_eq!(trimmed_qual, b"AIII");
        assert_eq!(group_label, "Fbar_rc__Rbar_fw");

        annotations[0].strand = Strand::Fwd;
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, false, true);
        let (trimmed_seq, trimmed_qual, group_label, _) = &results[0];
        println!("trimmed_seq: {}", String::from_utf8_lossy(trimmed_seq));
        assert_eq!(trimmed_seq, b"AGGC");
        assert_eq!(trimmed_qual, b"IIIA");
        assert_eq!(group_label, "Fbar_fw__Rbar_fw");

        // Chaning the strand in Fbar match should give original seq and qual again
    }
}
