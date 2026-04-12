use crate::annotate::barcodes::BarcodeType;
use crate::annotate::searcher::BarbellMatch;
use crate::config::{OutputFormat, TrimConfig};
use crate::filter::pattern::{Cut, CutDirection};
use crate::io::io::{open_reads, ReadSource};
use crate::progress::progress::{ProgressTracker, TRIM_PROGRESS_SPECS};
use anyhow::anyhow;
use csv;
use flate2::Compression;
use flate2::write::GzEncoder;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam as sam;
use noodles_sam::alignment::io::Write as _;
use sassy::Strand;
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

impl TrimConfig {
    pub fn get_label_config(&self) -> LabelConfig {
        LabelConfig::new(
            self.add_labels,
            self.add_orientation,
            self.add_flank,
            self.sort_labels,
            self.only_side,
        )
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

pub struct TrimmedRead {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub label: String,
    pub suffix: String,
    pub original_record: Option<bam::Record>,
    pub start: usize,
    pub end: usize,
}

pub fn process_read_and_anno(
    seq: &[u8],
    qual: &[u8],
    annotations: &[BarbellMatch],
    label_config: &LabelConfig,
    skip_trim: bool,
    flip: bool,
    original_bam_record: Option<&bam::Record>,
) -> Vec<TrimmedRead> {
    let mut results = Vec::new();
    let seq_len = seq.len();

    // Preprocess cuts to get complete slices
    let slices = preprocess_cuts(annotations, seq_len);

    // Group slices by cut group ID
    for (slice_count, slice) in slices.iter().enumerate() {
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
        let read_suffix = if slice_count == 0 {
            "".to_string()
        } else {
            format!("_{slice_count}")
        };

        let original_record = original_bam_record.cloned();

        results.push(TrimmedRead {
            seq: trimmed_seq,
            qual: trimmed_qual,
            label: group_label,
            suffix: read_suffix,
            original_record,
            start: slice.start,
            end: slice.end,
        });
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

fn build_record_buf_from_bam(
    original_record: &bam::Record,
    seq: &[u8],
    qual: &[u8],
    start: usize,
    end: usize,
    flip: bool,
) -> anyhow::Result<sam::alignment::record_buf::RecordBuf> {
    let mut builder = sam::alignment::RecordBuf::builder();

    if let Some(name) = original_record.name() {
        builder = builder.set_name(name.to_owned());
    }

    builder = builder.set_flags(original_record.flags());

    if let Some(reference_sequence_id) = original_record.reference_sequence_id().transpose()? {
        builder = builder.set_reference_sequence_id(reference_sequence_id);
    }

    if let Some(alignment_start) = original_record.alignment_start().transpose()? {
        builder = builder.set_alignment_start(alignment_start);
    }

    if let Some(mapping_quality) = original_record.mapping_quality() {
        builder = builder.set_mapping_quality(mapping_quality);
    }

    let cigar = sam::alignment::record_buf::Cigar::try_from(original_record.cigar())?;
    builder = builder.set_cigar(cigar);

    if let Some(mate_reference_sequence_id) = original_record.mate_reference_sequence_id().transpose()? {
        builder = builder.set_mate_reference_sequence_id(mate_reference_sequence_id);
    }

    if let Some(mate_alignment_start) = original_record.mate_alignment_start().transpose()? {
        builder = builder.set_mate_alignment_start(mate_alignment_start);
    }

    builder = builder.set_template_length(original_record.template_length());
    builder = builder.set_sequence(seq.to_vec().into());
    builder = builder.set_quality_scores(qual.to_vec().into());
    let mut data = sam::alignment::record_buf::Data::try_from(original_record.data())?;

    if start != 0 || end != original_record.sequence().len() || flip {
        update_mm_ml_tags(
            &mut data,
            original_record.sequence().as_ref(),
            seq,
            start,
            end,
            flip,
        )?;
    }

    builder = builder.set_data(data);

    Ok(builder.build())
}

fn update_mm_ml_tags(
    data: &mut sam::alignment::record_buf::Data,
    original_seq: &[u8],
    trimmed_seq: &[u8],
    start: usize,
    end: usize,
    flip: bool,
) -> anyhow::Result<()> {
    use sam::alignment::record::data::field::Tag;
    use sam::alignment::record_buf::data::field::Value;
    use sam::alignment::record_buf::data::field::value::Array;
    use sam::alignment::record_buf::Sequence;
    use sam::record::data::field::value::base_modifications::BaseModifications;
    use sam::record::data::field::value::base_modifications::group::{Group, Strand};

    if start == 0 && end == original_seq.len() && !flip {
        return Ok(());
    }

    if let Some(Value::String(mm_raw)) = data.get(&Tag::BASE_MODIFICATIONS) {
        let mm_bytes = mm_raw.as_ref();
        let original_sequence: Sequence = original_seq.to_vec().into();
        let mut base_modifications = BaseModifications::parse(mm_bytes, false, &original_sequence)?;

        let original_ml: Option<Vec<u8>> = match data.get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
            Some(Value::Array(Array::UInt8(values))) => {
                let mut bytes = Vec::new();
                for byte in values.iter() {
                    bytes.push(*byte);
                }
                Some(bytes)
            }
            _ => None,
        };

        let trimmed_len = trimmed_seq.len();
        let mut new_groups = Vec::new();
        let mut new_ml = Vec::new();
        let mut ml_index = 0;
        let has_ml = original_ml.is_some();

        for group in base_modifications.as_mut().iter() {
            let mut entries: Vec<(usize, Option<u8>)> = Vec::new();

            for &pos in group.positions() {
                if pos >= start && pos < end {
                    let ml_byte = original_ml.as_ref().map(|ml_bytes| ml_bytes[ml_index]);
                    entries.push((pos - start, ml_byte));
                }

                if has_ml {
                    ml_index += 1;
                }
            }

            if entries.is_empty() {
                continue;
            }

            let (positions, kept_ml) = if flip {
                let mut transformed: Vec<(usize, Option<u8>)> = entries
                    .into_iter()
                    .map(|(pos, ml)| (trimmed_len - 1 - pos, ml))
                    .collect();
                transformed.sort_unstable_by_key(|(pos, _)| *pos);
                let positions = transformed.iter().map(|(pos, _)| *pos).collect();
                let kept_ml = if has_ml {
                    transformed
                        .into_iter()
                        .map(|(_, ml)| ml.expect("ML length mismatch"))
                        .collect()
                } else {
                    Vec::new()
                };
                (positions, kept_ml)
            } else {
                let positions = entries.iter().map(|(pos, _)| *pos).collect();
                let kept_ml = if has_ml {
                    entries
                        .into_iter()
                        .map(|(_, ml)| ml.expect("ML length mismatch"))
                        .collect()
                } else {
                    Vec::new()
                };
                (positions, kept_ml)
            };

            let transformed_group = if flip {
                let reversed_strand = match group.strand() {
                    Strand::Forward => Strand::Reverse,
                    Strand::Reverse => Strand::Forward,
                };
                let complemented_base = group.unmodified_base().complement();
                Group::new(
                    complemented_base,
                    reversed_strand,
                    group.modifications().to_vec(),
                    group.status(),
                    positions,
                )
            } else {
                Group::new(
                    group.unmodified_base(),
                    group.strand(),
                    group.modifications().to_vec(),
                    group.status(),
                    positions,
                )
            };

            if original_ml.is_some() {
                new_ml.extend(kept_ml);
            }
            new_groups.push(transformed_group);
        }

        if new_groups.is_empty() {
            data.remove(&Tag::BASE_MODIFICATIONS);
            data.remove(&Tag::BASE_MODIFICATION_PROBABILITIES);
        } else {
            let new_mm = format_base_modifications(trimmed_seq, &new_groups)?;
            data.insert(Tag::BASE_MODIFICATIONS, Value::from(new_mm));

            if let Some(_) = original_ml {
                data.insert(Tag::BASE_MODIFICATION_PROBABILITIES, Value::from(new_ml));
            }
        }

        if data.get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH).is_some() {
            data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::from(trimmed_len as u32));
        }
    }

    Ok(())
}

fn format_base_modifications(
    sequence: &[u8],
    groups: &[sam::record::data::field::value::base_modifications::group::Group],
) -> anyhow::Result<String> {
    let mut output = String::new();

    for group in groups {
        let group_string = format_base_modifications_group(sequence, group)?;
        output.push_str(&group_string);
    }

    Ok(output)
}

fn format_base_modifications_group(
    sequence: &[u8],
    group: &sam::record::data::field::value::base_modifications::group::Group,
) -> anyhow::Result<String> {
    use sam::record::data::field::value::base_modifications::group::{Modification, Status, Strand};

    let base_char = u8::from(group.unmodified_base());
    let strand_char = match group.strand() {
        Strand::Forward => '+',
        Strand::Reverse => '-',
    };

    let modifications = group
        .modifications()
        .iter()
        .map(|modification| match modification {
            Modification::Code(c) => Ok((*c as char).to_string()),
            Modification::ChebiId(id) => Ok(id.to_string()),
        })
        .collect::<Result<Vec<_>, anyhow::Error>>()?
        .join("");

    let status_marker = match group.status() {
        Some(Status::Implicit) => ".",
        Some(Status::Explicit) => "?",
        None => "",
    };

    let counts = encode_mm_skip_counts(sequence, group.unmodified_base(), group.positions())?;
    let counts_string = counts
        .iter()
        .map(|count| count.to_string())
        .collect::<Vec<_>>()
        .join(",");

    let prefix = format!(
        "{}{}{}{}",
        base_char as char,
        strand_char,
        modifications,
        status_marker
    );

    Ok(format!("{prefix},{counts_string};"))
}

fn encode_mm_skip_counts(
    sequence: &[u8],
    unmodified_base: sam::record::data::field::value::base_modifications::group::UnmodifiedBase,
    positions: &[usize],
) -> anyhow::Result<Vec<usize>> {
    let base_byte = u8::from(unmodified_base);
    let matching_positions: Vec<usize> = sequence
        .iter()
        .enumerate()
        .filter_map(|(i, &b)| if b == base_byte { Some(i) } else { None })
        .collect();

    let mut counts = Vec::new();
    let mut previous_rank = 0;
    for (index, &pos) in positions.iter().enumerate() {
        let rank = matching_positions
            .iter()
            .position(|&match_pos| match_pos == pos)
            .ok_or_else(|| anyhow::anyhow!("MM group position {} not found", pos))?;

        if index == 0 {
            counts.push(rank);
        } else {
            counts.push(rank - previous_rank - 1);
        }
        previous_rank = rank;
    }

    Ok(counts)
}

fn should_flip(annotations: &[BarbellMatch]) -> bool {
    // If we matched an Ftag in rc we flip
    annotations
        .iter()
        .any(|anno| anno.match_type == BarcodeType::Ftag && anno.strand == Strand::Rc)
}

struct FastqTrimWriters {
    writers: HashMap<String, Box<dyn Write>>,
    output_format: OutputFormat,
    write_full_header: bool,
}

impl FastqTrimWriters {
    fn new(output_format: OutputFormat, write_full_header: bool) -> Self {
        Self {
            writers: HashMap::new(),
            output_format,
            write_full_header,
        }
    }

    fn write_record(
        &mut self,
        output_folder: &str,
        group: &str,
        read_id: &str,
        read_suffix: &str,
        desc: &str,
        trimmed_seq: &[u8],
        trimmed_qual: &[u8],
    ) -> anyhow::Result<()> {
        if !self.writers.contains_key(group) {
            let output_file = match self.output_format {
                OutputFormat::FastqGz => {
                    format!("{output_folder}/{group}.trimmed.fastq.gz")
                }
                _ => format!("{output_folder}/{group}.trimmed.fastq"),
            };

            let file = File::create(&output_file)
                .map_err(|err| anyhow!(format_output_file_error(&output_file, &err)))?;

            let writer: Box<dyn Write> = if self.output_format == OutputFormat::FastqGz {
                Box::new(GzEncoder::new(file, Compression::default()))
            } else {
                Box::new(BufWriter::new(file))
            };

            self.writers.insert(group.to_string(), writer);
        }

        let writer = self
            .writers
            .get_mut(group)
            .expect("writer should exist after insertion");

        if self.write_full_header {
            writeln!(writer, "@{read_id}{read_suffix} {desc}").expect("Failed to write header");
        } else {
            writeln!(writer, "@{read_id}{read_suffix}").expect("Failed to write header");
        }

        writeln!(writer, "{}", String::from_utf8_lossy(trimmed_seq))
            .expect("Failed to write sequence");
        writeln!(writer, "+").expect("Failed to write separator");
        writeln!(writer, "{}", String::from_utf8_lossy(trimmed_qual))
            .expect("Failed to write quality");

        Ok(())
    }

    fn finish(&mut self) -> anyhow::Result<()> {
        for (_, writer) in self.writers.iter_mut() {
            writer.flush().expect("Failed to flush output");
        }
        Ok(())
    }
}

struct BamTrimWriters {
    writers: HashMap<String, bam::io::Writer<bgzf::io::Writer<File>>>,
    header: sam::Header,
    flip: bool,
}

impl BamTrimWriters {
    fn new(header: sam::Header, flip: bool) -> Self {
        Self {
            writers: HashMap::new(),
            header,
            flip,
        }
    }

    fn write_record(
        &mut self,
        output_folder: &str,
        group: &str,
        original_record: Option<&bam::Record>,
        trimmed_seq: &[u8],
        trimmed_qual: &[u8],
        start: usize,
        end: usize,
    ) -> anyhow::Result<()> {
        if !self.writers.contains_key(group) {
            let output_file = format!("{output_folder}/{group}.trimmed.bam");
            let mut writer = bam::io::Writer::new(
                File::create(&output_file)
                    .map_err(|err| anyhow!(format_output_file_error(&output_file, &err)))?,
            );
            writer.write_header(&self.header).map_err(|err| {
                anyhow!(format!("Failed to write BAM header to '{output_file}': {err}"))
            })?;
            self.writers.insert(group.to_string(), writer);
        }

        let writer = self
            .writers
            .get_mut(group)
            .expect("BAM writer should exist after insertion");

        let original_record = original_record
            .as_ref()
            .ok_or_else(|| anyhow!("BAM output requires BAM input to preserve tags"))?;

        let record_buf = build_record_buf_from_bam(
            original_record,
            trimmed_seq,
            trimmed_qual,
            start,
            end,
            self.flip,
        )?;

        writer
            .write_alignment_record(&self.header, &record_buf)
            .map_err(|err| anyhow!(format!("Failed to write BAM record for group '{group}': {err}")))?;

        Ok(())
    }

    fn finish(&mut self) -> anyhow::Result<()> {
        for (_, writer) in self.writers.iter_mut() {
            writer
                .finish(&self.header)
                .expect("Failed to flush BAM output");
        }
        Ok(())
    }
}

enum TrimOutputWriter {
    Fastq(FastqTrimWriters),
    Bam(BamTrimWriters),
}

impl TrimOutputWriter {
    fn write_record(
        &mut self,
        output_folder: &str,
        group: &str,
        read_id: &str,
        read_suffix: &str,
        desc: &str,
        original_record: Option<&bam::Record>,
        trimmed_seq: &[u8],
        trimmed_qual: &[u8],
        start: usize,
        end: usize,
    ) -> anyhow::Result<()> {
        match self {
            Self::Fastq(writer) => writer.write_record(
                output_folder,
                group,
                read_id,
                read_suffix,
                desc,
                trimmed_seq,
                trimmed_qual,
            ),
            Self::Bam(writer) => writer.write_record(
                output_folder,
                group,
                original_record,
                trimmed_seq,
                trimmed_qual,
                start,
                end,
            ),
        }
    }

    fn finish(&mut self) -> anyhow::Result<()> {
        match self {
            Self::Fastq(writer) => writer.finish(),
            Self::Bam(writer) => writer.finish(),
        }
    }
}

pub fn trim_matches(
    filtered_match_file: &str,
    read_fastq_file: &str,
    output_folder: &str,
    config: &TrimConfig,
) -> anyhow::Result<()> {
    // Create output folder if it doesn't exist
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
    }

    // Label formatting config
    let label_config = config.get_label_config();

    if config.sort_labels && config.only_side.is_some() {
        return Err(anyhow!(
            "Cannot enable only keeping left/right label and sorting; this is ambiguous"
        ));
    }

    // Read all annotations and group by read ID
    let mut annotations_by_read: HashMap<String, Vec<BarbellMatch>> = HashMap::new();

    // Create progress bars
    let progress = if config.verbose {
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

    // If there is a failed trimmed writer, create it
    let mut failed_trimmed_writer = config
        .failed_trimmed_writer
        .as_ref()
        .map(|failed_trimmed_writer_path| {
            BufWriter::new(File::create(failed_trimmed_writer_path).unwrap())
        });

    // Process reads
    let mut reader = open_reads(read_fastq_file)?;
    let mut bam_header = None;
    if let ReadSource::Bam(ref mut bam_reader) = reader {
        bam_header = Some(
            bam_reader
                .read_header()
                .map_err(|err| anyhow!("Failed to read BAM header from '{read_fastq_file}': {err}"))?,
        );
    }

    if config.output_format == OutputFormat::Bam && bam_header.is_none() {
        return Err(anyhow!(
            "BAM output requires BAM input so that a BAM header and tags can be preserved"
        ));
    }

    let mut output_writer = match config.output_format {
        OutputFormat::Bam => TrimOutputWriter::Bam(BamTrimWriters::new(
            bam_header
                .as_ref()
                .expect("BAM header should exist for BAM output")
                .clone(),
            config.flip,
        )),
        format => TrimOutputWriter::Fastq(FastqTrimWriters::new(format, config.write_full_header)),
    };

    while let Some(record) = reader.next_record() {
        let record = record.expect("Error reading record");
        let read_id = record.id.clone();
        let desc = "";
        progress.inc(TOTAL_IDX);

        // Check if this read has annotations
        if let Some(annotations) = annotations_by_read.get(&read_id) {
            let results: Vec<TrimmedRead> = process_read_and_anno(
                &record.seq,
                &record.qual,
                annotations,
                &label_config,
                config.skip_trim,
                config.flip,
                record.original_bam_record.as_ref(),
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

            for trimmed in results {
                let TrimmedRead {
                    seq: trimmed_seq,
                    qual: trimmed_qual,
                    label: group,
                    suffix: read_suffix,
                    original_record,
                    start,
                    end,
                } = trimmed;

                if let Err(err) = output_writer.write_record(
                    output_folder,
                    &group,
                    &read_id,
                    &read_suffix,
                    desc,
                    original_record.as_ref(),
                    &trimmed_seq,
                    &trimmed_qual,
                    start,
                    end,
                ) {
                    progress.print_error(err.to_string());
                    progress.clear();
                    return Err(err);
                }
            }
        }

        progress.refresh();
    }

    output_writer.finish()?;

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
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, false, false, None);

        assert_eq!(results.len(), 1);
        let TrimmedRead { seq: trimmed_seq, qual: trimmed_qual, label: group_label, .. } = &results[0];
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
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, false, false, None);

        assert_eq!(results.len(), 2);

        let TrimmedRead { seq: trimmed_seq1, qual: trimmed_qual1, label: label1, .. } = &results[0];
        assert_eq!(trimmed_seq1, b"AAAAAAAAAAAA");
        assert_eq!(trimmed_qual1, b"IIIIIIIIIIII");
        assert_eq!(label1, "F1_fw__R1_fw");

        let TrimmedRead { seq: trimmed_seq2, qual: trimmed_qual2, label: label2, .. } = &results[1];
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
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, true, false, None);

        assert_eq!(results.len(), 1);
        let TrimmedRead { seq: trimmed_seq, qual: trimmed_qual, label: group_label, .. } = &results[0];
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
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, false, true, None);

        assert_eq!(results.len(), 1);
        let TrimmedRead { seq: trimmed_seq, qual: trimmed_qual, label: group_label, .. } = &results[0];
        println!("trimmed_seq: {}", String::from_utf8_lossy(trimmed_seq));
        assert_eq!(trimmed_seq, b"GCCT");
        assert_eq!(trimmed_qual, b"AIII");
        assert_eq!(group_label, "Fbar_rc__Rbar_fw");

        annotations[0].strand = Strand::Fwd;
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, false, true, None);
        let TrimmedRead { seq: trimmed_seq, qual: trimmed_qual, label: group_label, .. } = &results[0];
        println!("trimmed_seq: {}", String::from_utf8_lossy(trimmed_seq));
        assert_eq!(trimmed_seq, b"AGGC");
        assert_eq!(trimmed_qual, b"IIIA");
        assert_eq!(group_label, "Fbar_fw__Rbar_fw");

        // Changing the strand in Fbar match should give original seq and qual again
    }

    #[test]
    fn test_mm_ml_tags_are_trimmed() {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value;

        let mut data = sam::alignment::record_buf::Data::default();
        data.insert(Tag::BASE_MODIFICATIONS, Value::from("C+m,1,1;"));
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::from(vec![50u8, 80u8]),
        );
        data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::from(6u32));

        update_mm_ml_tags(&mut data, b"ACCCCC", b"CCCC", 2, 6, false)
            .expect("Failed to update MM/ML tags");

        let mm = data
            .get(&Tag::BASE_MODIFICATIONS)
            .expect("MM tag missing");
        assert_eq!(mm, &Value::from("C+m,0,1;"));

        let ml = data
            .get(&Tag::BASE_MODIFICATION_PROBABILITIES)
            .expect("ML tag missing");
        assert_eq!(ml, &Value::from(vec![50u8, 80u8]));

        let mn = data
            .get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH)
            .expect("MN tag missing");
        assert_eq!(mn, &Value::from(4u32));
    }

    #[test]
    fn test_mm_ml_tags_are_removed_when_no_modifications_survive() {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value;

        let mut data = sam::alignment::record_buf::Data::default();
        data.insert(Tag::BASE_MODIFICATIONS, Value::from("C+m,0,1;"));
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::from(vec![50u8, 80u8]),
        );
        data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::from(5u32));

        update_mm_ml_tags(&mut data, b"ACCCC", b"A", 0, 1, false)
            .expect("Failed to update MM/ML tags");

        assert!(data.get(&Tag::BASE_MODIFICATIONS).is_none());
        assert!(data.get(&Tag::BASE_MODIFICATION_PROBABILITIES).is_none());
        let mn = data
            .get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH)
            .expect("MN tag missing");
        assert_eq!(mn, &Value::from(1u32));
    }

    #[test]
    fn test_mm_ml_tags_support_multiple_modification_types() {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value;

        let mut data = sam::alignment::record_buf::Data::default();
        data.insert(Tag::BASE_MODIFICATIONS, Value::from("C+m,0,0;G+h,0;"));
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::from(vec![10u8, 20u8, 30u8]),
        );
        data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::from(5u32));

        update_mm_ml_tags(&mut data, b"ACCGT", b"CCG", 1, 4, false)
            .expect("Failed to update MM/ML tags");

        let mm = data
            .get(&Tag::BASE_MODIFICATIONS)
            .expect("MM tag missing");
        assert_eq!(mm, &Value::from("C+m,0,0;G+h,0;"));

        let ml = data
            .get(&Tag::BASE_MODIFICATION_PROBABILITIES)
            .expect("ML tag missing");
        assert_eq!(ml, &Value::from(vec![10u8, 20u8, 30u8]));

        let mn = data
            .get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH)
            .expect("MN tag missing");
        assert_eq!(mn, &Value::from(3u32));
    }

    #[test]
    fn test_mm_ml_tags_skip_trim_noop() {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value;

        let mut data = sam::alignment::record_buf::Data::default();
        data.insert(Tag::BASE_MODIFICATIONS, Value::from("C+m,1,1;"));
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::from(vec![50u8, 80u8]),
        );
        data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::from(6u32));

        update_mm_ml_tags(&mut data, b"ACCCCC", b"ACCCCC", 0, 6, false)
            .expect("Failed to update MM/ML tags");

        let mm = data
            .get(&Tag::BASE_MODIFICATIONS)
            .expect("MM tag missing");
        assert_eq!(mm, &Value::from("C+m,1,1;"));

        let ml = data
            .get(&Tag::BASE_MODIFICATION_PROBABILITIES)
            .expect("ML tag missing");
        assert_eq!(ml, &Value::from(vec![50u8, 80u8]));

        let mn = data
            .get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH)
            .expect("MN tag missing");
        assert_eq!(mn, &Value::from(6u32));
    }

    #[test]
    fn test_mm_ml_tags_preserve_status_sentinel() {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value;

        let mut data = sam::alignment::record_buf::Data::default();
        data.insert(Tag::BASE_MODIFICATIONS, Value::from("C+m?,1,1;"));
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::from(vec![50u8, 80u8]),
        );
        data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::from(6u32));

        update_mm_ml_tags(&mut data, b"ACCCCC", b"CCCC", 2, 6, false)
            .expect("Failed to update MM/ML tags");

        let mm = data
            .get(&Tag::BASE_MODIFICATIONS)
            .expect("MM tag missing");
        assert_eq!(mm, &Value::from("C+m?,0,1;"));
    }
}
