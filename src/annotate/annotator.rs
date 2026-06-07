use crate::annotate::barcodes::{BarcodeGroup, BarcodeType};
use crate::annotate::edit_model::get_edit_cut_off;
use crate::annotate::searcher::{BarbellMatch, Demuxer};
use crate::config::AnnotateConfig;
use crate::io::io::{open_fastq_collection, split_fastq_header};
use crate::progress::progress::{ProgressTracker, ANNOTATION_PROGRESS_SPECS};
use anyhow::anyhow;
use paraseq::prelude::{ParallelProcessor, Record};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

#[inline(always)]
fn write_annotation_batch(
    writer: &Arc<Mutex<csv::Writer<std::fs::File>>>,
    record_set_results: &[BarbellMatch],
) -> anyhow::Result<()> {
    let mut writer = writer
        .lock()
        .map_err(|err| anyhow!("Annotation writer lock failed: {err}"))?;
    for annotation in record_set_results {
        writer
            .serialize(annotation)
            .map_err(|err| anyhow!("Failed to write annotation output: {err}"))?;
    }
    Ok(())
}

struct DemuxProcessor {
    query_groups: Vec<BarcodeGroup>,
    alpha: f32,
    verbose: bool,
    min_score: f64,
    min_score_diff: f64,
    writer: Arc<Mutex<csv::Writer<std::fs::File>>>,
    progress: Arc<ProgressTracker>,
    demuxer: Option<Demuxer>,
    annotations: Vec<BarbellMatch>,
    total_count: usize,
    found_count: usize,
    thread_id: usize,
}

impl Clone for DemuxProcessor {
    fn clone(&self) -> Self {
        Self {
            query_groups: self.query_groups.clone(),
            alpha: self.alpha,
            verbose: self.verbose,
            min_score: self.min_score,
            min_score_diff: self.min_score_diff,
            writer: Arc::clone(&self.writer),
            progress: Arc::clone(&self.progress),
            demuxer: None,
            annotations: Vec::new(),
            total_count: 0,
            found_count: 0,
            thread_id: self.thread_id,
        }
    }
}

impl DemuxProcessor {
    fn new(
        query_groups: Vec<BarcodeGroup>,
        alpha: f32,
        verbose: bool,
        min_score: f64,
        min_score_diff: f64,
        writer: Arc<Mutex<csv::Writer<std::fs::File>>>,
        progress: Arc<ProgressTracker>,
    ) -> Self {
        Self {
            query_groups,
            alpha,
            verbose,
            min_score,
            min_score_diff,
            writer,
            progress,
            demuxer: None,
            annotations: Vec::new(),
            total_count: 0,
            found_count: 0,
            thread_id: 0,
        }
    }

    fn demuxer(&mut self) -> &mut Demuxer {
        self.demuxer.get_or_insert_with(|| {
            let mut demux = Demuxer::new(
                self.alpha,
                self.verbose,
                self.min_score,
                self.min_score_diff,
            );
            for query_group in &self.query_groups {
                demux.add_query_group(query_group.clone());
            }
            demux
        })
    }

    fn flush_batch(&mut self) -> anyhow::Result<()> {
        if !self.annotations.is_empty() {
            write_annotation_batch(&self.writer, &self.annotations)?;
            self.annotations.clear();
        }

        if self.total_count > 0 {
            self.progress.add(0, self.total_count);
            self.progress.add(1, self.found_count);
            self.progress.add(2, self.total_count - self.found_count);
            self.total_count = 0;
            self.found_count = 0;
        }

        self.progress.refresh();
        Ok(())
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for DemuxProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        self.total_count += 1;
        let (read_id, _) = split_fastq_header(record.id())?;
        let seq = record.seq();
        let matches = self.demuxer().demux(read_id, seq.as_ref());

        if !matches.is_empty() {
            self.found_count += 1;
            self.annotations.extend(matches);
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        self.flush_batch().map_err(Into::into)
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        self.flush_batch().map_err(Into::into)
    }

    fn set_thread_id(&mut self, thread_id: usize) {
        self.thread_id = thread_id;
    }

    fn get_thread_id(&self) -> usize {
        self.thread_id
    }
}

// used by custom experiments (direct annotate call)
pub fn annotate_with_files(
    read_files: &[PathBuf],
    query_files: &[PathBuf],
    query_types: Vec<BarcodeType>,
    out_file: &str,
    config: &AnnotateConfig,
) -> anyhow::Result<()> {
    if query_files.len() != query_types.len() {
        return Err(anyhow!(
            "Expected the same number of query files and barcode types, got {} query file(s) and {} barcode type(s)",
            query_files.len(),
            query_types.len()
        ));
    }

    // Get query groups
    let mut query_groups = Vec::new();
    for (query_file, query_type) in query_files.iter().zip(query_types.iter()) {
        let mut query_group = BarcodeGroup::new_from_fasta(
            query_file.to_str().ok_or_else(|| {
                anyhow!(
                    "Query FASTA path is not valid UTF-8: {}",
                    query_file.display()
                )
            })?,
            query_type.clone(),
        );
        if let Some(max_flank_errors) = config.max_flank_errors {
            query_group.set_flank_threshold(max_flank_errors);
        } else {
            // Determine based on formula
            let edit_cut_off = get_edit_cut_off(query_group.get_effective_len());
            println!("Auto edit flank cut off: {edit_cut_off}");
            query_group.set_flank_threshold(edit_cut_off);
        }
        query_groups.push(query_group);
    }
    annotate(read_files, out_file, query_groups, config)
}

// we could maybe just discard annotate_with_groups and only have kit or fasta
pub fn annotate_with_kit(
    read_files: &[PathBuf],
    out_file: &str,
    kit: &str,
    config: &AnnotateConfig,
) -> anyhow::Result<()> {
    let query_groups: Vec<BarcodeGroup> = BarcodeGroup::new_from_kit(kit, config.use_extended);
    annotate_with_groups(read_files, out_file, query_groups, config)
}

// used by kit
pub fn annotate_with_groups(
    read_files: &[PathBuf],
    out_file: &str,
    query_groups: Vec<BarcodeGroup>,
    config: &AnnotateConfig,
) -> anyhow::Result<()> {
    // Hmm not sure: fixme: think about where flank error should be set
    // Cannot mutate query_groups because it's not mutable (Vec<BarcodeGroup> is not mutable when passed by value and iterated by reference).
    // Instead, create a new Vec with updated groups.
    let query_groups: Vec<BarcodeGroup> = query_groups
        .into_iter()
        .map(|mut query_group| {
            if let Some(max_flank_errors) = config.max_flank_errors {
                query_group.set_flank_threshold(max_flank_errors);
            } else {
                // Determine based on formula
                let edit_cut_off = get_edit_cut_off(query_group.get_effective_len());
                println!("Auto edit flank cut off: {edit_cut_off}");
                query_group.set_flank_threshold(edit_cut_off);
            }
            query_group
        })
        .collect();
    annotate(read_files, out_file, query_groups, config)
}

pub fn annotate(
    read_files: &[PathBuf],
    out_file: &str,
    query_groups: Vec<BarcodeGroup>,
    config: &AnnotateConfig,
) -> anyhow::Result<()> {
    let alpha = config.alpha;
    let n_threads = config.n_threads;
    let verbose = config.verbose;
    let min_score = config.min_score;
    let min_score_diff = config.min_score_diff;

    let reader = open_fastq_collection(read_files)?;
    let writer = Arc::new(Mutex::new(
        csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(out_file)
            .map_err(|e| anyhow!("Failed to create annotation output file '{out_file}': {e}"))?,
    ));

    // Dispaly to user
    for (i, query_group) in query_groups.iter().enumerate() {
        println!("{}: {}", query_group.barcode_type.as_str(), i);
        query_group.display(5);
    }

    let progress = Arc::new(if config.verbose {
        let log_dir = Path::new(out_file)
            .parent()
            .unwrap_or_else(|| Path::new("."));
        ProgressTracker::new_with_logging(&ANNOTATION_PROGRESS_SPECS, "annotate", log_dir)
    } else {
        ProgressTracker::new(&ANNOTATION_PROGRESS_SPECS)
    });

    let mut processor = DemuxProcessor::new(
        query_groups,
        alpha,
        verbose,
        min_score,
        min_score_diff,
        writer,
        Arc::clone(&progress),
    );

    reader
        .process_parallel(&mut processor, n_threads as usize, None)
        .map_err(|err| anyhow!("Input FASTQ parsing failed: {err}"))?;

    progress.finish("records");

    Ok(())
}
