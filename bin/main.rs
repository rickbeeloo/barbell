use anyhow::Result;
use barbell::filter::filter::filter_from_text_file;
use barbell::inspect::inspect;
use barbell::search::barcodes::{BarcodeGroup, BarcodeType};
use barbell::search::searcher::Demuxer;
use colored::Colorize;
use csv::WriterBuilder;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use needletail::{FastxReader, Sequence, parse_fastx_file};
use rand::Rng;
use seq_io::fastq::Reader as FastqReader;
use seq_io_parallel::{MinimalRefRecord, ParallelProcessor, ParallelReader};
use serde_json;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread_local;
use std::time::Instant;
use std::{
    path::PathBuf,
    sync::{Arc, Mutex},
};

thread_local! {
    static DEMUXER: std::cell::RefCell<Option<Demuxer>> = std::cell::RefCell::new(None);
}

#[derive(Clone)]
pub struct ParallelAnnotator {
    demuxer: Arc<Demuxer>,
    writer: Arc<Mutex<csv::Writer<BufWriter<std::fs::File>>>>,
    total: Arc<AtomicUsize>,
    found: Arc<AtomicUsize>,
    missed: Arc<AtomicUsize>,
    total_bar: Arc<ProgressBar>,
    found_bar: Arc<ProgressBar>,
    missed_bar: Arc<ProgressBar>,
}

impl ParallelAnnotator {
    pub fn new(demuxer: Demuxer) -> Self {
        let multi_progress = MultiProgress::new();
        let total_bar = multi_progress.add(ProgressBar::new_spinner());
        let found_bar = multi_progress.add(ProgressBar::new_spinner());
        let missed_bar = multi_progress.add(ProgressBar::new_spinner());

        // Style the progress bars
        total_bar.set_style(
            ProgressStyle::with_template("{spinner:.blue} {prefix:<12} {msg:>6} {elapsed_precise}")
                .unwrap()
                .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
        );
        found_bar.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} {prefix:<12} {msg:>6} {elapsed_precise}",
            )
            .unwrap()
            .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
        );
        missed_bar.set_style(
            ProgressStyle::with_template("{spinner:.red} {prefix:<12} {msg:>6} {elapsed_precise}")
                .unwrap()
                .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
        );

        total_bar.set_prefix("Total:");
        found_bar.set_prefix("Found:");
        missed_bar.set_prefix("Missed:");

        Self {
            demuxer: Arc::new(demuxer),
            writer: Arc::new(Mutex::new(csv::Writer::from_writer(BufWriter::new(
                std::fs::File::create("/dev/null").unwrap(),
            )))),
            total: Arc::new(AtomicUsize::new(0)),
            found: Arc::new(AtomicUsize::new(0)),
            missed: Arc::new(AtomicUsize::new(0)),
            total_bar: Arc::new(total_bar),
            found_bar: Arc::new(found_bar),
            missed_bar: Arc::new(missed_bar),
        }
    }

    pub fn process_fastq(
        &mut self,
        fastq_file: &str,
        output_file: &str,
        threads: usize,
    ) -> Result<()> {
        let start_time = Instant::now();

        println!("\n{}", "Configuration".bold().underline());
        println!("  • Input:  {}", fastq_file.bold());
        println!("  • Output: {}", output_file.bold());
        println!("  • Threads: {}", threads.to_string().bold());

        // Create output file and initialize writer
        let output_handle = std::fs::File::create(output_file)?;
        let writer = BufWriter::new(output_handle);
        let csv_writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(writer);
        self.writer = Arc::new(Mutex::new(csv_writer));

        // Process reads
        let reader = FastqReader::from_path(fastq_file)?;
        reader.process_parallel(self.clone(), threads)?;

        // Finish progress bars
        self.total_bar.finish_with_message("Done!");
        self.found_bar.finish_with_message("Done!");
        self.missed_bar.finish_with_message("Done!");

        println!("\n{}", "Summary".bold().underline());
        println!(
            "  • Time: {} seconds",
            start_time.elapsed().as_secs().to_string().bold()
        );
        println!(
            "  • Total reads: {}",
            self.total.load(Ordering::Relaxed).to_string().bold()
        );
        println!(
            "  • Tagged reads: {}",
            self.found
                .load(Ordering::Relaxed)
                .to_string()
                .green()
                .bold()
        );
        println!(
            "  • Missed reads: {}\n",
            self.missed.load(Ordering::Relaxed).to_string().red().bold()
        );

        // Flush before closing
        self.writer.lock().unwrap().flush()?;

        Ok(())
    }
}

impl ParallelProcessor for ParallelAnnotator {
    fn process_record<'a, Rf: MinimalRefRecord<'a>>(
        &mut self,
        record: Rf,
        _record_set_idx: usize,
        _record_idx: usize,
    ) -> Result<()> {
        // Parse the read ID and sequence
        let read_id = record.ref_id().unwrap().split_whitespace().next().unwrap();
        let read = record.ref_seq();

        // Initialize demuxer for this thread if not already done
        DEMUXER.with(|cell| {
            if cell.borrow().is_none() {
                *cell.borrow_mut() = Some(self.demuxer.as_ref().clone());
            }
        });

        // Annotate the read
        let matches = DEMUXER.with(|cell| {
            if let Some(ref mut demuxer) = *cell.borrow_mut() {
                demuxer.demux(read_id, read)
            } else {
                vec![]
            }
        });

        // Update total counter
        let total_count = self.total.fetch_add(1, Ordering::Relaxed);

        // Get the writer lock
        let mut writer = self.writer.lock().unwrap();

        if !matches.is_empty() {
            // Update found counter
            self.found.fetch_add(1, Ordering::Relaxed);

            // Write matches using CSV writer (which uses serde)
            for m in &matches {
                writer.serialize(m)?;
            }
        } else {
            self.missed.fetch_add(1, Ordering::Relaxed);
        }

        // Release the lock when done
        drop(writer);

        // Update progress bars periodically
        if total_count % 100 == 0 {
            self.total_bar.set_message(format!("{}", total_count));
            self.found_bar
                .set_message(format!("{}", self.found.load(Ordering::Relaxed)));
            self.missed_bar
                .set_message(format!("{}", self.missed.load(Ordering::Relaxed)));
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<()> {
        // Optionally flush after each batch
        // self.writer.lock().unwrap().flush()?;
        Ok(())
    }
}

fn random_fasta_sample(fasta_file: &str, n_reads: usize) -> Vec<Vec<u8>> {
    let mut reader = parse_fastx_file(fasta_file).unwrap();
    let mut rng = rand::thread_rng();
    let mut sample = Vec::new();
    let mut count = 0;

    while let Some(Ok(record)) = reader.next() {
        count += 1;
        let seq = record.seq().to_vec();

        if count <= n_reads {
            // Fill the sample initially
            sample.push(seq);
        } else {
            // Reservoir sampling: replace with probability n_reads/count
            let j = rng.gen_range(0..count);
            if j < n_reads {
                sample[j] = seq;
            }
        }
    }

    println!("Sampled {} reads from {} total reads", sample.len(), count);
    sample
}

fn annotate(
    read_file: &str,
    query_files: Vec<&str>,
    query_types: Vec<BarcodeType>,
    out_file: &str,
) {
    // Set up a demuxer
    let mut demuxer = Demuxer::new(0.5);

    // Get random sample of reads
    let random_reads = random_fasta_sample(read_file, 10_000);

    // Write random reads to fasta file
    let mut writer = BufWriter::new(File::create("random_reads.fasta").unwrap());
    for (i, seq) in random_reads.iter().enumerate() {
        writer.write_all(format!(">{}\n", i).as_bytes()).unwrap();
        writer.write_all(seq).unwrap();
        writer.write_all(b"\n").unwrap();
    }

    // Create query groups from fasta files
    for (query_file, query_type) in query_files.iter().zip(query_types.iter()) {
        let mut query_group = BarcodeGroup::new_from_fasta(query_file, query_type.clone());
        // Tune it
        //query_group.tune_group(100_000, 0.01, 0.5, 5); // 0.01% false positives
        query_group.tune_group_with_fasta(&["random_reads.fasta"], 0.001, 0.5, 15); // 0.1% false positives

        //let barcodes = query_group.barcodes.len();
        //let barcode_cut_offs = vec![5; barcodes];
        //query_group.tune_group_manual(41, barcode_cut_offs);
        demuxer.add_query_group(query_group);
    }

    // Use parallel processor
    let mut annotator = ParallelAnnotator::new(demuxer);
    annotator.process_fastq(read_file, out_file, 10).unwrap();
}

fn main() {
    annotate(
        "/home/solprof/Downloads/sub.fastq",
        vec!["/home/solprof/PhD/barbell-sassy-rewrite/examples/rapid_bars.fasta"],
        vec![BarcodeType::Fbar],
        "examples/test.tsv",
    );
    inspect::inspect("examples/test.tsv", 10).expect("Failed to parse results");

    filter_from_text_file(
        "examples/test.tsv",
        "filter.txt",
        "examples/test_filtered.tsv",
    )
    .expect("Failed to filter");
}
