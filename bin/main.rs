use barbell::filter::filter::filter_from_text_file;
use barbell::inspect::inspect;
use barbell::progress::{ProgressTracker, print_header, print_summary_stats};
use barbell::search::barcodes::{BarcodeGroup, BarcodeType};
use barbell::search::searcher::Demuxer;
use barbell::trim::trim::trim_matches;
use clap::{Parser, Subcommand};
use colored::*;
use csv::WriterBuilder;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use needletail::{FastxReader, Sequence, parse_fastx_file};
use rand::Rng;
use seq_io::fastq::Reader as FastqReader;
use seq_io_parallel::{MinimalRefRecord, ParallelProcessor, ParallelReader};
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
    ) -> anyhow::Result<()> {
        let mut progress = ProgressTracker::new();
        print_header("Parallel Annotation");

        progress.step("Configuration");
        progress.indent();
        progress.substep(&format!("Input: {}", fastq_file));
        progress.substep(&format!("Output: {}", output_file));
        progress.substep(&format!("Threads: {}", threads));
        progress.dedent();

        // Create output file and initialize writer
        progress.step("Initializing output");
        progress.indent();
        let output_handle = std::fs::File::create(output_file)?;
        let writer = BufWriter::new(output_handle);
        let csv_writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(writer);
        self.writer = Arc::new(Mutex::new(csv_writer));
        progress.success("Output file initialized");
        progress.dedent();

        // Process reads
        progress.step("Processing reads");
        progress.indent();
        let reader = FastqReader::from_path(fastq_file)?;
        reader.process_parallel(self.clone(), threads)?;
        progress.dedent();

        // Finish progress bars
        self.total_bar.finish_with_message("Done!");
        self.found_bar.finish_with_message("Done!");
        self.missed_bar.finish_with_message("Done!");

        // Print summary
        let total = self.total.load(Ordering::Relaxed);
        let found = self.found.load(Ordering::Relaxed);
        let missed = self.missed.load(Ordering::Relaxed);

        print_summary_stats(
            total,
            found,
            found, // For annotation, mapped = found
            found, // For annotation, trimmed = found
            progress.elapsed(),
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
    ) -> anyhow::Result<()> {
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
            self.total_bar.set_message(total_count.to_string());
            self.found_bar
                .set_message(self.found.load(Ordering::Relaxed).to_string());
            self.missed_bar
                .set_message(self.missed.load(Ordering::Relaxed).to_string());
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> anyhow::Result<()> {
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
) -> anyhow::Result<()> {
    let mut progress = ProgressTracker::new();
    print_header("Barcode Annotation");

    progress.step("Initializing demuxer");
    progress.indent();
    // Set up a demuxer
    let mut demuxer = Demuxer::new(0.5);
    progress.substep("Created demuxer with alpha=0.5");
    progress.dedent();

    // Get random sample of reads
    progress.step("Sampling reads for tuning");
    progress.indent();
    let random_reads = random_fasta_sample(read_file, 10_000);
    progress.substep(&format!("Sampled {} reads for tuning", random_reads.len()));

    // Write random reads to fasta file
    progress.substep("Writing sampled reads to temporary file");
    let mut writer = BufWriter::new(File::create("random_reads.fasta").unwrap());
    for (i, seq) in random_reads.iter().enumerate() {
        writer.write_all(format!(">{}\n", i).as_bytes()).unwrap();
        writer.write_all(seq).unwrap();
        writer.write_all(b"\n").unwrap();
    }
    progress.success("Temporary FASTA file created");
    progress.dedent();

    // Create query groups from fasta files
    progress.step("Processing query groups");
    progress.indent();
    for (i, (query_file, query_type)) in query_files.iter().zip(query_types.iter()).enumerate() {
        progress.substep(&format!("Processing group {}: {}", i + 1, query_file));
        progress.indent();

        let mut query_group = BarcodeGroup::new_from_fasta(query_file, query_type.clone());
        progress.info(&format!("Loaded {} barcodes", query_group.barcodes.len()));

        progress.substep("Tuning barcode group");
        query_group.tune_group_random_sequences(1000, 0.001, 0.5, 15);
        //query_group.tune_group_with_fasta(&["random_reads.fasta"], 0.01, 0.5, 15, 10000);

        demuxer.add_query_group(query_group);
        progress.success(&format!("Group {} processed successfully", i + 1));
        progress.dedent();
    }
    progress.dedent();

    // Use parallel processor
    progress.step("Starting parallel annotation");
    progress.indent();
    let mut annotator = ParallelAnnotator::new(demuxer);
    progress.substep("Created parallel annotator");

    annotator.process_fastq(read_file, out_file, 10)?;
    progress.dedent();

    progress.success("Annotation completed successfully");
    progress.print_elapsed();

    Ok(())
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Annotate FASTQ files with barcode information
    Annotate {
        /// Input FASTQ file
        #[arg(short = 'i', long)]
        input: String,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 10)]
        threads: usize,

        /// Output file path
        #[arg(short = 'o', long, default_value = "output.tsv")]
        output: String,

        /// Query files (comma-separated paths)
        #[arg(short = 'q', long)]
        queries: String,

        /// Barcode types (comma-separated: Fbar,Rbar,Fflank,Rflank)
        #[arg(short = 'b', long, default_value = "Fbar")]
        barcode_types: String,
    },
    /// Filter annotation files based on pattern
    Filter {
        /// Input annotation file
        #[arg(short = 'i', long, required = true)]
        input: String,

        /// Output filtered file path
        #[arg(short = 'o', long, required = true)]
        output: String,

        /// File containing patterns to filter by
        #[arg(short = 'f', long, required = true)]
        file: String,
    },
    /// Trim and sort reads based on filtered annotations
    Trim {
        /// Input filtered annotation file
        #[arg(short = 'i', long)]
        input: String,

        /// Read FASTQ file
        #[arg(short = 'r', long)]
        reads: String,

        /// Output folder path for trimmed reads
        #[arg(short = 'o', long)]
        output: String,

        /// Disable label in output filenames
        #[arg(long, default_value_t = false)]
        no_label: bool,

        /// Disable orientation in output filenames
        #[arg(long, default_value_t = false)]
        no_orientation: bool,

        /// Disable flank in output filenames
        #[arg(long, default_value_t = false)]
        no_flanks: bool,

        /// Sort barcode labels in output filenames
        #[arg(long, default_value_t = false)]
        sort_labels: bool,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 1)]
        threads: usize,
    },

    /// View most common patterns in annotation
    Inspect {
        /// Input filtered annotation file
        #[arg(short = 'i', long)]
        input: String,

        /// Top N
        #[arg(short = 'n', long, default_value_t = 10)]
        top_n: usize,
    },
}

fn main() {
    print_banner();

    let cli = Cli::parse();

    match &cli.command {
        Commands::Annotate {
            input,
            threads,
            output,
            queries,
            barcode_types,
        } => {
            println!("{}", "Starting annotation...".green());

            // Split comma-separated query paths into Vec<String>
            let query_files: Vec<String> =
                queries.split(',').map(|s| s.trim().to_string()).collect();

            let query_files_refs: Vec<&str> = query_files.iter().map(|s| s.as_str()).collect();

            // Parse barcode types
            let barcode_types_vec: Vec<BarcodeType> = barcode_types
                .split(',')
                .map(|s| match s.trim() {
                    "Fbar" => BarcodeType::Fbar,
                    "Rbar" => BarcodeType::Rbar,
                    "Fflank" => BarcodeType::Fflank,
                    "Rflank" => BarcodeType::Rflank,
                    _ => panic!("Unknown barcode type: {}", s),
                })
                .collect();

            match annotate(input, query_files_refs, barcode_types_vec, output) {
                Ok(_) => println!("{}", "Annotation complete!".green()),
                Err(e) => println!("{} {}", "Error during processing:".red(), e),
            }
        }

        Commands::Filter {
            input,
            output,
            file,
        } => {
            println!("{}", "Starting filtering...".green());

            match filter_from_text_file(input, file, output) {
                Ok(_) => println!("{}", "Filtering successful!".green()),
                Err(e) => println!("{} {}", "Filtering failed:".red(), e),
            }
        }

        Commands::Trim {
            input,
            reads,
            output,
            no_label,
            no_orientation,
            no_flanks,
            sort_labels,
            threads,
        } => {
            println!("{}", "Starting trimming...".green());
            trim_matches(
                input,
                reads,
                output,
                !no_label,
                !no_orientation,
                !no_flanks,
                *sort_labels,
                *threads,
            );
        }

        Commands::Inspect { input, top_n } => {
            println!("{}", "Inspecting...".green());

            match inspect::inspect(input, *top_n) {
                Ok(_) => println!("{}", "Inspection complete!".green()),
                Err(e) => println!("{} {}", "Inspection failed:".red(), e),
            }
        }
    }
}

fn print_banner() {
    println!(
        "{}",
        r#"
    ██████╗  █████╗ ██████╗ ██████╗ ███████╗██╗     ██╗     
    ██╔══██╗██╔══██╗██╔══██╗██╔══██╗██╔════╝██║     ██║     
    ██████╔╝███████║██████╔╝██████╔╝█████╗  ██║     ██║     
    ██╔══██╗██╔══██║██╔══██╗██╔══██╗██╔══╝  ██║     ██║     
    ██████╔╝██║  ██║██║  ██║██████╔╝███████╗███████╗███████╗
    ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═════╝ ╚══════╝╚══════╝╚══════╝
    "#
        .blue()
    );
    println!(
        "{}",
        "        [===]------------------------------------------[===]        ".bright_yellow()
    );
}
