use barbell::annotate::annotate::*;
use barbell::annotate::barcodes::{BarcodeGroup, BarcodeType};
use barbell::annotate::searcher::Demuxer;
use barbell::filter::filter::filter_from_text_file;
use barbell::inspect::inspect;
use barbell::progress::{ProgressTracker, print_header, print_summary_stats};
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

        /// Max errors, percentage of length
        #[arg(short = 'e', long, default_value = "0.1")]
        max_error_perc: f32,
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
            max_error_perc,
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

            match annotate(
                input,
                query_files_refs,
                barcode_types_vec,
                output,
                *max_error_perc,
            ) {
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
