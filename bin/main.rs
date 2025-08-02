use barbell::annotate::annotator::*;
use barbell::annotate::barcodes::{BarcodeGroup, BarcodeType};
use barbell::annotate::searcher::Demuxer;
use barbell::filter::filter::filter_from_text_file;
use barbell::inspect::inspect;
use barbell::preset::presets::PresetName;
use barbell::preset::presets::use_preset;
use barbell::trim::trim::trim_matches;
use barbell::tune::tune::{TargetSide, tune};
use clap::{Parser, Subcommand};
use colored::*;

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

        /// Maximum errors as percentage of flank / barcode length (mutually exclusive with specific error counts)
        #[arg(short = 'e', long, conflicts_with_all = ["flank_max_errors", "barcode_max_errors"])]
        max_error_perc: Option<f32>,

        /// Maximum absolute errors allowed in flanks (mutually exclusive with -e/--max-error-perc)
        #[arg(
            long = "flank-max-errors",
            value_name = "INT",
            conflicts_with = "max_error_perc"
        )]
        flank_max_errors: Option<usize>,

        /// Maximum absolute errors allowed in barcodes (mutually exclusive with -e/--max-error-perc)
        #[arg(
            long = "barcode-max-errors",
            value_name = "INT",
            conflicts_with = "max_error_perc"
        )]
        barcode_max_errors: Option<usize>,
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

        /// Write dropped read annotation to this file
        #[arg(long)]
        dropped: Option<String>,
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
    },

    /// View most common patterns in annotation
    Inspect {
        /// Input filtered annotation file
        #[arg(short = 'i', long)]
        input: String,

        /// Top N
        #[arg(short = 'n', long, default_value_t = 10)]
        top_n: usize,

        /// Write pattern for each read to this file (optional)
        #[arg(short = 'o', long)]
        read_pattern_out: Option<String>,
    },

    /// Run a preset
    Preset {
        /// Preset to use
        #[arg(short = 'p', long)]
        preset: PresetName,

        /// Input FASTQ file
        #[arg(short = 'i', long)]
        input: String,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 10)]
        threads: usize,

        /// Output folder
        #[arg(short = 'o', long)]
        output: String,

        /// Maximum errors allowed in flanks
        #[arg(short = 'f', long, default_value_t = 26)]
        flank_max_errors: usize,

        /// Maximum errors allowed in barcodes
        #[arg(short = 'b', long, default_value_t = 6)]
        barcode_max_errors: usize,
    },

    /// Tune the parameters for a given query file
    Tune {
        /// Input FASTQ file
        #[arg(short = 'i', long)]
        input: String,

        /// Query file
        #[arg(short = 'q', long)]
        query: String,

        /// Target side
        #[arg(short = 't', long)]
        target_side: TargetSide,
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
            flank_max_errors,
            barcode_max_errors,
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
                    "Ftag" => BarcodeType::Ftag,
                    "Rtag" => BarcodeType::Rtag,
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
                max_error_perc.clone(),
                flank_max_errors.clone(),
                barcode_max_errors.clone(),
                0.5,
                *threads as u32,
            ) {
                Ok(_) => println!("{}", "Annotation complete!".green()),
                Err(e) => println!("{} {}", "Error during processing:".red(), e),
            }
        }

        Commands::Filter {
            input,
            output,
            file,
            dropped,
        } => {
            println!("{}", "Starting filtering...".green());

            match filter_from_text_file(input, file, output, dropped.as_deref()) {
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
            );
        }

        Commands::Inspect {
            input,
            top_n,
            read_pattern_out,
        } => {
            println!("{}", "Inspecting...".green());

            match inspect::inspect(input, *top_n, read_pattern_out.clone()) {
                Ok(_) => println!("{}", "Inspection complete!".green()),
                Err(e) => println!("{} {}", "Inspection failed:".red(), e),
            }
        }

        Commands::Preset {
            preset,
            input,
            threads,
            output,
            flank_max_errors,
            barcode_max_errors,
        } => {
            use_preset(
                preset.clone(),
                input,
                *threads,
                output,
                *flank_max_errors,
                *barcode_max_errors,
            );
        }

        Commands::Tune {
            input,
            query,
            target_side,
        } => {
            tune(input, query, target_side.clone());
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
