use barbell::annotate::annotator::*;
use barbell::annotate::barcodes::BarcodeType;
use barbell::filter::filter::filter_from_text_file;
use barbell::inspect::inspect;
use barbell::kits::use_kit::demux_using_kit;
use barbell::trim::trim::{LabelSide, trim_matches};
use clap::{Parser, Subcommand};
use colored::*;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum UtilsCommands {
    /// Group raw reads by label without trimming
    Pull {
        /// Input filtered annotation file
        #[arg(short = 'i', long)]
        input: String,

        /// Read FASTQ file
        #[arg(short = 'r', long)]
        reads: String,

        /// Output folder path for grouped reads
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

        /// Only keep left or right label in output filenames
        #[arg(long, conflicts_with = "sort_labels")]
        only_side: Option<LabelSide>,

        /// Write ids of reads with no annotations to this file
        #[arg(long)]
        failed_out: Option<String>,
    },
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
        #[arg(short = 'q', long, required_unless_present = "kit")]
        queries: Option<String>,

        /// Barcode types (comma-separated: Ftag,Rtag) matching your query file (-q)
        #[arg(short = 'b', long, default_value = "Ftag")]
        barcode_types: String,

        /// Kit name (e.g. SQK-RBK114-24). Conflicts with --queries/--barcode-types
        #[arg(long, conflicts_with = "queries", conflicts_with = "barcode_types")]
        kit: Option<String>,

        /// Flank maximum erors in flank, ONLY set manually when you know what you are doing
        #[arg(long = "flank-max-errors", value_name = "INT")]
        flank_max_errors: Option<usize>,

        /// Enable verbose output for debugging
        #[arg(long, default_value_t = false)]
        verbose: bool,

        /// Barcode: fraction compared to 'perfect' match score for top candidate
        #[arg(long = "min-score", default_value_t = 0.2)]
        min_score: f64,

        /// Barcode: fraction difference between top 2 candidates
        #[arg(long = "min-score-diff", default_value_t = 0.1)]
        min_score_diff: f64,

        /// Also use extended templates (if using kit), i.e. detect fusions, breaks, etc. (slower)
        #[arg(long, default_value_t = false)]
        use_extended: bool,

        /// Edit cost beyond read boundaries
        #[arg(long = "alpha", default_value_t = 0.4)]
        alpha: f32,
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

        /// Only keep left or right label in output filenames
        #[arg(long, conflicts_with = "sort_labels")]
        only_side: Option<LabelSide>,

        /// Write ids of failed trimmed reads to this file
        #[arg(long)]
        failed_out: Option<String>,
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

        /// To summarize results we uses "buckets", such that matches 100 and 103 from the start end up in the same bucket
        #[arg(short = 's', long = "bucket-size", default_value_t = 250)]
        bucket_size: usize,
    },

    /// Run a preset
    Kit {
        /// Kit to use (e.g. SQK-RBK114-24)
        #[arg(short = 'k', long)]
        kit: String,

        /// Input FASTQ file
        #[arg(short = 'i', long)]
        input: String,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 10)]
        threads: usize,

        /// Output folder
        #[arg(short = 'o', long)]
        output: String,

        /// Add more 'risky' patterns to demuxing to maximize assigned reads
        #[arg(long, default_value_t = false)]
        maximize: bool,

        /// Enable verbose output for debugging
        #[arg(long, default_value_t = false)]
        verbose: bool,

        /// Fraction compared to 'perfect' match score for top candidate
        #[arg(long = "min-score", default_value_t = 0.2)]
        min_score: f64,

        /// Fraction difference between top 2 candidates
        #[arg(long = "min-score-diff", default_value_t = 0.1)]
        min_score_diff: f64,

        /// Flank maximum erors in flank, ONLY set manually when you know what you are doing
        #[arg(long = "flank-max-errors", value_name = "INT")]
        flank_max_errors: Option<usize>,

        /// Write ids of failed trimmed reads to this file
        #[arg(long)]
        failed_out: Option<String>,

        /// Also use extended templates (if using kit), i.e. detect fusions, breaks, etc. (slower)
        #[arg(long, default_value_t = false)]
        use_extended: bool,

        /// Edit cost beyond read boundaries
        #[arg(long = "alpha", default_value_t = 0.4)]
        alpha: f32,
    },
}

fn main() {
    // Make sure that AVX2 is supported by the current CPU if it was used during compilation.
    ensure_simd::ensure_simd();

    print_banner();

    let cli = Cli::parse();

    match &cli.command {
        Commands::Annotate {
            input,
            threads,
            output,
            queries,
            barcode_types,
            kit,
            flank_max_errors,
            verbose,
            min_score,
            min_score_diff,
            use_extended,
            alpha,
        } => {
            println!("{}", "Starting annotation...".green());

            if let Some(kit_name) = kit.as_ref() {
                match annotate_with_kit(
                    input,
                    output,
                    kit_name.as_str(),
                    *flank_max_errors,
                    *alpha,
                    *threads as u32,
                    *verbose,
                    *min_score,
                    *min_score_diff,
                    *use_extended,
                ) {
                    Ok(_) => println!("{}", "Annotation complete!".green()),
                    Err(e) => println!("{} {}", "Error during processing:".red(), e),
                }
                return;
            }

            // Split comma-separated query paths into Vec<String>
            let queries_value = queries
                .as_ref()
                .expect("--queries is required unless --kit is provided");
            let query_files: Vec<String> = queries_value
                .split(',')
                .map(|s| s.trim().to_string())
                .collect();

            let query_files_refs: Vec<&str> = query_files.iter().map(|s| s.as_str()).collect();

            // Parse barcode types
            let barcode_types_vec: Vec<BarcodeType> = barcode_types
                .split(',')
                .map(|s| match s.trim() {
                    "Ftag" => BarcodeType::Ftag,
                    "Rtag" => BarcodeType::Rtag,
                    _ => {
                        panic!("Unknown barcode type: {s}, use one of: Ftag, Rtag")
                    }
                })
                .collect();

            match annotate_with_files(
                input,
                query_files_refs,
                barcode_types_vec,
                output,
                *flank_max_errors,
                *alpha,
                *threads as u32,
                *verbose,
                *min_score,
                *min_score_diff,
            ) {
                // Convert fractions to raw scores
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
            only_side,
            failed_out,
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
                *only_side,
                failed_out.clone(),
                true, // Maybe make this optional but dont see a reason why you would not want this
            );
        }

        Commands::Inspect {
            input,
            top_n,
            read_pattern_out,
            bucket_size,
        } => {
            println!("{}", "Inspecting...".green());

            match inspect::inspect(input, *top_n, read_pattern_out.clone(), *bucket_size) {
                Ok(_) => println!("{}", "Inspection complete!".green()),
                Err(e) => println!("{} {}", "Inspection failed:".red(), e),
            }
        }

        Commands::Kit {
            kit,
            input,
            threads,
            output,
            maximize,
            verbose,
            min_score,
            min_score_diff,
            flank_max_errors,
            failed_out,
            use_extended,
            alpha,
        } => {
            demux_using_kit(
                kit.as_str(),
                input,
                *threads,
                output,
                *maximize,
                *verbose,
                *min_score,
                *min_score_diff,
                *flank_max_errors,
                failed_out.clone(),
                *use_extended,
                *alpha,
            );
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
