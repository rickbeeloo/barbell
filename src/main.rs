use clap::{Parser, Subcommand};
use colored::*;
use nano_bar::annotate::mutations::ErrorRates;
use nano_bar::annotate::search::BarMan;
use nano_bar::annotate::*;
use nano_bar::filter::filter::*;
use nano_bar::parallel::ParallelAnnotator;
use nano_bar::read::reader::*;
use nano_bar::inspect::inspect::*;
use nano_bar::trim::trim::*;
use std::env;

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
        #[arg(short = 't', long, default_value_t = 5)]
        threads: usize,

        /// Output file path
        #[arg(short = 'o', long, default_value = "output.csv")]
        output: String,

        /// Query files (comma-separated paths)
        #[arg(short = 'q', long)]
        queries: String,

        /// Tuning
        #[arg(long, default_value_t = false)]
        tune: bool,

        /// Insertion error rate
        #[arg(long, default_value_t = 0.03)]
        ins_error: f64,

        /// Deletion error rate
        #[arg(long, default_value_t = 0.1)]
        del_error: f64,

        /// Substitution error rate
        #[arg(long, default_value_t = 0.03)]
        sub_error: f64,

        /// Log probability cutoff (use '=' for negative values, e.g. --log-prob-cutoff=-21.0)
        #[arg(long, allow_negative_numbers = true, default_value_t = -10.0)]
        log_prob_cutoff: f64,

        /// Target false positive rate
        #[arg(long, default_value_t = 0.00001)]
        fp_target: f64,
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

       /// Output file path for sorted reads
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

    // View most common patterns in annotation
    Inspect {
        /// Input filtered annotation file
        #[arg(short = 'i', long)]
        input: String,

        /// Top N
        #[arg(short = 'n', long, default_value_t = 10)]
        top_n: usize,
    }
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
            ins_error,
            del_error,
            sub_error,
            log_prob_cutoff,
            tune,
            fp_target
        } => {
            println!("{}", "Starting annotation...".green());
            
            // Split comma-separated query paths into Vec<String>
            let query_files: Vec<String> = queries.split(',')
                .map(|s| s.trim().to_string())
                .collect();
            
            let query_files_refs: Vec<&str> = query_files.iter().map(|s| s.as_str()).collect();
            let groups = read_queries(query_files_refs, None);
            
            let errors = ErrorRates::new(*ins_error, *del_error, *sub_error);
            let mut bar_searcher = BarMan::new(
                groups, 
                0.4, 
                errors,
                0.4, // Will be tuned when --tune
                0.5,
                *log_prob_cutoff, 
                0.9,
                *fp_target
            );
            if *tune {
                bar_searcher.auto_tune_parmas();
            }

            //REMOVE
            //bar_searcher.posterior_threshold = -21.0;
            
            let mut parallel_annotator = ParallelAnnotator::new(bar_searcher);
            
            match parallel_annotator.process_fastq(input, output, *threads) {
                Ok(_) => println!("{}", "Annotation complete!".green()),
                Err(e) => println!("{} {}", "Error during processing:".red(), e),
            }
        },

        Commands::Filter { input, output, file } => {
            println!("{}", "Starting filtering...".green());
            
            match filter_from_text_file(input, file, output) {
                Ok(_) => println!("{}", "Filtering successful!".green()),
                Err(e) => println!("{} {}", "Filtering failed:".red(), e),
            }
        },

        Commands::Trim { input, reads, output, no_label, no_orientation, no_flanks, sort_labels } => {
            trim_matches(input, reads, output, !no_label, !no_orientation, !no_flanks, *sort_labels);
        },

        Commands::Inspect { input, top_n } => {
            println!("{}", "Inspecting...".green());
            
            match inspect(input, *top_n) {
                Ok(_) => println!("{}", "Inspection complete!".green()),
                Err(e) => println!("{} {}", "Inspection failed:".red(), e),
            }
        }
    }
}

fn print_banner() {
    println!("{}", r#"
    ███╗   ██╗ █████╗ ███╗   ██╗ ██████╗     ██████╗  █████╗ ██████╗ 
    ████╗  ██║██╔══██╗████╗  ██║██╔═══██╗    ██╔══██╗██╔══██╗██╔══██╗
    ██╔██╗ ██║███████║██╔██╗ ██║██║   ██║    ██████╔╝███████║██████╔╝
    ██║╚██╗██║██╔══██║██║╚██╗██║██║   ██║    ██╔══██╗██╔══██║██╔══██╗
    ██║ ╚████║██║  ██║██║ ╚████║╚██████╔╝    ██████╔╝██║  ██║██║  ██║
    ╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝  ╚═══╝ ╚═════╝     ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝
    "#.blue());
    println!("{}",
        "        [===]------------------------------------------[===]        "
        .bright_yellow());
}
