use clap::{Parser, Subcommand};
use colored::*;
mod barbell;
use barbell::reader::*;
use barbell::annotater::*;
use barbell::annotate_strategy::*;
use barbell::filter;
use barbell::trim;

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
        threads: u32,

        /// Enable autotuning
        #[arg(long)]
        tune: bool,

        /// Output file path
        #[arg(short = 'o', long)]
        output: String,

        /// Query files (comma-separated paths)
        #[arg(short = 'q', long)]
        queries: String,
    },
    /// Filter annotation files based on pattern
    Filter {
        /// Input annotation file
        #[arg(short = 'i', long, required = true)]
        input: String,

        /// Output filtered file path
        #[arg(short = 'o', long, required = true)]
        output: String,

        /// Pattern string to filter by
        #[arg(short = 'p', long, conflicts_with = "file", required_unless_present = "file")]
        pattern: Option<String>,

        /// File containing patterns to filter by
        #[arg(short = 'f', long, conflicts_with = "pattern", required_unless_present = "pattern")]
        file: Option<String>,
    },
    /// Trim and sort reads based on filtered annotations
    Trimm {
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
    },
    /// Plot results (not implemented yet)
    Plot,
}

fn main() {
    // Add debug prints
    eprintln!("Starting program");
    eprintln!("Current directory: {:?}", std::env::current_dir().unwrap());
    
    print_banner();
    eprintln!("Banner printed");
    
    let cli = Cli::parse();
    eprintln!("CLI parsed");

    match &cli.command {
        Commands::Annotate { input, threads, tune, output, queries } => {
            println!("{}", "Starting annotation...".green());
            
            // Split comma-separated query paths into Vec<String>
            let query_paths: Vec<&str> = queries.split(',')
                .map(|s| s.trim())
                .collect();
            
            let groups = read_queries(query_paths, None);
            
            let mut demuxer = Demuxer::new(SimpleStrategy::default());
            demuxer.demux_fastq(input, &groups, *tune, *threads, output);
            
            println!("{}", "Annotation complete!".green());
        }

        Commands::Filter { input, output, pattern, file } => {
            println!("{}", "Starting filtering...".green());
            
            let result = if let Some(pattern) = pattern {
                barbell::filter::filter_from_pattern_str(input, &pattern, output)
            } else if let Some(file) = file {
                barbell::filter::filter_from_text_file(input, &file, output)
            } else {
                unreachable!("Clap ensures either pattern or file is provided")
            };

            match result {
                Ok(_) => println!("{}", "Filtering complete!".green()),
                Err(e) => eprintln!("{} {}", "Filtering failed:".red(), e),
            }
        }

        Commands::Trimm { input, reads, output, no_label, no_orientation, no_flanks } => {
            println!("{}", "Starting read trimming and sorting...".green());
            
            barbell::trim::trim_matches(
                input, 
                reads, 
                output,
                !no_label,           // Invert the flags since the function expects positive logic
                !no_orientation,
                !no_flanks
            );
            
            println!("{}", "Trimming complete!".green());
        }

        Commands::Plot => {
            println!("{}", "Plot functionality not implemented yet".yellow());
        }
    }
}

fn print_banner() {
    println!("{}", r#"
    ██████╗  █████╗ ██████╗ ██████╗ ███████╗██╗     ██╗     
    ██╔══██╗██╔══██╗██╔══██╗██╔══██╗██╔════╝██║     ██║     
    ██████╔╝███████║██████╔╝██████╔╝█████╗  ██║     ██║     
    ██╔══██╗██╔══██║██╔══██╗██╔══██╗██╔══╝  ██║     ██║     
    ██████╔╝██║  ██║██║  ██║██████╔╝███████╗███████╗███████╗
    ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═════╝ ╚══════╝╚══════╝╚══════╝
    "#.blue());
    println!("{}",
        "        [===]------------------------------------------[===]        "
        .bright_yellow());
}
