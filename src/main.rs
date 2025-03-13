use clap::{Parser, Subcommand};
use colored::*;
mod barbell;
use barbell::reader::*;
use barbell::demux::*;
use crate::barbell::strategy::*;

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
    /// Plot results (not implemented yet)
    Plot,
}

fn main() {
    print_banner();
    
    let cli = Cli::parse();

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
    "#.bright_green());
    println!("{}",
        "        [===]------------------------------------------[===]        "
        .bright_yellow());
}
