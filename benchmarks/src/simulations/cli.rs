use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Output directory for results
    #[arg(short, long, default_value = "output")]
    pub output_dir: String,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Run simulations
    #[command(name = "sim")]
    Simulate(SimulationArgs),
}

#[derive(Parser)]
pub struct SimulationArgs {
    #[arg(short = 'o', long, default_value = ".")]
    pub output_dir: String,

    #[arg(short = 'n', long, default_value = "1000")]
    pub n: usize,

    #[arg(short = 'b', long, default_value = "data/rapid_bars.txt")]
    pub barcode_file: String,
}
