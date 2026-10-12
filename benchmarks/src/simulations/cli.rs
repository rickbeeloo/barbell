use clap::{Parser, Subcommand};

#[derive(Parser)]
pub struct SimulationArgs {
    #[arg(short = 'o', long, default_value = ".")]
    pub output_dir: String,

    #[arg(short = 'n', long, default_value = "1000")]
    pub n: usize,

    #[arg(short = 'b', long, default_value = "data/rapid_bars.txt")]
    pub barcode_file: String,

    #[arg(short = 'r', long, default_value = "0.0")]
    pub rc_frac: f64,
}
