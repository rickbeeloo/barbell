use clap::{Parser, Subcommand};

#[derive(Parser)]
pub struct CompareArgs {
    #[arg(short = 'o', long, default_value = ".")]
    pub output_dir: String,

    #[arg(short = 'r', long, default_value = "data/rapid_bars.txt")]
    pub fastq_file: String,

    #[arg(short = 't', long, default_value = "1")]
    pub threads: usize,

    #[arg(short = 'e', long, default_value = "data/flexiplex_barcodes.flex")]
    pub extra_file: String,

    #[arg(short = 'd', long, default_value = "dorado")]
    pub dorado_exec_path: String,

    #[arg(short = 'b', long, default_value = "barbell")]
    pub barbell_exec_path: String,

    #[arg(short = 'f', long, default_value = "flexiplex")]
    pub flexiplex_exec_path: String,
}
