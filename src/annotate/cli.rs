#[derive(clap::Parser, Clone)]
pub struct AnnotateArgs {
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
}
