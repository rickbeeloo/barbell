use crate::trim::trim::LabelSide;

#[derive(Debug, Clone)]
pub struct AnnotateConfig {
    pub max_flank_errors: Option<usize>,
    pub alpha: f32,
    pub n_threads: u32,
    pub verbose: bool,
    pub min_score: f64,
    pub min_score_diff: f64,
    pub use_extended: bool,
}

#[derive(Debug, Clone)]
pub struct FilterConfig {
    pub verbose: bool,
}

/// Configuration for trimming the Nanopore ligation adapters (kit14 / `SQK-LSK114`
/// and friends) from the read ends. See [`crate::trim::adapter`].
///
/// This is only used by the standalone `barbell trim-adapters` command; the
/// barcode-based `barbell trim` / `barbell kit` pipelines do not trim adapters.
#[derive(Debug, Clone)]
pub struct AdapterTrimConfig {
    /// Minimum number of aligned bases before a hit is considered an adapter.
    pub min_match_len: usize,
    /// Minimum identity of the aligned part of the adapter.
    pub min_identity: f64,
    /// How far the outer edge of an adapter may be from the very start/end of the
    /// read and still be considered a terminal adapter.
    pub end_slack: usize,
    /// Number of threads to use.
    pub threads: usize,
    /// Write output FASTQ files as gzip.
    pub gzip: bool,
}

impl Default for AdapterTrimConfig {
    fn default() -> Self {
        Self {
            min_match_len: crate::trim::adapter::ADAPTER_MIN_MATCH_LEN,
            min_identity: crate::trim::adapter::ADAPTER_MIN_IDENTITY,
            end_slack: crate::trim::adapter::ADAPTER_END_SLACK,
            threads: crate::trim::adapter::DEFAULT_ADAPTER_THREADS,
            gzip: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct TrimConfig {
    pub add_labels: bool,
    pub add_orientation: bool,
    pub add_flank: bool,
    pub sort_labels: bool,
    pub only_side: Option<LabelSide>,
    pub failed_trimmed_writer: Option<String>,
    pub write_full_header: bool,
    pub skip_trim: bool,
    pub flip: bool,
    pub verbose: bool,
    pub gzip: bool,
}

#[derive(Debug, Clone)]
pub struct KitConfig {
    pub kit_name: String,
    pub threads: usize,
    pub output_folder: String,
    pub maximize: bool,
    pub verbose: bool,
    pub min_score: f64,
    pub min_score_diff: f64,
    pub max_flank_errors: Option<usize>,
    pub failed_out: Option<String>,
    pub use_extended: bool,
    pub alpha: f32,
    pub gzip: bool,
}
