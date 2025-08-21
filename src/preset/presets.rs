use crate::preset::native_barcodes::demux_native_barcodes;
use crate::preset::rapid_barcodes::demux_rapid_barcodes;

#[derive(Debug, Clone, PartialEq, clap::ValueEnum)]
pub enum PresetName {
    Rapid,
    Native,
}

pub fn use_preset(
    preset: PresetName,
    fastq_file: &str,
    threads: usize,
    output_folder: &str,
    maximize: bool,
    verbose: bool,
    min_score: f64,
    min_score_diff: f64,
    max_flank_errros: Option<usize>,
    failed_out: Option<String>,
) {
    match preset {
        PresetName::Rapid => demux_rapid_barcodes(
            fastq_file,
            threads,
            output_folder,
            maximize,
            verbose,
            min_score,
            min_score_diff,
            max_flank_errros,
            failed_out,
        ),
        PresetName::Native => demux_native_barcodes(
            fastq_file,
            threads,
            output_folder,
            maximize,
            verbose,
            min_score,
            min_score_diff,
            max_flank_errros,
            failed_out,
        ),
    }
}
