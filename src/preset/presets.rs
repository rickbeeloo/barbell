use crate::preset::rapid_barcodes::demux_rapid_barcodes;

#[derive(Debug, Clone, PartialEq, clap::ValueEnum)]
pub enum PresetName {
    Rapid,
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
        ),
        _ => panic!("Unsupported preset, use one of: {:?}", PresetName::Rapid),
    }
}
