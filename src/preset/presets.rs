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
    flank_max_errors: usize,
    barcode_max_errors: usize,
) {
    match preset {
        PresetName::Rapid => demux_rapid_barcodes(
            fastq_file,
            threads,
            output_folder,
            flank_max_errors,
            barcode_max_errors,
        ),
        _ => panic!("Unsupported preset, use one of: {:?}", PresetName::Rapid),
    }
}
