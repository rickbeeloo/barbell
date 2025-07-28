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
    max_edit_perc: f32,
) {
    match preset {
        PresetName::Rapid => {
            demux_rapid_barcodes(fastq_file, threads, output_folder, max_edit_perc)
        }
        _ => panic!("Unsupported preset, use one of: {:?}", PresetName::Rapid),
    }
}
