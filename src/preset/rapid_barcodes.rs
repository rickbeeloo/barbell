use std::io::Write;
use std::path::Path;

use crate::annotate::{annotator::annotate, barcodes::BarcodeType};
use crate::filter::filter::filter;
use crate::inspect::inspect::inspect;
use crate::pattern_from_str;
use crate::trim::trim::trim_matches;
use colored::*;
use tempfile::NamedTempFile;

pub fn demux_rapid_barcodes(
    fastq_file: &str,
    threads: usize,
    output_folder: &str,
    max_error_perc: Option<f32>,
) {
    // Create temp file for string content
    let mut tmp_query_file = NamedTempFile::new().expect("Failed to create temporary file");
    tmp_query_file
        .write_all(crate::RAPID_BARS_CONTENT.as_bytes())
        .expect("Failed to write to temporary file");
    tmp_query_file
        .flush()
        .expect("Failed to flush temporary file");

    // Create output folder if not exists yet
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
    }

    // If the default values are

    println!("\n{}", "Running annotation".purple().bold());
    annotate(
        fastq_file,
        vec![tmp_query_file.path().to_str().unwrap()],
        vec![BarcodeType::Ftag],
        format!("{output_folder}/annotation.tsv").as_str(),
        None,
        Some(20),
        Some(7),
        0.5, // Overhang alpha
        threads as u32,
    )
    .expect("Annotation failed");

    // // After annotating we show inspect
    println!("\n{}", "Top 10 most common patterns".purple().bold());
    let pattern_per_read_out = format!("{output_folder}/pattern_per_read.tsv");
    inspect(
        format!("{output_folder}/annotation.tsv").as_str(),
        10,
        Some(pattern_per_read_out),
    )
    .expect("Inspect failed");

    // Filter
    println!("\n{}", "Running annotation filter".purple().bold());
    let pattern1 = pattern_from_str!("Ftag[fw, *, @left(0..250), >>]");
    let pattern2 =
        pattern_from_str!("Ftag[fw, *, @left(0..250)]__Ftag[fw, *, @prev_left(0..250), >>]");
    filter(
        format!("{output_folder}/annotation.tsv").as_str(),
        format!("{output_folder}/filtered.tsv").as_str(),
        None,
        vec![pattern1, pattern2],
    )
    .expect("Filter failed");

    // Trimming
    println!("\n{}", "Running trimming".purple().bold());
    trim_matches(
        format!("{output_folder}/filtered.tsv").as_str(),
        fastq_file,
        output_folder,
        true,
        true,
        true,
        true,
    );
}
