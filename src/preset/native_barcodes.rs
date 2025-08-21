use std::io::Write;
use std::path::Path;

use crate::annotate::{annotator::annotate, barcodes::BarcodeType};
use crate::filter::filter::filter;
use crate::inspect::inspect::inspect;
use crate::pattern_from_str;
use crate::trim::trim::{LabelSide, trim_matches};
use colored::*;
use tempfile::NamedTempFile;

pub fn demux_native_barcodes(
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
    // Create temp file for string content
    let mut tmp_query_file = NamedTempFile::new().expect("Failed to create temporary file");
    tmp_query_file
        .write_all(crate::NATIVE_BARS_CONTENT.as_bytes())
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
        max_flank_errros,
        0.5,
        threads as u32,
        false,
        min_score,
        min_score_diff,
    )
    .expect("Annotation failed");

    // // After annotating we show inspect
    println!("\n{}", "Top 15 most common patterns".purple().bold());
    let pattern_per_read_out = format!("{output_folder}/pattern_per_read.tsv");
    inspect(
        format!("{output_folder}/annotation.tsv").as_str(),
        15,
        Some(pattern_per_read_out),
    )
    .expect("Inspect failed");

    // Filter
    println!("\n{}", "Running annotation filter".purple().bold());

    // The "safe" patterns
    // ideal - single barcode on the left side

    // Barcode on just the left side (fw)
    let pattern1 = pattern_from_str!("Ftag[fw, *, @left(0..250), >>]");

    // Barcode on just the right side (rc)
    let pattern2 = pattern_from_str!("Ftag[<<, rc, *, @right(0..250)]");

    // Barcode on both sides, with same label (?1 wildcard)
    let pattern3 =
        pattern_from_str!("Ftag[fw, ?1, @left(0..250), >>]__Ftag[<<, rc, ?1, @right(0..250)]");

    // The less safe patterns

    // Two barcodes on the left, where inner left and right side have the same label (?1 wildcard)
    let pattern4 = pattern_from_str!(
        "Ftag[fw, *, @left(0..250)]__Ftag[fw, ?1, @prev_left(0..250), >>]__Ftag[<<, rc, ?1, @right(0..250)]"
    );

    // Barcode on the left, missing right
    let pattern5 =
        pattern_from_str!("Ftag[fw, *, @left(0..250), >>]__Fflank[<<, rc, *, @right(0..250)]");

    // Barcode on the right, missing left
    let pattern7 =
        pattern_from_str!("Fflank[fw, *, @left(0..250), >>]__Ftag[<<, rc, *, @right(0..250)]");

    // Two barcodes on the left side
    let pattern6 =
        pattern_from_str!("Ftag[fw, *, @left(0..250)]__Ftag[fw, *, @prev_left(0..250), >>]");

    // Two barcodes on the right side, cut inner
    let pattern8 = pattern_from_str!(
        "Ftag[<<, fw, ?1, @left(0..250), >>]__Ftag[<<, fw, ?1, @right(0..250)]__Ftag[rc, *, @right(0..250)]"
    );

    // Three barcodes left side, one on the right should be the same as inner left
    let pattern9 = pattern_from_str!(
        "Ftag[fw, *, @left(0..250)]__Ftag[rc, *, @prev_left(0..250)]__Ftag[fw, ?1, @prev_left(0..250), >>]__Ftag[<<, rc, ?1, @right(0..250)]"
    );

    let patterns = if maximize {
        vec![pattern1, pattern2, pattern3]
    } else {
        vec![
            pattern1, pattern2, pattern3, pattern4, pattern5, pattern6, pattern7, pattern8,
            pattern9,
        ]
    };

    filter(
        format!("{output_folder}/annotation.tsv").as_str(),
        format!("{output_folder}/filtered.tsv").as_str(),
        None,
        patterns,
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
        false,
        false,
        Some(LabelSide::Left),
        failed_out,
    );
}
