use crate::kits::kits::*;
use std::io::Write;
use std::path::Path;

use crate::annotate::annotator::annotate_with_groups;
use crate::annotate::barcodes::BarcodeGroup;
use crate::annotate::{annotator::annotate, barcodes::BarcodeType};
use crate::filter::filter::filter;
use crate::inspect::inspect::inspect;
use crate::pattern_from_str;
use crate::trim::trim::{LabelSide, trim_matches};
use colored::*;

pub fn demux_using_kit(
    kit_name: &str,
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
    // Create output folder if not exists yet
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
    }

    let kit_info = get_kit_info(kit_name);
    let query_groups = BarcodeGroup::new_from_kit(kit_name);

    // If the default values are
    println!("\n{}", "Running annotation".purple().bold());
    annotate_with_groups(
        fastq_file,
        format!("{output_folder}/annotation.tsv").as_str(),
        query_groups,
        max_flank_errros,
        0.5,
        threads as u32,
        false,
        min_score,
        min_score_diff,
    )
    .expect("Annotation failed");

    // // After annotating we show inspect
    println!("\n{}", "Top 20 most common patterns".purple().bold());
    let pattern_per_read_out = format!("{output_folder}/pattern_per_read.tsv");
    inspect(
        format!("{output_folder}/annotation.tsv").as_str(),
        20,
        Some(pattern_per_read_out),
        250,
    )
    .expect("Inspect failed");

    // Filter
    println!("\n{}", "Running annotation filter".purple().bold());

    let patterns = if maximize {
        (kit_info.maximize_patterns)()
    } else {
        (kit_info.safe_patterns)()
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
