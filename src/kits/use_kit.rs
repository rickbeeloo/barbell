use crate::annotate::annotator::annotate_with_kit;
use crate::filter::filter::filter;
use crate::inspect::inspect::inspect;
use crate::kits::kits::*;
use crate::trim::trim::{LabelSide, trim_matches};
use colored::*;
use std::path::Path;

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
    use_extended: bool,
    alpha: f32,
) {
    // Create output folder if not exists yet
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
    }

    let kit_info = get_kit_info(kit_name);

    // Print some kit info
    println!("\n{}", "Kit info".purple().bold());
    println!("Kit name: {}", kit_info.name);
    println!("Kit type: {}", if maximize { "Maximize" } else { "Safe" });
    for tmpl in kit_info.templates {
        println!("Barcodes: {} - {}", tmpl.barcodes.from, tmpl.barcodes.to);
    }

    // If the default values are
    println!("\n{}", "Annotating reads...".purple().bold());
    annotate_with_kit(
        fastq_file,
        format!("{output_folder}/annotation.tsv").as_str(),
        kit_name,
        max_flank_errros,
        alpha,
        threads as u32,
        verbose,
        min_score,
        min_score_diff,
        use_extended,
    )
    .expect("Annotation failed");

    // // After annotating we show inspect
    println!("\n{}", "Top 10 most common patterns".purple().bold());
    let pattern_per_read_out = format!("{output_folder}/pattern_per_read.tsv");
    inspect(
        format!("{output_folder}/annotation.tsv").as_str(),
        10,
        Some(pattern_per_read_out),
        250,
    )
    .expect("Inspect failed");
    println!(
        "Want to see more patterns? Run: `barbell inspect {output_folder}/annotation.tsv -n 100`"
    );

    // Filter
    println!("\n{}", "Filtering reads...".purple().bold());

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
    println!("\n{}", "Trimming reads...".purple().bold());
    trim_matches(
        format!("{output_folder}/filtered.tsv").as_str(),
        fastq_file,
        output_folder,
        true,
        false, // probalby nobody using kits uses this (only for custom)
        false,
        false,
        Some(LabelSide::Left),
        failed_out,
        true,
    );

    println!("\n{}", "Done!".green().bold());
}
