use crate::annotate::annotator::annotate_with_kit;
use crate::config::{AnnotateConfig, FilterConfig, KitConfig, TrimConfig};
use crate::filter::filter::filter;
use crate::inspect::inspect::inspect;
use crate::kits::kits::*;
use crate::trim::trim::{LabelSide, trim_matches};
use anyhow::anyhow;
use colored::*;
use std::path::Path;

pub fn demux_using_kit(fastq_file: &str, config: &KitConfig) -> anyhow::Result<()> {
    let kit_name = config.kit_name.as_str();
    let output_folder = config.output_folder.as_str();
    // Create output folder if not exists yet
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder)?;
    }

    let kit_info = get_kit_info(kit_name);

    // Print some kit info
    println!("\n{}", "Kit info".purple().bold());
    println!("Kit name: {}", kit_info.name);
    println!(
        "Kit type: {}",
        if config.maximize { "Maximize" } else { "Safe" }
    );
    for tmpl in kit_info.templates {
        println!("Barcodes: {} - {}", tmpl.barcodes.from, tmpl.barcodes.to);
    }

    // If the default values are
    println!("\n{}", "Annotating reads...".purple().bold());
    let annotate_config = AnnotateConfig {
        max_flank_errors: config.max_flank_errors,
        alpha: config.alpha,
        n_threads: config.threads as u32,
        verbose: config.verbose,
        min_score: config.min_score,
        min_score_diff: config.min_score_diff,
        use_extended: config.use_extended,
    };
    annotate_with_kit(
        fastq_file,
        format!("{output_folder}/annotation.tsv").as_str(),
        kit_name,
        &annotate_config,
    )?;

    // // After annotating we show inspect
    println!("\n{}", "Top 10 most common patterns".purple().bold());
    let pattern_per_read_out = format!("{output_folder}/pattern_per_read.tsv");
    inspect(
        format!("{output_folder}/annotation.tsv").as_str(),
        10,
        Some(pattern_per_read_out),
        250,
    )
    .map_err(|e| anyhow!("{e}"))?;
    println!(
        "Want to see more patterns? Run: `barbell inspect {output_folder}/annotation.tsv -n 100`"
    );

    // Filter
    println!("\n{}", "Filtering reads...".purple().bold());

    let patterns = if config.maximize {
        (kit_info.maximize_patterns)()
    } else {
        (kit_info.safe_patterns)()
    };
    let filter_config = FilterConfig {
        verbose: config.verbose,
    };

    filter(
        format!("{output_folder}/annotation.tsv").as_str(),
        format!("{output_folder}/filtered.tsv").as_str(),
        None,
        patterns,
        &filter_config,
    )
    .map_err(|e| anyhow!("{e}"))?;

    // Trimming
    println!("\n{}", "Trimming reads...".purple().bold());
    let trim_config = TrimConfig {
        add_labels: true,
        add_orientation: false,
        add_flank: false,
        sort_labels: false,
        only_side: Some(LabelSide::Left),
        failed_trimmed_writer: config.failed_out.clone(),
        write_full_header: true,
        skip_trim: false,
        flip: false,
        verbose: config.verbose,
        output_format: config.output_format,
    };
    trim_matches(
        format!("{output_folder}/filtered.tsv").as_str(),
        fastq_file,
        output_folder,
        &trim_config,
    )?;

    println!("\n{}", "Done!".green().bold());
    Ok(())
}
