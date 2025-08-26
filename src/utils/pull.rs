use crate::annotate::searcher::BarbellMatch;
use crate::trim::trim::LabelConfig;
use csv;
use indicatif::MultiProgress;
use indicatif::{ProgressBar, ProgressStyle};
use seq_io::fastq::{Reader, Record};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Duration;

fn create_progress_bar() -> (ProgressBar, ProgressBar, ProgressBar, ProgressBar) {
    let multi_progress = MultiProgress::new();
    let total_bar = multi_progress.add(ProgressBar::new_spinner());
    let grouped_bar = multi_progress.add(ProgressBar::new_spinner());
    let written_bar = multi_progress.add(ProgressBar::new_spinner());
    let failed_bar = multi_progress.add(ProgressBar::new_spinner());

    total_bar.set_style(
        ProgressStyle::with_template(
            "{spinner:.cyan} {prefix:.bold.white:<8} {msg:.bold.cyan:>6} {elapsed:.dim}",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    grouped_bar.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} {prefix:.bold.white:<8} {msg:.bold.green:>6} {elapsed:.dim}",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    written_bar.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} {prefix:.bold.white:<8} {msg:.bold.red:>6} {elapsed:.dim}",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    failed_bar.set_style(
        ProgressStyle::with_template(
            "{spinner:.red} {prefix:.bold.white:<8} {msg:.bold.red:>6} {elapsed:.dim}",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );

    total_bar.enable_steady_tick(Duration::from_millis(100));
    grouped_bar.enable_steady_tick(Duration::from_millis(120));
    written_bar.enable_steady_tick(Duration::from_millis(140));
    failed_bar.enable_steady_tick(Duration::from_millis(160));

    total_bar.set_prefix("Total:");
    grouped_bar.set_prefix("Grouped:");
    written_bar.set_prefix("Written:");
    failed_bar.set_prefix("Failed:");

    (total_bar, grouped_bar, written_bar, failed_bar)
}

pub fn pull_reads_without_trimming(
    filtered_match_file: &str,
    read_fastq_file: &str,
    output_folder: &str,
    add_labels: bool,
    add_orientation: bool,
    add_flank: bool,
    sort_labels: bool,
    only_side: Option<crate::trim::trim::LabelSide>,
    failed_out: Option<String>,
) {
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
    }

    let label_config = LabelConfig::new(
        add_labels,
        add_orientation,
        add_flank,
        sort_labels,
        only_side,
    );

    let mut annotations_by_read: HashMap<String, Vec<BarbellMatch>> = HashMap::new();

    let (total_bar, grouped_bar, written_bar, failed_bar) = create_progress_bar();

    let mut matches_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(filtered_match_file)
        .expect("Failed to open matches file");

    for result in matches_reader.deserialize() {
        let anno: BarbellMatch = result.expect("Failed to parse annotation line");
        annotations_by_read
            .entry(anno.read_id.clone())
            .or_default()
            .push(anno);
    }

    let mut writers: HashMap<String, BufWriter<File>> = HashMap::new();

    let mut failed_writer = failed_out.map(|path| BufWriter::new(File::create(path).unwrap()));

    let mut reader = Reader::from_path(read_fastq_file).expect("Failed to open FASTQ file");

    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        // Unlike trimming we do not clean the id as we need to move info
        let (read_id, desc) = record.id_desc().unwrap();
        let read_id = read_id.to_string();
        let desc = desc.unwrap_or_default();
        let full_header = format!("{read_id} {desc}");
        total_bar.inc(1);

        if let Some(annotations) = annotations_by_read.get(&read_id) {
            grouped_bar.inc(1);

            let label = label_config.create_label(annotations);
            let group = label;

            let writer = writers.entry(group.clone()).or_insert_with(|| {
                let output_file = format!("{output_folder}/{group}.raw.fastq");
                BufWriter::new(
                    File::create(&output_file)
                        .unwrap_or_else(|_| panic!("Failed to create output file: {output_file}")),
                )
            });

            writeln!(writer, "@{}", full_header).expect("Failed to write header");
            writeln!(writer, "{}", String::from_utf8_lossy(record.seq()))
                .expect("Failed to write sequence");
            writeln!(writer, "+").expect("Failed to write separator");
            writeln!(writer, "{}", String::from_utf8_lossy(record.qual()))
                .expect("Failed to write quality");
            written_bar.inc(1);
        } else {
            failed_bar.inc(1);
            if let Some(ref mut failed_writer) = failed_writer {
                writeln!(failed_writer, "{}", full_header).expect("Failed to write failed id");
            }
        }

        total_bar.set_message(total_bar.position().to_string());
        grouped_bar.set_message(grouped_bar.position().to_string());
        written_bar.set_message(written_bar.position().to_string());
        failed_bar.set_message(failed_bar.position().to_string());
    }

    for (_, writer) in writers.iter_mut() {
        writer.flush().expect("Failed to flush output");
    }

    let total_count = total_bar.position();
    let grouped_count = grouped_bar.position();
    let written_count = written_bar.position();
    let failed_count = failed_bar.position();

    total_bar.set_message(total_count.to_string());
    grouped_bar.set_message(grouped_count.to_string());
    written_bar.set_message(written_count.to_string());
    failed_bar.set_message(failed_count.to_string());

    total_bar.finish_with_message(format!("{total_count} reads"));
    grouped_bar.finish_with_message(format!("{grouped_count} reads"));
    written_bar.finish_with_message(format!("{written_count} reads"));
    failed_bar.finish_with_message(format!("{failed_count} reads"));
}
