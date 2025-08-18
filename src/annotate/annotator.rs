use crate::annotate::barcodes::{BarcodeGroup, BarcodeType};
use crate::annotate::edit_model::get_edit_cut_off;
use crate::annotate::searcher::Demuxer;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use seq_io::fastq::{Reader, Record};
use seq_io::parallel::read_parallel;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread_local;
use std::time::Duration;

fn create_progress_bar() -> (ProgressBar, ProgressBar, ProgressBar) {
    // Create multiprogress bar
    let multi_progress = MultiProgress::new();
    let total_bar = multi_progress.add(ProgressBar::new_spinner());
    let found_bar = multi_progress.add(ProgressBar::new_spinner());
    let missed_bar = multi_progress.add(ProgressBar::new_spinner());

    // Style the progress bars with colors - more compact layout
    total_bar.set_style(
        ProgressStyle::with_template(
            "{spinner:.cyan} {prefix:.bold.white:<8} {msg:.bold.cyan:>6} {elapsed:.dim}",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    found_bar.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} {prefix:.bold.white:<8} {msg:.bold.green:>6} {elapsed:.dim}",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    missed_bar.set_style(
        ProgressStyle::with_template(
            "{spinner:.red} {prefix:.bold.white:<8} {msg:.bold.red:>6} {elapsed:.dim}",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );

    total_bar.enable_steady_tick(Duration::from_millis(100));
    found_bar.enable_steady_tick(Duration::from_millis(120)); // Slightly different timing for visual variety
    missed_bar.enable_steady_tick(Duration::from_millis(140));

    total_bar.set_prefix("Total:");
    found_bar.set_prefix("Found:");
    missed_bar.set_prefix("Missed:");

    (total_bar, found_bar, missed_bar)
}

pub fn annotate(
    read_file: &str,
    query_files: Vec<&str>,
    query_types: Vec<BarcodeType>,
    out_file: &str,
    max_flank_errors: Option<usize>,
    alpha: f32,
    n_threads: u32,
    verbose: bool,
    min_score: f64,
    min_score_diff: f64,
) -> anyhow::Result<()> {
    let reader = Reader::from_path(read_file).unwrap();
    let writer = Arc::new(Mutex::new(
        csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(out_file)
            .unwrap(),
    ));

    // Get query groups
    let mut query_groups = Vec::new();
    for (query_file, query_type) in query_files.iter().zip(query_types.iter()) {
        let mut query_group = BarcodeGroup::new_from_fasta(query_file, query_type.clone());
        if let Some(max_flank_errors) = max_flank_errors {
            query_group.set_flank_threshold(max_flank_errors);
        } else {
            // Determine based on formula
            let edit_cut_off = get_edit_cut_off(query_group.get_effective_len());
            println!("Auto edit flank cut off: {}", edit_cut_off);
            query_group.set_flank_threshold(edit_cut_off);
        }
        query_groups.push(query_group);
    }
    // Dispaly to user
    for (i, query_group) in query_groups.iter().enumerate() {
        println!("{}: {}", query_group.barcode_type.as_str(), i);
        query_group.display();
    }

    // Create progress bars
    let (total_bar, found_bar, missed_bar) = create_progress_bar();

    // Track counts, we need atomics here as multiple threads update the counters
    // we then use these to update the progressbars again
    let total = Arc::new(AtomicUsize::new(0));
    let found = Arc::new(AtomicUsize::new(0));
    let missed = Arc::new(AtomicUsize::new(0));

    read_parallel(
        reader,
        n_threads,
        1000,
        |record_set| {
            // Create thread local demuxer if not init for current thread yet
            thread_local! {
                static DEMUXER: std::cell::RefCell<Option<Demuxer>> = std::cell::RefCell::new(None);
            }
            DEMUXER.with(|cell| {
                if cell.borrow().is_none() {
                    let mut demux = Demuxer::new(alpha, verbose, min_score, min_score_diff);
                    for query_group in query_groups.iter() {
                        demux.add_query_group(query_group.clone());
                    }

                    *cell.borrow_mut() = Some(demux);
                }
            });

            // Go over the
            let mut record_set_annotations = Vec::new();
            // this function does the heavy work
            // total, fount we track here
            let mut found = 0;
            for record in record_set.into_iter() {
                // Use the demuxer through thread-local storage
                let matches: Vec<crate::annotate::searcher::BarbellMatch> = DEMUXER.with(|cell| {
                    if let Some(ref mut demuxer) = *cell.borrow_mut() {
                        demuxer.demux(record.id().unwrap(), record.seq())
                    } else {
                        panic!("Demuxer not initialized");
                    }
                });

                if !matches.is_empty() {
                    found += 1;
                    record_set_annotations.extend(matches);
                }
            }
            Some((record_set_annotations, found))
        },
        |record_sets| {
            // This function runs in the main thread. It provides a streaming iterator over
            // record sets and the corresponding return values from the worker function
            // (not necessarily in the same order as in the file)
            while let Some(result) = record_sets.next() {
                let (record_set, found_result) = result.unwrap();
                let n_records = record_set.len();
                total.fetch_add(n_records, Ordering::Relaxed);
                if let Some((record_set_results, found_count)) = found_result {
                    // Write to output file, we lock for entire record set
                    let mut writer = writer.lock().unwrap();
                    for m in &record_set_results {
                        writer.serialize(m).unwrap();
                    }
                    drop(writer);

                    // Update found count
                    found.fetch_add(found_count, Ordering::Relaxed);
                    missed.fetch_add(n_records - found_count, Ordering::Relaxed);
                }

                // Update progress bars
                let total_count = total.load(Ordering::Relaxed);
                let found_count = found.load(Ordering::Relaxed);
                let missed_count = missed.load(Ordering::Relaxed);

                total_bar.set_message(total_count.to_string());
                found_bar.set_message(found_count.to_string());
                missed_bar.set_message(missed_count.to_string());
            }
        },
    );

    // Print final summary
    let final_total = total.load(Ordering::Relaxed);
    let final_found = found.load(Ordering::Relaxed);
    let final_missed = missed.load(Ordering::Relaxed);

    // Finish progress bars
    total_bar.finish_with_message(format!("Done: {final_total} records"));
    found_bar.finish_with_message(format!("Found: {final_found} records"));
    missed_bar.finish_with_message(format!("Missed: {final_missed} records"));

    Ok(())
}
