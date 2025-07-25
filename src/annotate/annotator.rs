use crate::annotate::barcodes::{BarcodeGroup, BarcodeType};
use crate::annotate::searcher::Demuxer;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use once_cell::sync::Lazy;
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::{Reader, Record};
use seq_io::parallel::read_parallel;
use std::cell::RefCell;
use std::io::BufWriter;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread_local;
use thread_id;

// Macro to create thread-local demuxer with lazy initialization
macro_rules! thread_local_demuxer {
    ($name:ident, $alpha:expr, $query_files:expr, $query_types:expr, $max_error_perc:expr) => {
        thread_local! {
            static $name: std::cell::RefCell<Option<Demuxer>> = std::cell::RefCell::new(None);
        }

        $name.with(|cell| {
            if cell.borrow().is_none() {
                println!("Thread {:?} initializing demuxer", thread_id::get());
                *cell.borrow_mut() = Some(Demuxer::new($alpha));
                for (query_file, query_type) in $query_files.iter().zip($query_types.iter()) {
                    let mut query_group =
                        BarcodeGroup::new_from_fasta(query_file, query_type.clone());
                    query_group.tune_set_perc_threshold($max_error_perc);
                    cell.borrow_mut()
                        .as_mut()
                        .unwrap()
                        .add_query_group(query_group);
                }
            }
        });
    };
}

// Move thread_local outside the function
thread_local! {
    static DEMUXER: std::cell::RefCell<Option<Demuxer>> = std::cell::RefCell::new(None);
}

pub fn annotate(
    read_file: &str,
    query_files: Vec<&str>,
    query_types: Vec<BarcodeType>,
    out_file: &str,
    max_error_perc: f32,
    alpha: f32,
    n_threads: u32,
) -> anyhow::Result<()> {
    let reader = Reader::from_path(read_file).unwrap();
    let writer = Arc::new(Mutex::new(
        csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(out_file)
            .unwrap(),
    ));

    read_parallel(
        reader,
        n_threads,
        10,
        |record_set| {
            // Create thread local demuxer if not init for current thread yet
            thread_local_demuxer!(DEMUXER, alpha, query_files, query_types, max_error_perc);

            // Go over the
            let mut record_set_annotations = Vec::new();
            // this function does the heavy work
            for (i, record) in record_set.into_iter().enumerate() {
                // Use the demuxer through thread-local storage
                let matches = DEMUXER.with(|cell| {
                    if let Some(ref mut demuxer) = *cell.borrow_mut() {
                        demuxer.demux(record.id().unwrap(), record.seq())
                    } else {
                        panic!("Demuxer not initialized");
                    }
                });

                if !matches.is_empty() {
                    record_set_annotations.extend(matches);
                }
            }
            Some(record_set_annotations)
        },
        |record_sets| {
            // This function runs in the main thread. It provides a streaming iterator over
            // record sets and the corresponding return values from the worker function
            // (not necessarily in the same order as in the file)
            while let Some(result) = record_sets.next() {
                let (record_set, found) = result.unwrap();
                if let Some(record_set_results) = found {
                    // Write to output file, we lock for entire record set
                    let mut writer = writer.lock().unwrap();
                    for m in &record_set_results {
                        writer.serialize(m).unwrap();
                    }
                    drop(writer);
                }
            }
        },
    );

    Ok(())
}
