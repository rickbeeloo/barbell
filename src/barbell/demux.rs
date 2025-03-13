use crate::barbell::strategy::*;
use crate::barbell::reader::QueryGroup;
use crate::barbell::seq::Match;
use seq_io::fastq::{Reader as FastqReader, Record}; // both Fasta and Fastq import same name struct
use seq_io::parallel::read_parallel;
use std::io::Write;
use indicatif::ProgressStyle;
use std::time::Instant;
use std::sync::atomic::Ordering;
use colored::Colorize;

pub struct Demuxer<S: Strategy + Send + Sync> {
    strategy: S,
}

impl<S: Strategy + Send + Sync> Demuxer<S> {
    pub fn new(strategy: S) -> Self {
        Self { strategy }
    }

    pub fn demux_read(&self, read: &[u8], query_groups: &[QueryGroup]) -> Vec<Match> {
        // First we annotate the read 
        let annotations = self.strategy.annotate(query_groups, read);

        // Then we final assign the read
        self.strategy.final_assignment(&annotations)
    }

    pub fn auto_tune_parmas(&mut self, query_groups: &[QueryGroup]) {
        self.strategy.auto_tune_parmas(query_groups);
    }


    pub fn demux_fastq(&mut self, read_file: &str, query_groups: &[QueryGroup], auto_tune: bool, threads: u32, output_file: &str) {
        let start_time = Instant::now();

        println!("\n{}", "Configuration".bold().underline());
        println!("  • Input:  {}", read_file.bold());
        println!("  • Output: {}", output_file.bold());
        println!("  • Threads: {}", threads.to_string().bold());
        println!("  • Auto-tune: {}\n", if auto_tune { "yes".green() } else { "no".red() });

        if auto_tune {
            println!("{}", "Parameter Tuning".bold().underline());
            println!("  • Range: {} - {}", "0.05".dimmed(), "0.50".dimmed());
            println!("  • Test sequences: {}\n", "10,000".dimmed());
            
            self.auto_tune_parmas(query_groups);
        }

        // Create output file
        let output_file = std::fs::File::create(output_file).expect("Failed to create output file");
        let mut writer = std::io::BufWriter::new(output_file);
        writeln!(writer, "read\tlabel\tstart\tend\tedits\tdist.to.end\tread.len").expect("Failed to write header");

        // Setup progress bars
        let multi_progress = indicatif::MultiProgress::new();
        let total_bar = multi_progress.add(indicatif::ProgressBar::new_spinner());
        let found_bar = multi_progress.add(indicatif::ProgressBar::new_spinner());
        let missed_bar = multi_progress.add(indicatif::ProgressBar::new_spinner());

        // Style the progress bars
        total_bar.set_style(
            ProgressStyle::with_template("{spinner:.blue} {prefix:<12} {msg:>6} {elapsed_precise}")
                .unwrap()
                .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
        );
        found_bar.set_style(
            ProgressStyle::with_template("{spinner:.green} {prefix:<12} {msg:>6} {elapsed_precise}")
                .unwrap()
                .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
        );
        missed_bar.set_style(
            ProgressStyle::with_template("{spinner:.red} {prefix:<12} {msg:>6} {elapsed_precise}")
                .unwrap()
                .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
        );

        total_bar.set_prefix("Total:");
        found_bar.set_prefix("Found:");
        missed_bar.set_prefix("Missed:");

        // Counters
        let total = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let found = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let missed = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));

        let reader = FastqReader::from_path(read_file).expect("Failed to open read file");
        read_parallel(reader, threads, 1000, |record_set| {
            // this function does the heavy work
            let mut results = Vec::with_capacity(record_set.len());
            for record in record_set.into_iter() {
                let read_id = record.id().unwrap().to_string();
                let read = record.seq();
                let matches = self.demux_read(read, query_groups);
                results.push((read_id, matches, read.len()));
            }
            
            results
        }, |record_sets| {
            // This function runs in the main thread. It provides a streaming iterator over
            // record sets and the corresponding return values from the worker function
            // (not necessarily in the same order as in the file)
            while let Some(result) = record_sets.next() {
                let (_, found_results) = result.unwrap();
                for (read_id, matches, read_len) in found_results {
                    let total_count = total.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    total_bar.set_message(format!("{}", total_count));
                    
                    if !matches.is_empty() {
                        let found_count = found.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        found_bar.set_message(format!("{}", found_count));
                        
                        for m in matches {
                            writeln!(
                                writer,
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            read_id, m.some_id, m.start, m.end, m.edits, m.rel_dist_to_end, read_len
                            ).expect("Failed to write to output file");
                        }
                    } else {
                        let missed_count = missed.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        missed_bar.set_message(format!("{}", missed_count));
                        
                        writeln!(
                            writer,
                            "{}\tno_tag\t{}\t{}\t{}\t{}\t{}",
                            read_id, 0, 0, 0, 0, read_len
                        ).expect("Failed to write to output file");
                    }
                }
            }
        });

        // Finish progress bars
        total_bar.finish_with_message("Done!");
        found_bar.finish_with_message("Done!");
        missed_bar.finish_with_message("Done!");

        println!("\n{}", "Summary".bold().underline());
        println!("  • Time: {} seconds", start_time.elapsed().as_secs().to_string().bold());
        println!("  • Total reads: {}", total.load(Ordering::Relaxed).to_string().bold());
        println!("  • Tagged reads: {}", found.load(Ordering::Relaxed).to_string().green().bold());
        println!("  • Missed reads: {}\n", missed.load(Ordering::Relaxed).to_string().red().bold());
    }
}