use crate::barbell::annotate_strategy::*;
use crate::barbell::reader::QueryGroup;
use seq_io::fastq::{Reader as FastqReader, Record}; // both Fasta and Fastq import same name struct
use seq_io::parallel::read_parallel;
use std::io::Write;
use indicatif::ProgressStyle;
use std::time::Instant;
use std::sync::atomic::Ordering;
use colored::Colorize;
use seq_io_parallel::{MinimalRefRecord, ParallelProcessor, ParallelReader};
use std::sync::{Arc, Mutex};
use std::io::BufWriter;
use anyhow::Result;
use crate::barbell::merge_sort::merge_sort_files;
use crate::barbell::pattern_assign::Match;

#[derive(Clone)]
pub struct Demuxer<S: Strategy + Send + Sync + Clone> {
    strategy: S,
    writer: Arc<Mutex<BufWriter<std::fs::File>>>,
    query_groups: Arc<Vec<QueryGroup>>,
    total: Arc<std::sync::atomic::AtomicUsize>,
    found: Arc<std::sync::atomic::AtomicUsize>,
    missed: Arc<std::sync::atomic::AtomicUsize>,
    total_bar: Arc<indicatif::ProgressBar>,
    found_bar: Arc<indicatif::ProgressBar>,
    missed_bar: Arc<indicatif::ProgressBar>,
}

impl<S: Strategy + Send + Sync + Clone> Demuxer<S> {
    pub fn new(strategy: S) -> Self {
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

        Self { 
            strategy,
            writer: Arc::new(Mutex::new(BufWriter::new(std::fs::File::create("/dev/null").unwrap()))),
            query_groups: Arc::new(Vec::new()),
            total: Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            found: Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            missed: Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            total_bar: Arc::new(total_bar),
            found_bar: Arc::new(found_bar),
            missed_bar: Arc::new(missed_bar),
        }
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

    pub fn demux_fastq(&mut self, read_file: &str, query_groups: &[QueryGroup], auto_tune: bool, threads: u32, output_file_path: &str) {
        let start_time = Instant::now();

        println!("\n{}", "Configuration".bold().underline());
        println!("  • Input:  {}", read_file.bold());
        println!("  • Output: {}", output_file_path.bold());
        println!("  • Threads: {}", threads.to_string().bold());
        println!("  • Auto-tune: {}\n", if auto_tune { "yes".green() } else { "no".red() });

        if auto_tune {           
            self.auto_tune_parmas(query_groups);
        }

        // Create output file and initialize writer
        let output_file = std::fs::File::create(output_file_path).expect("Failed to create output file");
        let writer = BufWriter::new(output_file);
        self.writer = Arc::new(Mutex::new(writer));
        self.query_groups = Arc::new(query_groups.to_vec());

        // Write header
        let mut writer = self.writer.lock().unwrap();
        writeln!(writer, "read\tlabel\tstart\tend\tedits\tdist.to.end\tread.len\trecord_set_idx\trecord_idx")
            .expect("Failed to write header");
        drop(writer);

        // Process reads
        let reader = FastqReader::from_path(read_file).expect("Failed to open read file");
        reader.process_parallel(self.clone(), threads as usize).expect("Failed to process reads");

        // Finish progress bars
        self.total_bar.finish_with_message("Done!");
        self.found_bar.finish_with_message("Done!");
        self.missed_bar.finish_with_message("Done!");


        println!("\n{}", "Summary".bold().underline());
        println!("  • Time: {} seconds", start_time.elapsed().as_secs().to_string().bold());
        println!("  • Total reads: {}", self.total.load(Ordering::Relaxed).to_string().bold());
        println!("  • Tagged reads: {}", self.found.load(Ordering::Relaxed).to_string().green().bold());
        println!("  • Missed reads: {}\n", self.missed.load(Ordering::Relaxed).to_string().red().bold());

        //  Run merge sort on output file to get original order
        println!("  • Sorting output file (merge sort)");
        merge_sort_files(output_file_path).expect("Failed to merge sort output file");
    }
}

impl<S: Strategy + Send + Sync + Clone> ParallelProcessor for Demuxer<S> {
    fn process_record<'a, Rf: MinimalRefRecord<'a>>(&mut self, record: Rf, record_set_idx: usize, record_idx: usize) -> Result<()> {
        let read_id = record.ref_id().unwrap();
        let read = record.ref_seq();
        
        let matches = self.demux_read(read, &self.query_groups);
        
        let total_count = self.total.fetch_add(1, Ordering::Relaxed);
        
        let mut writer = self.writer.lock().unwrap();
        
        if !matches.is_empty() {
            let found_count = self.found.fetch_add(1, Ordering::Relaxed);
            
            for m in matches {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    read_id,
                    m.match_str.stringify(),
                    m.start,
                    m.end,
                    m.edits,
                    m.rel_dist_to_end,
                    read.len(),
                    record_set_idx,
                    record_idx
                )?;
            }
        } else {
            let missed_count = self.missed.fetch_add(1, Ordering::Relaxed);
        }
        
        // Update progress bars
        if total_count % 100 == 0 {  // Update every 100 reads to avoid too frequent updates
            self.total_bar.set_message(format!("{}", total_count));
            self.found_bar.set_message(format!("{}", self.found.load(Ordering::Relaxed)));
            self.missed_bar.set_message(format!("{}", self.missed.load(Ordering::Relaxed)));
        }
        
        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<()> {
        Ok(())
    }
}