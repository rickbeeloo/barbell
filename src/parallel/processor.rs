use crate::annotate::search::BarMan;
use crate::annotate::merge_sort::merge_sort_files;

use seq_io::fastq::{Reader as FastqReader};
use seq_io_parallel::{MinimalRefRecord, ParallelProcessor, ParallelReader};
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::io::{Write, BufWriter};
use std::time::Instant;
use indicatif::{ProgressBar, ProgressStyle, MultiProgress};
use colored::Colorize;
use anyhow::Result;

#[derive(Clone)]
pub struct ParallelAnnotator {
    bar_man: Arc<BarMan>,
    writer: Arc<Mutex<BufWriter<std::fs::File>>>,
    total: Arc<AtomicUsize>,
    found: Arc<AtomicUsize>,
    missed: Arc<AtomicUsize>,
    total_bar: Arc<ProgressBar>,
    found_bar: Arc<ProgressBar>,
    missed_bar: Arc<ProgressBar>,
}

impl ParallelAnnotator {
    pub fn new(bar_man: BarMan) -> Self {
        let multi_progress = MultiProgress::new();
        let total_bar = multi_progress.add(ProgressBar::new_spinner());
        let found_bar = multi_progress.add(ProgressBar::new_spinner());
        let missed_bar = multi_progress.add(ProgressBar::new_spinner());

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
            bar_man: Arc::new(bar_man),
            writer: Arc::new(Mutex::new(BufWriter::new(std::fs::File::create("/dev/null").unwrap()))),
            total: Arc::new(AtomicUsize::new(0)),
            found: Arc::new(AtomicUsize::new(0)),
            missed: Arc::new(AtomicUsize::new(0)),
            total_bar: Arc::new(total_bar),
            found_bar: Arc::new(found_bar),
            missed_bar: Arc::new(missed_bar),
        }
    }

    pub fn process_fastq(&mut self, fastq_file: &str, output_file: &str, threads: usize) -> Result<()> {
        let start_time = Instant::now();

        println!("\n{}", "Configuration".bold().underline());
        println!("  • Input:  {}", fastq_file.bold());
        println!("  • Output: {}", output_file.bold());
        println!("  • Threads: {}", threads.to_string().bold());

        // Create output file and initialize writer
        let output_handle = std::fs::File::create(output_file)?;
        let writer = BufWriter::new(output_handle);
        self.writer = Arc::new(Mutex::new(writer));

        // Write header
        let mut writer = self.writer.lock().unwrap();
        writeln!(writer, "read\tlabel\tstart\tend\tlog_prob\tedit_dist\tread_len\trel_dist_to_end\trecord_set_idx\trecord_idx\tcuts")?;
        drop(writer);

        // Process reads
        let reader = FastqReader::from_path(fastq_file)?;
        reader.process_parallel(self.clone(), threads)?;

        // Finish progress bars
        self.total_bar.finish_with_message("Done!");
        self.found_bar.finish_with_message("Done!");
        self.missed_bar.finish_with_message("Done!");

        println!("\n{}", "Summary".bold().underline());
        println!("  • Time: {} seconds", start_time.elapsed().as_secs().to_string().bold());
        println!("  • Total reads: {}", self.total.load(Ordering::Relaxed).to_string().bold());
        println!("  • Tagged reads: {}", self.found.load(Ordering::Relaxed).to_string().green().bold());
        println!("  • Missed reads: {}\n", self.missed.load(Ordering::Relaxed).to_string().red().bold());

        // Flush before closing
        self.writer.lock().unwrap().flush()?;

        // Merge sort output file
        println!("\n{}", "Merge sort".bold().underline());
        let start_time = Instant::now();
        merge_sort_files(output_file);
        println!("  • Time: {} seconds", start_time.elapsed().as_secs().to_string().bold());

        Ok(())
    }
}

impl ParallelProcessor for ParallelAnnotator {
    fn process_record<'a, Rf: MinimalRefRecord<'a>>(&mut self, record: Rf, record_set_idx: usize, record_idx: usize) -> Result<()> {
        // Parse the read ID and sequence
        let read_id = record.ref_id().unwrap().split_whitespace().next().unwrap();
        let read = record.ref_seq();
        
        // Annotate the read
        let matches = self.bar_man.annotate(read);
        
        // Update total counter
        let total_count = self.total.fetch_add(1, Ordering::Relaxed);
        
        // Get the writer lock FIRST, exactly like in your old code
        let mut writer = self.writer.lock().unwrap();
        
        if !matches.is_empty() {
            // Update found counter
            self.found.fetch_add(1, Ordering::Relaxed);
            
            // Directly write each match using writeln! just like your old code
            for m in &matches {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t-",
                    read_id,
                    m.label.stringify(),
                    m.start,
                    m.end,
                    m.log_prob.unwrap_or(0.0),
                    m.edit_dist.unwrap_or(0),
                    read.len(),
                    m.rel_dist_to_end,
                    record_set_idx,
                    record_idx
                )?;
            }
        } else {
            self.missed.fetch_add(1, Ordering::Relaxed);
        }
        
        // Release the lock when done
        drop(writer);
        
        // Update progress bars periodically
        if total_count % 100 == 0 {
            self.total_bar.set_message(format!("{}", total_count));
            self.found_bar.set_message(format!("{}", self.found.load(Ordering::Relaxed)));
            self.missed_bar.set_message(format!("{}", self.missed.load(Ordering::Relaxed)));
        }
        
        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<()> {
        // Optionally flush after each batch
        // self.writer.lock().unwrap().flush()?;
        Ok(())
    }
} 