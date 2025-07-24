use colored::*;
use std::time::Instant;

pub struct ProgressTracker {
    start_time: Instant,
    indent_level: usize,
}

impl ProgressTracker {
    pub fn new() -> Self {
        Self {
            start_time: Instant::now(),
            indent_level: 0,
        }
    }

    pub fn step(&mut self, message: &str) {
        let indent = "  ".repeat(self.indent_level);
        println!("{}{} {}", indent, "•".blue(), message.bold());
    }

    pub fn substep(&mut self, message: &str) {
        let indent = "  ".repeat(self.indent_level + 1);
        println!("{}{} {}", indent, "◦".cyan(), message);
    }

    pub fn info(&mut self, message: &str) {
        let indent = "  ".repeat(self.indent_level + 2);
        println!("{}{} {}", indent, "→".yellow(), message.dimmed());
    }

    pub fn success(&mut self, message: &str) {
        let indent = "  ".repeat(self.indent_level);
        println!("{}{} {}", indent, "✓".green(), message.green().bold());
    }

    pub fn warning(&mut self, message: &str) {
        let indent = "  ".repeat(self.indent_level);
        println!("{}{} {}", indent, "⚠".yellow(), message.yellow());
    }

    pub fn error(&mut self, message: &str) {
        let indent = "  ".repeat(self.indent_level);
        println!("{}{} {}", indent, "✗".red(), message.red().bold());
    }

    pub fn indent(&mut self) {
        self.indent_level += 1;
    }

    pub fn dedent(&mut self) {
        if self.indent_level > 0 {
            self.indent_level -= 1;
        }
    }

    pub fn elapsed(&self) -> std::time::Duration {
        self.start_time.elapsed()
    }

    pub fn print_elapsed(&self) {
        let elapsed = self.elapsed();
        let indent = "  ".repeat(self.indent_level);
        println!(
            "{}{} Completed in {:.2}s",
            indent,
            "⏱".blue(),
            elapsed.as_secs_f64()
        );
    }
}

pub fn print_header(title: &str) {
    println!("\n{}", "=".repeat(60).blue());
    println!("{}", format!("  {}", title).blue().bold());
    println!("{}", "=".repeat(60).blue());
}

pub fn print_barcode_group_info(
    flank: &[u8],
    barcodes: &[crate::annotate::barcodes::Barcode],
    cutoffs: &[usize],
) {
    let indent = "  ".repeat(2);

    // Print flank information
    println!("{}{} {}", indent, "🔍".blue(), "Flank sequence:".bold());
    let flank_str = String::from_utf8_lossy(flank);
    let flank_display = if flank_str.len() > 80 {
        format!(
            "{}...{}",
            &flank_str[..40],
            &flank_str[flank_str.len() - 40..]
        )
    } else {
        flank_str.to_string()
    };
    println!("{}{} {}", indent, "  ".repeat(1), flank_display.cyan());
    println!(
        "{}{} Cutoff: {}",
        indent,
        "  ".repeat(1),
        cutoffs[0].to_string().green()
    );

    // Print barcode information (show first 2 and last 2)
    println!("{}{} {}", indent, "🧬".blue(), "Barcode sequences:".bold());

    let total_barcodes = barcodes.len();
    let show_count = total_barcodes.min(4);

    for i in 0..show_count {
        let barcode = &barcodes[i];
        let cutoff = cutoffs[i + 1]; // +1 because first cutoff is for flank

        if i == 2 && total_barcodes > 4 {
            println!(
                "{}{} ... and {} more barcodes",
                indent,
                "  ".repeat(1),
                total_barcodes - 4
            );
            break;
        }

        let barcode_str = String::from_utf8_lossy(&barcode.seq);
        println!(
            "{}{} {}: {} (cutoff: {})",
            indent,
            "  ".repeat(1),
            barcode.label.cyan(),
            barcode_str.yellow(),
            cutoff.to_string().green()
        );
    }
}

pub fn print_tuning_progress(sequences: &[&[u8]], fasta_files: &[&str], max_tune_seqs: usize) {
    let indent = "  ".repeat(2);

    println!(
        "{}{} {}",
        indent,
        "📊".blue(),
        "Tuning configuration:".bold()
    );
    println!(
        "{}{} Sequences to tune: {}",
        indent,
        "  ".repeat(1),
        sequences.len()
    );
    println!(
        "{}{} FASTA files: {}",
        indent,
        "  ".repeat(1),
        fasta_files.len()
    );
    println!(
        "{}{} Max tuning sequences: {}",
        indent,
        "  ".repeat(1),
        max_tune_seqs
    );
}

pub fn print_cutoffs_horizontal(cutoffs: &[usize], labels: &[&str]) {
    let indent = "  ".repeat(2);
    println!("{}{} {}", indent, "🎯".blue(), "Optimized cutoffs:".bold());

    let cutoff_str = cutoffs
        .iter()
        .enumerate()
        .map(|(i, &cutoff)| {
            if i < labels.len() {
                format!("{}: {}", labels[i].cyan(), cutoff.to_string().green())
            } else {
                format!("seq{}: {}", i, cutoff.to_string().green())
            }
        })
        .collect::<Vec<_>>()
        .join(" | ");

    println!("{}{} {}", indent, "  ".repeat(1), cutoff_str);
}

pub fn print_summary_stats(
    total_reads: usize,
    mapped_reads: usize,
    trimmed_reads: usize,
    trimmed_split_reads: usize,
    elapsed: std::time::Duration,
) {
    let indent = "  ".repeat(1);

    println!("{}{} {}", indent, "📈".blue(), "Summary:".bold());
    println!(
        "{}{} Total reads: {}",
        indent,
        "  ".repeat(1),
        total_reads.to_string().bold()
    );
    println!(
        "{}{} Mapped reads: {} ({:.1}%)",
        indent,
        "  ".repeat(1),
        mapped_reads.to_string().green().bold(),
        (mapped_reads as f64 / total_reads as f64) * 100.0
    );
    println!(
        "{}{} Trimmed reads: {} ({:.1}%)",
        indent,
        "  ".repeat(1),
        trimmed_reads.to_string().green().bold(),
        (trimmed_reads as f64 / mapped_reads as f64) * 100.0
    );
    println!(
        "{}{} Extra splits: {}",
        indent,
        "  ".repeat(1),
        (trimmed_split_reads - trimmed_reads)
            .to_string()
            .blue()
            .bold()
    );
    println!(
        "{}{} Time: {:.2}s",
        indent,
        "  ".repeat(1),
        elapsed.as_secs_f64().to_string().bold()
    );
}
