use crate::WIGGLE_ROOM;
use crate::annotate::barcodes::*;
use indicatif::{ProgressBar, ProgressStyle};
use needletail::{Sequence, parse_fastx_file};
use plotpy::{Histogram, Plot, StrError};
use rand::Rng;
use rayon::prelude::*;
use sassy::{profiles::Iupac, *};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;

#[derive(Clone, Debug, PartialEq, clap::ValueEnum)]
pub enum TargetSide {
    Left,
    Right,
}

fn random_dna_seq(length: usize) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let mut seq = Vec::with_capacity(length);
    let bases = b"ACGT";
    for _ in 0..length {
        seq.push(bases[rng.gen_range(0..bases.len())]);
    }
    seq
}

fn get_min_cost(matches: &[Match], q_len: usize) -> i32 {
    matches.iter().map(|m| m.cost).min().unwrap_or(q_len as i32)
}

fn get_edits(
    searcher: &mut Searcher<Iupac>,
    window: &[u8],
    random_window: &[u8],
    query: &[u8],
) -> (i32, i32) {
    // We just get all edits, then the lowest
    let k = query.len();

    // Do the real window
    let matches = searcher.search_all(query, &window, k);
    let lowest_edits = get_min_cost(&matches, query.len());

    // Do the random window
    let random_matches = searcher.search_all(query, &random_window, k);
    let random_lowest_edits = get_min_cost(&random_matches, query.len());

    (lowest_edits, random_lowest_edits)
}

fn find_best_cutoff(real: &[i32], random: &[i32], prevalence: f64) -> i32 {
    assert!(
        (0.0..=1.0).contains(&prevalence),
        "prevalence must be in [0,1]"
    );

    let max_cost = *real
        .iter()
        .chain(random.iter())
        .max()
        .expect("No costs available");

    let real_total = real.len() as f64;
    let random_total = random.len() as f64;

    let mut best_cutoff = 0;
    let mut best_error = f64::INFINITY;

    for t in 0..=max_cost {
        // False-negative rate among reads that truly contain the pattern
        let fn_rate = real.iter().filter(|&&c| c > t).count() as f64 / real_total;
        // False-positive rate among reads that truly do NOT contain the pattern
        let fp_rate = random.iter().filter(|&&c| c <= t).count() as f64 / random_total;

        // Expected overall misclassification probability given prevalence
        let error_rate = prevalence * fn_rate + (1.0 - prevalence) * fp_rate;

        // println!("Cut-off {t}: weighted error {error_rate:.4}");

        if error_rate < best_error {
            best_error = error_rate;
            best_cutoff = t;
        }
    }

    best_cutoff as i32
}

pub fn tune(fastq_file: &str, query_file: &str, target_side: TargetSide) {
    // Store edits for flank and barcodes
    let mut real_flank_edits = Vec::new();
    let mut random_flank_edits = Vec::new();
    let mut real_bar_edits = Vec::new();
    let mut random_bar_edits = Vec::new();

    // Create query group
    let barcode_group = BarcodeGroup::new_from_fasta(query_file, BarcodeType::Ftag);
    // Then we can extract the flank queries and barcode queries
    let flank = barcode_group.flank;
    let barcodes = barcode_group.barcodes;

    // Reusable sassy searcher
    let mut searcher = Searcher::<Iupac>::new_rc_with_overhang(0.5);

    // Window size is based on flank length + some wiggle room
    let window_size = flank.len() + WIGGLE_ROOM;

    // Progress bar
    let pb = ProgressBar::new_spinner();
    pb.set_style(ProgressStyle::with_template("processed: {pos}").unwrap());

    // Then we go over the reads,
    let mut reader = parse_fastx_file(fastq_file).unwrap();
    while let Some(Ok(record)) = reader.next() {
        let seq = record.seq();
        let id = record.id();
        // Check if we have at least the window size of read length
        if seq.len() < window_size {
            continue;
        }
        // Then we can extract the window
        let window = if target_side == TargetSide::Left {
            seq[0..window_size].to_vec()
        } else {
            seq[seq.len() - window_size..].to_vec()
        };
        // Then we create an equally long random sequence
        let random_window = random_dna_seq(window_size);

        // Get flank edits for this read
        let (flank_real_edits, flank_rand_edits) =
            get_edits(&mut searcher, &window, &random_window, &flank);

        real_flank_edits.push(flank_real_edits);
        random_flank_edits.push(flank_rand_edits);

        // Parallel iterate over all barcodes and keep the lowest edit distance
        let (min_bar_real, min_bar_rand) = barcodes
            .par_iter()
            .map_init(
                || Searcher::<Iupac>::new_fwd_with_overhang(0.5),
                |s, barcode| {
                    let query = &barcode.seq;
                    let k = query.len();

                    let cost_real = {
                        let matches_real = s.search_all(query, &window, k);
                        get_min_cost(&matches_real, query.len())
                    };

                    let cost_rand = {
                        let matches_rand = s.search_all(query, &random_window, k);
                        get_min_cost(&matches_rand, query.len())
                    };

                    (cost_real, cost_rand)
                },
            )
            .reduce(|| (i32::MAX, i32::MAX), |a, b| (a.0.min(b.0), a.1.min(b.1)));

        real_bar_edits.push(min_bar_real);
        random_bar_edits.push(min_bar_rand);

        pb.inc(1);
    }

    pb.finish_and_clear();

    // ----- Write cost count files -----
    fn write_counts(path: &str, real: &[i32], random: &[i32]) -> std::io::Result<()> {
        let mut counts: HashMap<i32, (u32, u32)> = HashMap::new();
        for &v in real {
            counts.entry(v).or_insert((0, 0)).0 += 1;
        }
        for &v in random {
            counts.entry(v).or_insert((0, 0)).1 += 1;
        }
        // Write sorted by cost
        let mut keys: Vec<_> = counts.keys().cloned().collect();
        keys.sort();
        let mut file = File::create(path)?;
        writeln!(file, "cost\treal\trandom")?;
        for k in keys {
            let (r, rand) = counts[&k];
            writeln!(file, "{k}\t{r}\t{rand}")?;
        }
        Ok(())
    }

    if let Err(e) = write_counts("flank_costs.tsv", &real_flank_edits, &random_flank_edits) {
        eprintln!("Failed to write flank_costs.tsv: {e}");
    }
    if let Err(e) = write_counts("barcode_costs.tsv", &real_bar_edits, &random_bar_edits) {
        eprintln!("Failed to write barcode_costs.tsv: {e}");
    }

    // Prepare data: two series of values
    let values = vec![real_flank_edits.clone(), random_flank_edits.clone()];
    let labels = ["Real edits", "Random edits"];

    let mut histogram = Histogram::new();
    histogram
        .set_number_bins(30) // set number of bins
        .set_style("stepfilled") // or "bar", "step", etc.
        .set_stacked(false) // separate groups
        .set_colors(&["#1f77b4", "#ff7f0e"]) // custom colors
        .set_line_width(1.5);

    histogram.draw(&values, &labels);

    let mut plot = Plot::new();
    plot.add(&histogram)
        .grid_labels_legend("Edit Distance", "Frequency");
    plot.set_title("Edit Distance Distributions");

    plot.save("edit_hist.svg").expect("Plot save failed");
    plot.show("edit_hist.svg").unwrap();

    // Determine recommended cut-offs (assume ~80% of real reads contain a match)
    const PREVALENCE: f64 = 0.8;
    let flank_cutoff = find_best_cutoff(&real_flank_edits, &random_flank_edits, PREVALENCE);
    let barcode_cutoff = find_best_cutoff(&real_bar_edits, &random_bar_edits, PREVALENCE);

    println!("Recommended flank cut-off: {flank_cutoff}");
    println!("Recommended barcode cut-off: {barcode_cutoff}");
}
