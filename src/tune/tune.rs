use crate::annotate::barcodes::*;
use rand::{Rng, thread_rng};
use sassy::*;
use std::cmp::min;

/// Number of Monte Carlo samples to draw (increase for more precision).
const NUM_SAMPLES: usize = 1_000_000;
/// The alphabet to draw from.
const ALPHABET: &[u8] = b"ACGT";

/// Compute Levenshtein edit distance between two equal-length byte slices.
fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let len = a.len();
    let mut prev = (0..=len).collect::<Vec<_>>();
    let mut curr = vec![0; len + 1];

    for i in 0..len {
        curr[0] = i + 1;
        for j in 0..len {
            let cost = if a[i] == b[j] { 0 } else { 1 };
            curr[j + 1] = min(
                min(
                    prev[j + 1] + 1, // deletion
                    curr[j] + 1,
                ), // insertion
                prev[j] + cost, // substitution
            );
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[len]
}

/// Simulate `NUM_SAMPLES` random pairs of length-`n` sequences,
/// then return the edit-distance cutoff at quantile `conf_val`.
///
/// # Arguments
/// * `n` — length of each random sequence
/// * `conf_val` — desired confidence level in (0.0, 1.0), e.g. 0.999 for 99.9%
///
/// # Returns
/// * An integer edit-distance such that approximately `conf_val * 100%`
///   of simulated distances lie at or below it.
pub fn sim_edits(n: usize, conf_val: f64) -> usize {
    assert!((0.0..1.0).contains(&conf_val), "conf_val must be in (0,1)");

    let mut rng = thread_rng();
    let mut dists = Vec::with_capacity(NUM_SAMPLES);

    for _ in 0..NUM_SAMPLES {
        // generate two random sequences of length n
        let s1: Vec<u8> = (0..n)
            .map(|_| ALPHABET[rng.gen_range(0..ALPHABET.len())])
            .collect();
        let s2: Vec<u8> = (0..n)
            .map(|_| ALPHABET[rng.gen_range(0..ALPHABET.len())])
            .collect();
        let d = edit_distance(&s1, &s2);
        // println!("d: {}", d);
        dists.push(d);
    }

    dists.sort_unstable();
    let idx = ((NUM_SAMPLES as f64) * conf_val).ceil() as usize - 1;
    dists[idx.min(NUM_SAMPLES - 1)]
}

pub fn tune(
    _fastq_file: &str,
    query_file: &str,
    _max_fp_rate_cli: Option<f64>,
    _fp_budget_cli: Option<u64>,
    confidence_cli: Option<f64>,
) {
    // Create query group to get lengths
    let barcode_group = BarcodeGroup::new_from_fasta(query_file, BarcodeType::Ftag);

    // Get flank length (excluding N's)
    let flank_len = barcode_group.flank.iter().filter(|&&c| c != b'N').count();

    // Get barcode length (assuming all barcodes have same length)
    let barcode_len = barcode_group.barcodes[0].seq.len();

    // Use confidence level (default to 0.999 for 99.9%)
    let confidence = confidence_cli.unwrap_or(0.9999);

    println!("Simulating edit distances for:");
    println!("  Flank length (excluding N's): {}", flank_len);
    println!("  Barcode length: {}", barcode_len);
    println!("  Confidence level: {:.1}%", confidence * 100.0);
    println!();

    // Simulate flank cutoffs
    println!("Simulating flank edit distances...");
    let flank_cutoff = sim_edits(flank_len, confidence);
    println!(
        "Flank cutoff at {:.1}% confidence: {}",
        confidence * 100.0,
        flank_cutoff
    );

    // Simulate barcode cutoffs
    println!("Simulating barcode edit distances...");
    let barcode_cutoff = sim_edits(barcode_len, confidence);
    println!(
        "Barcode cutoff at {:.1}% confidence: {}",
        confidence * 100.0,
        barcode_cutoff
    );

    println!();
    println!("Recommended settings:");
    println!("  --flank-max-errors {}", flank_cutoff);
    println!("  --barcode-max-errors {}", barcode_cutoff);
}

// Backwards compatibility for existing CLI: call with default FP settings
pub fn tune_simple(fastq_file: &str, query_file: &str) {
    tune(fastq_file, query_file, None, None, None);
}
