use crate::progress::ProgressTracker;
use crate::search::distribution::*;
use colored::*;
use needletail::{FastxReader, Sequence, parse_fastx_file};
use rand::Rng;
use rayon::prelude::*;
use sassy::profiles::Iupac;
use sassy::search::Searcher;
use std::thread_local;

thread_local! {
    static TUNING_SEARCHER: std::cell::RefCell<Option<Searcher<Iupac>>> = std::cell::RefCell::new(None);
}

pub fn tune_single_sequence(seq: &[u8], n_iter: usize, fp_target: f32, alpha: f32) -> usize {
    let min_len = seq.len();
    let max_len = min_len + min_len;
    let mut costs = Vec::new();

    // Initialize thread-local searcher if not already done
    TUNING_SEARCHER.with(|cell| {
        if cell.borrow().is_none() {
            *cell.borrow_mut() = Some(Searcher::<Iupac>::new_rc_with_overhang(alpha));
        }
    });

    for _ in 0..n_iter {
        let random_seq = random_dna_seq(min_len, max_len);
        let matches = TUNING_SEARCHER.with(|cell| {
            if let Some(ref mut searcher) = *cell.borrow_mut() {
                searcher.search(seq, &random_seq, seq.len())
            } else {
                vec![]
            }
        });

        if matches.is_empty() {
            continue;
        }
        let lowest_cost = matches.iter().map(|m| m.cost).min().unwrap();
        costs.push(lowest_cost);
    }
    let cutoff = get_fp_threshold(costs, fp_target, TailSide::Left);
    cutoff as usize
}

pub fn random_dna_seq(min_len: usize, max_len: usize) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let len = rng.gen_range(min_len..max_len);
    let mut seq = Vec::new();
    for _ in 0..len {
        let base = match rng.gen_range(0..4) {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };
        seq.push(base);
    }
    seq
}

/// Tune multiple sequences in parallel and return their k-cutoffs
pub fn tune_sequences_parallel(
    sequences: &[&[u8]],
    n_iter: usize,
    fp_target: f32,
    alpha: f32,
    n_threads: usize,
) -> Vec<usize> {
    let mut progress = ProgressTracker::new();
    progress.step("Running parallel tuning with random sequences");
    progress.indent();
    progress.substep(&format!("Using {} threads", n_threads));
    progress.substep(&format!("Target false positive rate: {:.6}", fp_target));
    progress.substep(&format!("Random iterations per sequence: {}", n_iter));

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .unwrap();

    let cutoffs = pool.install(|| {
        sequences
            .par_iter()
            .map(|seq| tune_single_sequence(seq, n_iter, fp_target, alpha))
            .collect()
    });

    progress.dedent();
    progress.success("Parallel tuning completed");

    cutoffs
}

fn tune_single_sequence_with_real_data(
    seq: &[u8],
    real_sequences: &[Vec<u8>],
    fp_target: f32,
    alpha: f32,
) -> usize {
    let mut costs = Vec::new();

    // Initialize thread-local searcher if not already done
    TUNING_SEARCHER.with(|cell| {
        if cell.borrow().is_none() {
            *cell.borrow_mut() = Some(Searcher::<Iupac>::new_rc_with_overhang(alpha));
        }
    });

    for real_seq in real_sequences {
        // Skip if the real sequence is shorter than our query sequence
        if real_seq.len() < seq.len() {
            continue;
        }

        let matches = TUNING_SEARCHER.with(|cell| {
            if let Some(ref mut searcher) = *cell.borrow_mut() {
                searcher.search(seq, real_seq, seq.len())
            } else {
                vec![]
            }
        });

        if matches.is_empty() {
            continue;
        }
        let lowest_cost = matches.iter().map(|m| m.cost).min().unwrap();
        costs.push(lowest_cost);
    }

    if costs.is_empty() {
        return seq.len(); // Default to sequence length if no matches
    }

    let cutoff = get_fp_threshold(costs, fp_target, TailSide::Left);
    cutoff as usize
}

/// Tune multiple sequences in parallel using real sequences from FASTA files
pub fn tune_sequences_parallel_with_fasta(
    sequences: &[&[u8]],
    fasta_files: &[&str],
    fp_target: f32,
    alpha: f32,
    n_threads: usize,
    max_tune_seqs: usize,
) -> Vec<usize> {
    let mut progress = ProgressTracker::new();
    progress.step("Collecting real sequences from FASTA files");
    progress.indent();

    // The longest sequence to tune is the flank, so we should remove **at least** that many
    // bases, lets add 20% to that to make sure we strip it off
    let max_len = sequences.iter().map(|s| s.len()).max().unwrap();
    let flank_len = (max_len as f32 * 1.5) as usize;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .unwrap();

    // Collect all middle portions from FASTA files
    let mut middle_sequences: Vec<Vec<u8>> = Vec::new();

    for (i, fasta_file) in fasta_files.iter().enumerate() {
        progress.substep(&format!("Processing FASTA file {}: {}", i + 1, fasta_file));
        let mut reader = parse_fastx_file(fasta_file).expect("Failed to parse FASTA file");

        let mut file_count = 0;
        while let Some(Ok(record)) = reader.next() {
            let seq = record.seq();
            if seq.len() > flank_len * 2 {
                // Only use sequences longer than 400nt
                let end = seq.len().saturating_sub(flank_len);
                if end > flank_len {
                    middle_sequences.push(seq[flank_len..end].to_vec());
                    file_count += 1;
                }
                if middle_sequences.len() >= max_tune_seqs {
                    break;
                }
            }
        }
        progress.info(&format!(
            "Extracted {} sequences from {}",
            file_count, fasta_file
        ));
    }

    progress.substep(&format!(
        "Total real sequences collected: {}",
        middle_sequences.len()
    ));
    progress.dedent();

    progress.step("Running parallel tuning");
    progress.indent();
    progress.substep(&format!("Using {} threads", n_threads));
    progress.substep(&format!("Target false positive rate: {:.6}", fp_target));

    let cutoffs = pool.install(|| {
        sequences
            .par_iter()
            .map(|seq| {
                tune_single_sequence_with_real_data(seq, &middle_sequences, fp_target, alpha)
            })
            .collect()
    });

    progress.dedent();
    progress.success("Parallel tuning completed");

    cutoffs
}
