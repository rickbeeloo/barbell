use crate::annotate::distribution::*;
use crate::progress::ProgressTracker;
use colored::*;
use indicatif::ProgressBar;
use needletail::{FastxReader, Sequence, parse_fastx_file};
use rand::Rng;
use rand::RngCore;
use rayon::prelude::*;
use sassy::Searcher;
use sassy::profiles::Iupac;
use std::arch::x86_64::*;
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
        let start = std::time::Instant::now();
        let random_seq = random_dna_seq(min_len, max_len);
        let end = std::time::Instant::now();

        let start = std::time::Instant::now();
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

fn tune_single_seq_fast(seq: &[u8], n_iter: usize, fp_target: f32, alpha: f32) -> usize {
    let mut searcher = Searcher::<Iupac>::new_rc_with_overhang(alpha);

    // Pre-allocate a buffer large enough for the largest random sequence we will generate.
    // We always round the capacity up to the next multiple of 64 so we can use 64-byte SIMD
    // chunks when writing.
    let max_len = seq.len() * 2;
    let cap = max_len.next_multiple_of(64); // round-up to multiple of 64
    let mut buffer: Vec<u8> = vec![0; cap];

    let mut rng = rand::thread_rng();
    let mut costs = Vec::with_capacity(n_iter);

    // let pb = ProgressBar::new(n_iter as u64);
    for _ in 0..n_iter {
        // Determine the length of the random sequence for this iteration.
        //   • Always ≥64 so that at least one full 64-byte SIMD block is produced.
        //   • Always a multiple of 64 so the SIMD writer can run branch-free.
        let len = rng.gen_range(seq.len()..=max_len).next_multiple_of(64);

        // Fill the first `len` bytes of the buffer with a random DNA sequence.
        let start = std::time::Instant::now();
        fill_random_dna_simd(&mut buffer[..len]);
        let end = std::time::Instant::now();
        // println!("Time taken: {:?}", end.duration_since(start));

        // Search and record the lowest alignment cost (if any).
        let buffer_slice: &[u8] = &buffer[..len];
        let matches = searcher.search(seq, &buffer_slice, seq.len());
        if !matches.is_empty() {
            let lowest_cost = matches.iter().map(|m| m.cost).min().unwrap();
            costs.push(lowest_cost);
        }
        // pb.inc(1);
    }

    // If we never observed a hit, fall back to returning the sequence length
    if costs.is_empty() {
        return seq.len();
    }

    let cutoff = get_fp_threshold(costs, fp_target, TailSide::Left);
    cutoff as usize
}

/// Fill `dst` with a random DNA sequence consisting of the bytes `A`, `C`, `T` or `G`.
///
/// The length **must** be a multiple of 64 so we can operate in fixed-width SIMD blocks
/// without worrying about leftovers.  The implementation is fully portable – it relies on
/// writing 64-byte blocks at a time and uses bit-tricks instead of explicit SIMD intrinsics
/// so it compiles on any architecture, while still allowing the compiler to autovectorise.
#[inline(always)]
fn fill_random_dna_simd(dst: &mut [u8]) {
    debug_assert!(dst.len() % 64 == 0, "Length must be a multiple of 64");

    // Prepare our constant vectors:
    //  - mask: keep only the low 2 bits of each input byte
    //  - lut: a 16‑byte shuffle slot repeated four times → holds A,C,T,G pattern
    let mask = unsafe { _mm256_set1_epi8(0b11) };
    let lut = unsafe {
        _mm256_setr_epi8(
            b'A' as i8, b'C' as i8, b'T' as i8, b'G' as i8, b'A' as i8, b'C' as i8, b'T' as i8,
            b'G' as i8, b'A' as i8, b'C' as i8, b'T' as i8, b'G' as i8, b'A' as i8, b'C' as i8,
            b'T' as i8, b'G' as i8, // repeat for the high half
            b'A' as i8, b'C' as i8, b'T' as i8, b'G' as i8, b'A' as i8, b'C' as i8, b'T' as i8,
            b'G' as i8, b'A' as i8, b'C' as i8, b'T' as i8, b'G' as i8, b'A' as i8, b'C' as i8,
            b'T' as i8, b'G' as i8,
        )
    };

    let mut rng = rand::thread_rng();
    let mut rnd_block = [0u8; 64];

    for chunk in dst.chunks_exact_mut(64) {
        // 1) Fill 64 raw random bytes in one go
        rng.fill_bytes(&mut rnd_block);

        // 2) Process in two 32‑byte halves with AVX2
        for (off, ptr) in [(0, chunk.as_mut_ptr()), (32, chunk[32..].as_mut_ptr())] {
            unsafe {
                // load 32 random bytes
                let v = _mm256_loadu_si256(rnd_block[off..].as_ptr() as *const __m256i);
                // mask to 2 bits each
                let idx = _mm256_and_si256(v, mask);
                // shuffle LUT by those 2‑bit indices
                let dna = _mm256_shuffle_epi8(lut, idx);
                // store 32 DNA bytes
                _mm256_storeu_si256(ptr as *mut __m256i, dna);
            }
        }
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn random_tune() {
        let seq = b"GCTTGGGTGTTTAACCNNNNNNNNNNNNNNNNNNNNNNNNGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";
        let f = 0.0;
        let start_time = std::time::Instant::now();
        let cut_off = tune_single_seq_fast(seq, 10_000, f, 0.5);
        let end = std::time::Instant::now();
        println!("Time taken: {:?}", end.duration_since(start_time));
        println!("Cut off 1: {}", cut_off);

        // Check previous impl
        let start_time = std::time::Instant::now();
        let cut_off = tune_single_sequence(seq, 10_000, f, 0.5);
        let end = std::time::Instant::now();
        println!("Time taken: {:?}", end.duration_since(start_time));
        println!("Cutoff: {}", cut_off);
    }
}
