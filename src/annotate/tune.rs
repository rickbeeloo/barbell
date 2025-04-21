use crate::annotate::flank::*;
use rand::Rng;
use spinners::{Spinner, Spinners};
use colored::Colorize;
// use crate::annotate::mutations::*;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use pa_bitpacking::search::*;
use std::cmp::Ordering::Equal;
use indicatif::{ProgressBar, ProgressStyle};
use std::sync::atomic::{AtomicUsize, Ordering};
use crate::simdedits::simd::{TransposedQueries, simd_search};

pub fn generate_random_sequence(length: usize) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let bases = b"ACGT";
    (0..length).map(|_| bases[rng.gen_range(0..4)]).collect()
}


pub trait Number: Copy + PartialOrd + std::fmt::Debug {}

// Implement the trait for f64 and i32
impl Number for f64 {}
impl Number for i32 {}
impl Number for u8 {}

#[derive(PartialEq, Eq)]
pub enum TailSide {
    Left,
    Right,
}


fn get_fp_threshold<T: Number>(scores: Vec<T>, fp_target: f64, tail_side: TailSide) -> T {
    let mut sorted_scores = scores.clone();
    
    // Sort based on tail side
    if tail_side == TailSide::Left {
        // Sort from lowest to highest
        sorted_scores.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Equal));
    } else {
        // Sort from highest to lowest
        sorted_scores.sort_by(|a, b| b.partial_cmp(a).unwrap_or(Equal));
    }

    let target_index = (fp_target * sorted_scores.len() as f64) as usize;

    println!("Target index: {}", target_index);
    // println!("First 10 sorted scores: {:?}", &sorted_scores[0..10]);
    println!("Score length: {}", sorted_scores.len());

    if target_index < sorted_scores.len() {
        sorted_scores[target_index]
    } else {
        // Print yellow warning that we don't have enough data for this fp target
        println!("{}", "Warning: Not enough data for this fp target. Using max cut off (uncertain).".yellow());
        sorted_scores[sorted_scores.len() - 1]
    }
}

pub fn tune_edit_distance(flanks: &[FlankGroup], n: usize, fp_target: f64) -> Vec<i32> {
    let mut sp = Spinner::new(Spinners::OrangeBluePulse, "Tuning edit distance".into());

    let seq_len = if let Some(flank) = flanks.first() {
        flank.flank_seq.seq.len() * 2
    } else {
        panic!("{}", "Flanks are empty!".red());
    };

    // Collect queries from all flanks
    let queries: Vec<_> = flanks.iter().map(|f| f.flank_seq.seq.as_ref()).collect();

    // Use Rayon to chunk the workload and generate per-thread results
    let all_edit_distances: Vec<Vec<i32>> = queries.par_iter().map(|query| {
        (0..n)
            .into_par_iter()
            .map(|_| {
                let random_seq = generate_random_sequence(seq_len);
                let result = search(query, &random_seq, 0.4);
                let edit_dist = result.out.iter().min().unwrap_or(&i32::MAX);
                *edit_dist
            })
            .collect::<Vec<i32>>()
    }).collect();

    sp.stop();

    println!("Input sequence length: {}", seq_len);

    // Calculate the cutoff for each query
    all_edit_distances.into_iter().map(|distances| {
        get_fp_threshold(distances, fp_target, TailSide::Left)
    }).collect()
}



pub fn tune_max_edits(flanks: &Vec<FlankGroup>, n: usize, fp_target: f64) -> u8 {
    let mut sp = Spinner::new(Spinners::OrangeBluePulse, "Tuning max edits".into());

    // For now lets just test with the first flank
    let mask_query_slices = flanks[0]
        .flank_seq
        .mask_queries
        .iter()
        .map(|q| q.as_ref())
        .collect::<Vec<_>>();

    let mut mask_transposed_queries = Vec::with_capacity(mask_query_slices.len() / 32);
    for group in mask_query_slices.chunks(32) {
        mask_transposed_queries.push(TransposedQueries::new(group.to_vec()));
    }

    let min_diff: Vec<u8> = (0..n).into_par_iter()
    .map(|_| {
        let random_seq = generate_random_sequence(mask_query_slices[0].len() + 10);
        let mut edits = mask_transposed_queries
            .iter()
            .flat_map(|tq| simd_search(tq, &random_seq))
            .collect::<Vec<_>>();
        
        edits.sort_unstable(); // Sort to get the top two
        if edits.len() >= 2 {
            edits[1].saturating_sub(edits[0]) // difference between top 2
        } else {
            0 // fallback in case there's only one result
        }
    })
    .collect();

    let min_edit_diff = get_fp_threshold(min_diff, fp_target, TailSide::Right);
    sp.stop();

    min_edit_diff

    // // This loop can be parallelized
    // let min_edits: Vec<u8> = (0..n).into_par_iter()
    //     .map(|_| {
    //         let random_seq = generate_random_sequence(mask_query_slices[0].len() + 10);
    //         let edits = mask_transposed_queries
    //             .iter()
    //             .flat_map(|tq| simd_search(tq, &random_seq))
    //             .collect::<Vec<_>>();
    //         *edits.iter().min().unwrap()
    //     })
    //     .collect();

    // let min_edit = get_fp_threshold(min_edits, fp_target, TailSide::Left);
    // sp.stop();

    // min_edit
}

// pub fn tune_log_prob(flanks: &Vec<FlankGroup>, n: usize, fp_target: f64, error_rates: &ErrorRatesAffine) -> (f64, f64) {
//     let log_probs = Arc::new(Mutex::new(Vec::with_capacity(n)));
//     let top_two_pairs = Arc::new(Mutex::new(Vec::with_capacity(n)));

//     let per_flank = n / flanks.len().max(1); // avoid div by zero
//     let counter = Arc::new(AtomicUsize::new(0));

//     let pb = ProgressBar::new(n as u64);
//     pb.set_style(
//         ProgressStyle::default_bar()
//             .template("{spinner:.blue} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
//             .unwrap()
//             .progress_chars("##-"),
//     );

//     flanks.par_iter().for_each(|flank_group| {
//         let mut local_probs = Vec::with_capacity(per_flank);
//         let mut local_pairs = Vec::with_capacity(per_flank);

//         let queries: Vec<&[u8]> = flank_group.flank_seq.mask_queries.iter().map(|q| q.as_ref()).collect();
//         if queries.is_empty() {
//             return;
//         }

//         for _ in 0..per_flank {
//             let query_len = queries[0].len();
//             let random_seq = generate_random_sequence(query_len + 10);
//             let probs = optimal_metric(&queries, &random_seq, error_rates);

//             if probs.len() >= 2 {
//                 let mut sorted = probs.clone();
//                 sorted.sort_by(|a, b| b.partial_cmp(a).unwrap_or(Equal));
//                 let best = sorted[0];
//                 let second_best = sorted[1];

//                 local_probs.push(best);
//                 local_pairs.push((best, second_best));
//             }

//             let count = counter.fetch_add(1, Ordering::Relaxed);
//             if count < n {
//                 pb.inc(1);
//             }
//         }

//         log_probs.lock().unwrap().extend(local_probs);
//         top_two_pairs.lock().unwrap().extend(local_pairs);
//     });

//     pb.finish_with_message("Done");

//     let log_probs = Arc::try_unwrap(log_probs).unwrap().into_inner().unwrap();
//     let top_two_pairs = Arc::try_unwrap(top_two_pairs).unwrap().into_inner().unwrap();

//     let min_log_t = get_fp_threshold(log_probs.clone(), fp_target, TailSide::Right);
//     let diffs: Vec<f64> = top_two_pairs.into_iter()
//         .filter_map(|(best, second)| if best > min_log_t { Some(best - second) } else { None })
//         .collect();

//     let min_d_log = get_fp_threshold(diffs, fp_target, TailSide::Right);

//     println!("Minimum log probability: {}", min_log_t);
//     println!("Minimum log distance: {}", min_d_log);
//     (min_log_t, min_d_log)
// }

