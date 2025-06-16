use crate::annotate::flank::*;
use colored::Colorize;
use rand::Rng;
use spinners::{Spinner, Spinners};
// use crate::annotate::mutations::*;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use sassy::profiles::*;
use sassy::search::*;
use std::cmp::Ordering::Equal;
use std::thread_local;

thread_local! {
    static THREAD_SEARCHER: std::cell::RefCell<Option<Searcher<Iupac>>> = std::cell::RefCell::new(None);
}

pub fn generate_random_sequence(min_length: usize, max_length: usize) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let bases = b"ACGT";
    let length = rng.gen_range(min_length..max_length);
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
        println!(
            "{}",
            "Warning: Not enough data for this fp target. Using max cut off (uncertain).".yellow()
        );
        sorted_scores[sorted_scores.len() - 1]
    }
}

pub fn tune_edit_distance(
    searcher: &mut Searcher<Iupac>,
    flanks: &[FlankGroup],
    n: usize,
    fp_target: f64,
) -> Vec<i32> {
    let mut sp = Spinner::new(Spinners::OrangeBluePulse, "Tuning edit distance".into());

    let seq_len = if let Some(flank) = flanks.first() {
        flank.flank_seq.seq.len()
    } else {
        panic!("{}", "Flanks are empty!".red());
    };

    // Collect queries from all flanks
    let queries: Vec<_> = flanks.iter().map(|f| f.flank_seq.seq.as_ref()).collect();
    let half_flank_len = (seq_len as f32 / 2.0).round() as usize;

    // Create progress bar
    let pb = ProgressBar::new((queries.len() * n) as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.blue} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("##-"),
    );

    // Use Rayon to chunk the workload and generate per-thread results
    let all_edit_distances: Vec<Vec<i32>> = queries
        .par_iter()
        .map(|query| {
            (0..n)
                .into_par_iter()
                .map(|_| {
                    // Initialize thread-local searcher if not already done
                    THREAD_SEARCHER.with(|thread_searcher| {
                        if thread_searcher.borrow().is_none() {
                            *thread_searcher.borrow_mut() = Some(searcher.clone());
                        }
                    });

                    // Generate a random sequence
                    let random_seq = generate_random_sequence(seq_len / 2, seq_len * 2);
                    // Use Sassy search with thread-local searcher
                    let matches = THREAD_SEARCHER.with(|thread_searcher| {
                        thread_searcher.borrow_mut().as_mut().unwrap().search(
                            query,
                            &random_seq,
                            (query.len() as f32 / 2.0) as usize,
                        )
                    });
                    let edit_dist = matches.iter().min_by_key(|m| m.cost).unwrap().cost;
                    pb.inc(1);
                    edit_dist
                })
                .collect::<Vec<i32>>()
        })
        .collect();

    pb.finish_with_message("Done tuning edit distances");
    sp.stop();

    println!("Input sequence length: {}", seq_len);

    // Calculate the cutoff for each query
    all_edit_distances
        .into_iter()
        .map(|distances| get_fp_threshold(distances, fp_target, TailSide::Left))
        .collect()
}

pub fn tune_max_edits(
    searcher: &mut Searcher<Iupac>,
    flanks: &Vec<FlankGroup>,
    n: usize,
    fp_target: f64,
) -> i32 {
    let mut sp = Spinner::new(Spinners::OrangeBluePulse, "Tuning max edits".into());

    let half_query_len = (flanks[0].flank_seq.mask_queries[0].len() as f32 / 2.0).round() as usize;

    let mut min_edits = vec![];
    for i in 0..n {
        let random_seq = generate_random_sequence(half_query_len, half_query_len * 2);
        for flank in flanks {
            for flank_query in flank.flank_seq.mask_queries.iter() {
                let matches = searcher.search(&flank_query, &random_seq, half_query_len);
                let edit_dist = matches.iter().min_by_key(|m| m.cost).unwrap().cost;
                min_edits.push(edit_dist);
            }
        }
    }

    let min_edit = get_fp_threshold(min_edits, fp_target, TailSide::Left);
    sp.stop();

    min_edit
}

//     // let min_diff: Vec<u8> = (0..n).into_par_iter()
//     // .map(|_| {
//     //     let random_seq = generate_random_sequence(mask_query_slices[0].len() + 10);
//     //     let mut edits = mask_transposed_queries
//     //         .iter()
//     //         .flat_map(|tq| simd_search(tq, &random_seq))
//     //         .collect::<Vec<_>>();

//     //     edits.sort_unstable(); // Sort to get the top two
//     //     if edits.len() >= 2 {
//     //         edits[1].saturating_sub(edits[0]) // difference between top 2
//     //     } else {
//     //         0 // fallback in case there's only one result
//     //     }
//     // })
//     // .collect();

//     // let min_edit_diff = get_fp_threshold(min_diff, fp_target, TailSide::Right);
//     // sp.stop();

//     // min_edit_diff

//     // This loop can be parallelized
//     let min_edits: Vec<u8> = (0..n)
//         .into_par_iter()
//         .map(|_| {
//             let random_seq = generate_random_sequence(mask_query_slices[0].len() + 10);
//             let edits = mask_transposed_queries
//                 .iter()
//                 .flat_map(|tq| simd_search(tq, &random_seq))
//                 .collect::<Vec<_>>();
//             *edits.iter().min().unwrap()
//         })
//         .collect();

//     let min_edit = get_fp_threshold(min_edits, fp_target, TailSide::Left);
//     sp.stop();

//     min_edit
// }

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
