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
