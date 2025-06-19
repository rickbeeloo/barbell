use colored::Colorize;
use std::cmp::Ordering;

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

pub fn get_fp_threshold<T: Number>(scores: Vec<T>, fp_target: f32, tail_side: TailSide) -> T {
    let mut sorted_scores = scores.clone();

    // Sort based on tail side
    if tail_side == TailSide::Left {
        // Sort from lowest to highest
        sorted_scores.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    } else {
        // Sort from highest to lowest
        sorted_scores.sort_by(|a, b| b.partial_cmp(a).unwrap_or(Ordering::Equal));
    }
    println!("# scores: {:?}", sorted_scores.len());

    let target_index = (fp_target * sorted_scores.len() as f32) as usize;

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
