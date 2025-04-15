use crate::types::*;
use crate::annotate::flank::*;
use rand::Rng;
use spinners::{Spinner, Spinners};
use colored::Colorize;
use crate::annotate::mutations::*;


pub fn generate_random_sequence(length: usize) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let bases = b"ACGT";
    (0..length).map(|_| bases[rng.gen_range(0..4)]).collect()
}

pub fn tune_log_prob(flanks: &Vec<FlankGroup>, n: usize, fp_target: f64, error_rates: &ErrorRates) -> f64 {
    let mut sp = Spinner::new(Spinners::OrangeBluePulse, "Tuning log probability".into());

    let mut log_probs = Vec::with_capacity(n);

    // Number of mask queries per flank
    let total_comparisons = flanks.iter()
        .map(|f| f.flank_seq.mask_queries.len())
        .sum::<usize>();

    // Per comparison 
    let per_comp = n / total_comparisons;
    
    for flank_group in flanks {
        for mask_query in flank_group.flank_seq.mask_queries.iter() {
            for _ in 0..per_comp {
                let random_seq = generate_random_sequence(mask_query.len());
                let log_prob = metric(mask_query, &random_seq, error_rates);
                log_probs.push(log_prob);
            }
        }
    }
    
    // Sort the probabilities in descending order (best scores first)
    log_probs.sort_by(|a, b| b.partial_cmp(a).unwrap());

    sp.stop();

    // Create a combined list of all rates in proper order
    let mut all_rates = vec![0.0, 0.000001, 0.00001, 0.0001, 0.001, 0.01];
    for i in 2..=5 {
        all_rates.push(i as f64 / 100.0);
    }
    
    // Find closest rate that doesn't exceed target
    let closest_rate_idx = all_rates.iter()
        .enumerate()
        .filter(|(_, r)| **r <= fp_target)
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(idx, _)| idx)
        .unwrap_or(0);
    
    let closest_rate = all_rates[closest_rate_idx];
    
    // Get the cutoff value
    let cutoff_index = (closest_rate * log_probs.len() as f64) as usize;
    let cutoff_index = if closest_rate == 0.0 { 0 } else { cutoff_index.min(log_probs.len() - 1) };
    let target_cutoff = log_probs[cutoff_index];
    
    // Now print the table
    println!("\n\n{}", "Score Cutoffs for Different False Positive Rates".bold().underline());
    println!("{:<12} | {:<15} | {:<15}", 
        "FP Rate (%)".blue(), 
        "Score Cutoff".blue(), 
        "1 in X chance".blue());
    println!("{}", "-".repeat(48).dimmed());
    
    // Display all rates in order
    for (idx, &fp_rate) in all_rates.iter().enumerate() {
        let cutoff_index = (fp_rate * log_probs.len() as f64) as usize;
        
        // Handle edge cases
        let cutoff_index = if fp_rate == 0.0 { 
            0  // For 0%, use the very first (best) score
        } else {
            cutoff_index.min(log_probs.len() - 1)
        };
        
        let score_cutoff = log_probs[cutoff_index];
        
        // For 0%, the "1 in X" value is infinite, so handle specially
        let one_in_x = if fp_rate == 0.0 {
            "∞".to_string()  // Infinity symbol
        } else {
            format!("{:.0}", 1.0 / fp_rate)
        };
        
        // Format the rate percentage correctly based on magnitude
        let rate_str = if fp_rate < 0.01 {
            format!("{:.4}%", fp_rate * 100.0)
        } else {
            format!("{:.1}%", fp_rate * 100.0)
        };
        
        // ONLY highlight the row that exactly matches our selected cutoff
        let line_format = if idx == closest_rate_idx {
            format!("{:<12} | {:<15.6} | {:<15}", 
                rate_str.green().bold(), 
                format!("{:.6}", score_cutoff).green().bold(), 
                one_in_x.green().bold())
        } else {
            format!("{:<12} | {:<15.6} | {:<15}", 
                rate_str, 
                format!("{:.6}", score_cutoff), 
                one_in_x)
        };
        
        println!("{}", line_format);
    }
    
    // Calculate min and max for reporting
    let min_log_prob = log_probs.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max_log_prob = log_probs.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    
    println!("\n{}", "Overall Statistics".bold().underline());
    println!("  • Min log prob: {}", format!("{:.6}", min_log_prob).dimmed());
    println!("  • Max log prob: {}", format!("{:.6}", max_log_prob).dimmed());
    println!("  • Median log prob: {}", format!("{:.6}", log_probs[log_probs.len() / 2]).dimmed());
    println!("  • {} {}", "Selected cutoff:".green().bold(), format!("{:.6}", target_cutoff).green().bold());
    
    // Return the selected cutoff
    target_cutoff
}

