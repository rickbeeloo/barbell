pub const MAX_LEN: usize = 50;

use crate::annotate::simd::simd_metric;

#[derive(Debug, Clone)]
pub struct ErrorRates {
    pub substitution: f64,
    pub insertion: f64,
    pub deletion: f64,
    pub match_probability: f64,
}

impl ErrorRates {
    pub fn new(substitution: f64, insertion: f64, deletion: f64) -> Self {
        if substitution <= 0.0 || insertion <= 0.0 || deletion <= 0.0 {
            panic!("Probabilities must be positive");
        }
        if substitution + insertion + deletion >= 1.0 {
            panic!("Error probabilities should sum up to <= 1.0"); 
        }
        let match_probability = 1.0 - (substitution + insertion + deletion);
        ErrorRates {
            substitution: substitution.ln(),
            insertion: insertion.ln(),
            deletion: deletion.ln(),
            match_probability: match_probability.ln(),
        }
    }

    pub fn default() -> Self {
        ErrorRates::new(0.01, 0.01, 0.05)
    }
}

// See https://www.sciencedirect.com/science/article/pii/S2667325824001924
pub fn metric(query: &[u8], target: &[u8], error_rates: &ErrorRates) -> f64 {
    let m = query.len();
    let n = target.len();
    assert!(m <= MAX_LEN && n <= MAX_LEN, "Sequence length exceeds MAX_LEN");

    // Reuse, as it's static anyway
    let mut dp = [[0.0; MAX_LEN + 1]; 2];

    // First row is initialized to 0 for semi-global alignment
    // This allows the query to start at any position in the target without penalty
    for j in 0..=n {
        dp[0][j] = 0.0;
    }
    
    // Basic DP, maximizing log probability, considering mutation rates
    for i in 1..=m {
        let prev_row = (i - 1) & 1;
        let curr_row = i & 1;
        dp[curr_row][0] = dp[prev_row][0] + error_rates.deletion;

        for j in 1..=n {
            let match_prob = if query[i-1] == target[j-1] {
                error_rates.match_probability
            } else {
                error_rates.substitution
            };

            let diag = dp[prev_row][j-1] + match_prob;
            let up = dp[prev_row][j] + error_rates.deletion;
            let left = dp[curr_row][j-1] + error_rates.insertion;
            // Unlike edits, we now maximize (not minimize)
            dp[curr_row][j] = f64::max(diag, f64::max(up, left));
        }
    }
    
    // For semi-global alignment, find the best score in the last row or column
    // This allows the alignment to end at any position
    let mut max_score = dp[m & 1][n];
    for j in 0..=n {
        max_score = f64::max(max_score, dp[m & 1][j]);
    }
    
    max_score
}

pub fn optimal_metric(
    queries: &[&[u8]],
    target: &[u8],
    error_rates: &ErrorRates
) -> Vec<f64> {
    // Early return for empty queries
    if queries.is_empty() {
        return Vec::new();
    }
    
    // Pre-allocate result vector with the correct size
    let mut results = Vec::with_capacity(queries.len());
    
    #[cfg(target_feature = "avx2")]
    {
        if is_x86_feature_detected!("avx2") {
            // Process in chunks of 16 queries to maximize SIMD usage
            for chunk in queries.chunks(16) {
                // For each chunk, we need to create a properly sized slice
                unsafe {
                    let chunk_result = simd_metric(
                        chunk,
                        target,
                        chunk[0].len(), // All queries should have the same length
                        target.len(),
                        error_rates
                    );
                    
                    // Add this chunk's results to our final results vector
                    results.extend_from_slice(&chunk_result);
                }
            }
            
            // Return the combined results
            return results;
        }
    }
    
    // Fallback to scalar implementation if SIMD is not available
    queries.iter()
        .map(|q| metric(q, target, error_rates))
        .collect()
}


pub fn approximate_edit_distance(metric: f64, seq_len: usize, error_rates: &ErrorRates) -> f64 {
    let perfect_match = seq_len as f64 * error_rates.match_probability;
    let score_difference = perfect_match - metric;
    let average_error_prob = ((error_rates.substitution.exp() + error_rates.insertion.exp() + error_rates.deletion.exp()) / 3.0).ln();
    score_difference / average_error_prob.abs()
}

pub fn approx_edit_metric(read: &[u8], barcode: &[u8], error: &ErrorRates) -> usize {
    let metric = metric(read, barcode, error);
    approximate_edit_distance(metric, barcode.len(), error) as usize // round
}

pub fn posterior_probability_sorted(log_probs: &[f64]) -> f64 {
    if log_probs.is_empty() {
        return 0.0;
    }

    let max = log_probs[0]; // since they are sorted, the first is the max
    let mut sum_exp = 0.0;

    for &logp in log_probs {
        sum_exp += (logp - max).exp();
    }

    1.0 / sum_exp
}

mod test {
    use super::*;

    #[test]
    fn test_metric_identical() {
        let barcode =    b"AAGAAAGTTGTCGGTGTCTTTGTG"; // NB01 barcode
        let read_ident = b"AAGAAAGTTGTCGGTGTCTTTGTG"; 
        let error = ErrorRates::new(0.01, 0.01, 0.05);
        let metric_ident = metric(barcode, read_ident, &error);
        // Verify log likelihood, that is (1.0-0.07).ln() * 24
        let expected_log_likelihood = (1.0_f64 - 0.07).ln() * 24.0;
        assert!((metric_ident - expected_log_likelihood).abs() < 1e-10); //floating point precision
    }

    #[test]
    fn test_metric_mutated() {
        let barcode =    b"AAGAAAGTTGTCGGTGTCTTTGTG"; // NB01 barcode
        let read_mut =   b"AAGAAAGTTGGGTGTCTTTGTG";  // two chars deleted
        let error = ErrorRates::new(0.01, 0.01, 0.05);
        let metric_mut = metric(barcode, read_mut, &error);
        let matching_log_likelihood = (1.0_f64 - 0.07).ln() * 22.0;
        let deletion_log_likelihood = (0.05_f64).ln() * 2.0;
        let expected_log_likelihood = matching_log_likelihood + deletion_log_likelihood;
        assert!((metric_mut - expected_log_likelihood).abs() < 1e-10); //floating point precision
    }

    #[test]
    fn test_posterior_probability_equal() {
        let log_prob1 = (1.0_f64 - 0.07).ln() * 22.0;
        let log_prob2 = (1.0_f64 - 0.07).ln() * 22.0;
        let posterior = posterior_probability_sorted(&[log_prob1, log_prob2]);
        assert_eq!(posterior, 0.5);
    }

    #[test] 
    fn test_posterior_probability_different() {
        let log_prob1 = (1.0_f64 - 0.07).ln() * 24.0; // All match
        let log_prob2 = (1.0_f64 - 0.07).ln() * 22.0 + (0.05_f64).ln() * 2.0; // two deletions
        let posterior = posterior_probability_sorted(&[log_prob1, log_prob2]);
        assert!(posterior > 0.5);
    }

}

