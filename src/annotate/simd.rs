#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::annotate::mutations::MAX_LEN;
use crate::annotate::mutations::ErrorRates;

#[target_feature(enable = "avx2")]
pub unsafe fn simd_metric(
    queries: &[&[u8]],
    target: &[u8],
    query_length: usize,
    target_length: usize,
    error_rates: &ErrorRates
) -> Vec<f64> {
    assert!(query_length <= MAX_LEN, "Query length exceeds MAX_LEN");
    assert!(queries.len() <= 16, "Number of queries exceeds 16");
    
    // Find smallest absolute value to normalize by
    let min_abs = error_rates.match_probability.abs().min(error_rates.substitution.abs()).min(error_rates.insertion.abs()).min(error_rates.deletion.abs());
    
    // Calculate mathematically equivalent integer scores
    // We'll multiply by a factor and round to get clean integers
    // Should we just precompute these? kinda waste to do it over and over for static input 
    let factor = 1.0; // Adjust as needed for precision vs. magnitude
    let match_score_scaled = (error_rates.match_probability / min_abs * factor).round() as i16;
    let sub_score_scaled = (error_rates.substitution / min_abs * factor).round() as i16;
    let ins_score_scaled = (error_rates.insertion / min_abs * factor).round() as i16; 
    let del_score_scaled = (error_rates.deletion / min_abs * factor).round() as i16;
   
    // Store original scale factor for converting back at the end
    let scale_ratio = min_abs / factor;

    // Initialize two DP rows (curr, prev)
    let mut dp = [[_mm256_set1_epi16(0); MAX_LEN + 1]; 2];

    // Main DP loop
    for i in 1..=query_length {
        let prev = (i - 1) & 1;
        let curr = i & 1;

        // Initialize first column with gap penalties
        dp[curr][0] = _mm256_set1_epi16(i as i16 * del_score_scaled);

        for j in 1..=target_length {
            // Get target base at j-1
            let t = target[j - 1];

            let mut sub_vec = [sub_score_scaled; 16]; // Default to mismatch
            for q in 0..queries.len() {
                let q_base = queries[q][i - 1];
                // Branchless version using arithmetic
                let is_match = (q_base == t) as i16;
                sub_vec[q] = is_match * match_score_scaled + (1 - is_match) * sub_score_scaled;
            }

            let sub_simd = _mm256_loadu_si256(sub_vec.as_ptr() as *const __m256i);

            // Calculate three possible moves
            let diag = _mm256_add_epi16(dp[prev][j - 1], sub_simd);
            let up   = _mm256_add_epi16(dp[prev][j], _mm256_set1_epi16(del_score_scaled));
            let left = _mm256_add_epi16(dp[curr][j - 1], _mm256_set1_epi16(ins_score_scaled));

            // Take minimum of the three
            let min1 = _mm256_max_epi16(diag, up);
            dp[curr][j] = _mm256_max_epi16(min1, left);
        }
    }

    // Extract final values - for semi-global, this is the minimum in the last row
    let final_row = query_length & 1;
    let mut max_log_probs = vec![f64::MIN; queries.len()];

    for j in 0..=target_length {
        let scores = dp[final_row][j];
        let mut temp = [i16::MIN; 16];
        _mm256_storeu_si256(temp.as_mut_ptr() as *mut __m256i, scores);

        for q in 0..queries.len() {
            max_log_probs[q] = f64::max(max_log_probs[q], temp[q] as f64);
        }
    }

    // Fix the map function for Vec
    let mut max_scores_log = vec![f64::MIN; queries.len()];
    for (i, score) in max_log_probs.iter().enumerate() {
        max_scores_log[i] = score * scale_ratio;
    }
    max_scores_log
}


mod test {

    use super::*;

    #[test]
    fn test_simd_metric() {
        let queries = [b"ATTA", b"GGGG"];
        let target = b"AAAAAATTAAAAAA";
        let query_length = 4;
        let target_length = target.len();

        /*
            all matches:       math.log(0.93) * 4 -> -0.2902827713393415
            all gaps:          math.log(0.05) * 4 -> -11.982929094215963
         */

        let error_rates = ErrorRates::default();
        
        // Create a slice of references to pass to the function
        let query_refs: Vec<&[u8]> = queries.iter().map(|q| q.as_slice()).collect();
        
        let result = unsafe { 
            simd_metric(&query_refs, target, query_length, target_length, &error_rates) 
        };
        
        // Expected values from the calculation above
        let expected_matches = -0.2902827713393415;
        let expected_gaps = -11.982929094215963;
        
        // Check if results are within 0.1 of expected values (more lenient tolerance)
        assert!((result[0] - expected_matches).abs() < 0.1, 
                "Expected approx {}, got {}", expected_matches, result[0]);
        assert!((result[1] - expected_gaps).abs() < 0.1, 
                "Expected approx {}, got {}", expected_gaps, result[1]);
    }

    #[test]
    #[should_panic(expected = "Number of queries exceeds 16")]
    fn test_error_more_sequences() {
        let queries = [b"ATTA"].repeat(17);
        let target = b"AAAAAATTAAAAAA";
        let query_length = 4;
        let target_length = target.len();
        let error_rates = ErrorRates::default();
        let query_refs: Vec<&[u8]> = queries.iter().map(|q| q.as_slice()).collect();
        let _ = unsafe { 
            simd_metric(&query_refs, target, query_length, target_length, &error_rates) 
        };
    }

    #[test]
    #[should_panic(expected = "Query length exceeds MAX_LEN")]
    fn test_error_query_too_long() {
        let queries = [b"A"].repeat(51);
        let target = b"AAAAAATTAAAAAA";
        let query_length = 51;
        let target_length = target.len();
        let error_rates = ErrorRates::default();
        let query_refs: Vec<&[u8]> = queries.iter().map(|q| q.as_slice()).collect();
        let _ = unsafe { 
            simd_metric(&query_refs, target, query_length, target_length, &error_rates) 
        };
    }
}