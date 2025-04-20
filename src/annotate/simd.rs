#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::annotate::mutations::MAX_LEN;
use crate::annotate::mutations::ErrorRates;
use crate::annotate::mutations::ErrorRatesAffine;

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
    let mut dp = [[unsafe { _mm256_set1_epi16(0) }; MAX_LEN + 1]; 2];

    // Main DP loop
    for i in 1..=query_length {
        let prev = (i - 1) & 1;
        let curr = i & 1;

        // Initialize first column with gap penalties
        dp[curr][0] = unsafe { _mm256_set1_epi16(i as i16 * del_score_scaled) };

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

            let sub_simd = unsafe { _mm256_loadu_si256(sub_vec.as_ptr() as *const __m256i) };

            // Calculate three possible moves
            let diag = unsafe { _mm256_add_epi16(dp[prev][j - 1], sub_simd) };
            let up   = unsafe { _mm256_add_epi16(dp[prev][j], _mm256_set1_epi16(del_score_scaled)) };
            let left = unsafe { _mm256_add_epi16(dp[curr][j - 1], _mm256_set1_epi16(ins_score_scaled)) };

            // Take minimum of the three
            let min1 = unsafe { _mm256_max_epi16(diag, up) };
            dp[curr][j] = unsafe { _mm256_max_epi16(min1, left) }   ;
        }
    }

    // Extract final values - for semi-global, this is the minimum in the last row
    let final_row = query_length & 1;
    let mut max_log_probs = vec![f64::MIN; queries.len()];

    for j in 0..=target_length {
        let scores = dp[final_row][j];
        let mut temp = [i16::MIN; 16];
        unsafe { _mm256_storeu_si256(temp.as_mut_ptr() as *mut __m256i, scores) };

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





#[target_feature(enable = "avx2")]
pub unsafe fn simd_metric_affine(
    queries: &[&[u8]],
    target: &[u8],
    query_length: usize,
    target_length: usize,
    error_rates: &ErrorRatesAffine
) -> Vec<f64> {
    const MAX_LEN: usize = 60;
    assert!(query_length <= MAX_LEN, "Query length exceeds MAX_LEN");
    assert!(queries.len() <= 16, "Number of queries exceeds 16");

    // Find smallest absolute value to normalize by
    let min_abs = error_rates.match_probability.abs().min(error_rates.substitution.abs())
                  .min(error_rates.ins_open.abs()).min(error_rates.ins_extend.abs())
                  .min(error_rates.del_open.abs()).min(error_rates.del_extend.abs());

    // Calculate mathematically equivalent integer scores
    let factor = 1.0;
    let match_score_scaled = (error_rates.match_probability / min_abs * factor).round() as i16;
    let sub_score_scaled = (error_rates.substitution / min_abs * factor).round() as i16;
    let ins_open_score_scaled = (error_rates.ins_open / min_abs * factor).round() as i16;
    let ins_extend_score_scaled = (error_rates.ins_extend / min_abs * factor).round() as i16;
    let del_open_score_scaled = (error_rates.del_open / min_abs * factor).round() as i16;
    let del_extend_score_scaled = (error_rates.del_extend / min_abs * factor).round() as i16;

    // Store original scale factor for converting back at the end
    let scale_ratio = min_abs / factor;

    // Use an array-based approach instead of full matrices
    // Only need two rows at a time (previous and current)
    let mut dp_M = [unsafe { _mm256_set1_epi16(0) }; MAX_LEN + 1];
    let mut dp_I = [unsafe { _mm256_set1_epi16(i16::MIN) }; MAX_LEN + 1];
    let mut dp_D = [unsafe { _mm256_set1_epi16(i16::MIN) }; MAX_LEN + 1];
    
    let mut prev_dp_M = [unsafe { _mm256_set1_epi16(0) }; MAX_LEN + 1];
    let mut prev_dp_I = [unsafe { _mm256_set1_epi16(i16::MIN) }; MAX_LEN + 1];
    let mut prev_dp_D = [unsafe { _mm256_set1_epi16(i16::MIN) }; MAX_LEN + 1];

    // Initialize first column with gap penalties
    for j in 1..=target_length {
        dp_M[j] = unsafe { _mm256_set1_epi16(i16::MIN) };
        dp_I[j] = unsafe { _mm256_set1_epi16(ins_open_score_scaled + (j as i16 - 1) * ins_extend_score_scaled) };
        dp_D[j] = unsafe { _mm256_set1_epi16(i16::MIN) };
    }

    // Pre-calculate the match/mismatch scores for each query position
    // This avoids recalculating them in the inner loop
    let mut query_bases = [[0u8; MAX_LEN]; 16];
    for q in 0..queries.len() {
        for i in 0..query_length {
            query_bases[q][i] = queries[q][i];
        }
    }

    // Precompute SIMD constants
    let ins_open_simd = unsafe { _mm256_set1_epi16(ins_open_score_scaled) };
    let ins_extend_simd = unsafe { _mm256_set1_epi16(ins_extend_score_scaled) };
    let del_open_simd = unsafe { _mm256_set1_epi16(del_open_score_scaled) };
    let del_extend_simd = unsafe { _mm256_set1_epi16(del_extend_score_scaled) };

    // Temporary storage for match/mismatch scores
    let mut sub_vec = [0i16; 16];

    // Main DP loop
    for i in 1..=query_length {
        // Swap current and previous rows
        std::mem::swap(&mut dp_M, &mut prev_dp_M);
        std::mem::swap(&mut dp_I, &mut prev_dp_I);
        std::mem::swap(&mut dp_D, &mut prev_dp_D);
        
        // Initialize first column of current row
        dp_M[0] = unsafe { _mm256_set1_epi16(i16::MIN) };
        dp_I[0] = unsafe { _mm256_set1_epi16(i16::MIN) };
        dp_D[0] = unsafe { _mm256_set1_epi16(del_open_score_scaled + (i as i16 - 1) * del_extend_score_scaled) };

        for j in 1..=target_length {
            let t = target[j - 1];

            // Prepare match/mismatch scores for each query
            for q in 0..queries.len() {
                let q_base = query_bases[q][i - 1];
                sub_vec[q] = if q_base == t { match_score_scaled } else { sub_score_scaled };
            }
            let sub_simd = unsafe { _mm256_loadu_si256(sub_vec.as_ptr() as *const __m256i) };

            // Update M matrix (Match/mismatch)
            let diag_m = unsafe { _mm256_adds_epi16(prev_dp_M[j - 1], sub_simd) };
            let diag_i = unsafe { _mm256_adds_epi16(prev_dp_I[j - 1], sub_simd) };
            let diag_d = unsafe { _mm256_adds_epi16(prev_dp_D[j - 1], sub_simd) };
            dp_M[j] = unsafe { _mm256_max_epi16(diag_m, _mm256_max_epi16(diag_i, diag_d)) };

            // Update I matrix (Insertion - gap in target)
            let m_to_i = unsafe { _mm256_adds_epi16(dp_M[j - 1], ins_open_simd) };
            let i_to_i = unsafe { _mm256_adds_epi16(dp_I[j - 1], ins_extend_simd) };
            let d_to_i = unsafe { _mm256_adds_epi16(dp_D[j - 1], ins_open_simd) };
            dp_I[j] = unsafe { _mm256_max_epi16(m_to_i, _mm256_max_epi16(i_to_i, d_to_i)) };

            // Update D matrix (Deletion - gap in query)
            let m_to_d = unsafe { _mm256_adds_epi16(prev_dp_M[j], del_open_simd) };
            let i_to_d = unsafe { _mm256_adds_epi16(prev_dp_I[j], del_open_simd) };
            let d_to_d = unsafe { _mm256_adds_epi16(prev_dp_D[j], del_extend_simd) };
            dp_D[j] = unsafe { _mm256_max_epi16(m_to_d, _mm256_max_epi16(i_to_d, d_to_d)) };
        }
    }

    // Extract final values
    let mut max_scores = vec![i16::MIN; queries.len()];

    // Check all positions in the final row
    for j in 1..=target_length {
        // Extract scores from all three matrices
        let mut temp_m = [0i16; 16];
        let mut temp_i = [0i16; 16];
        let mut temp_d = [0i16; 16];
        
        unsafe {
            _mm256_storeu_si256(temp_m.as_mut_ptr() as *mut __m256i, dp_M[j]);
            _mm256_storeu_si256(temp_i.as_mut_ptr() as *mut __m256i, dp_I[j]);
            _mm256_storeu_si256(temp_d.as_mut_ptr() as *mut __m256i, dp_D[j]);
        }

        // Update max scores
        for q in 0..queries.len() {
            max_scores[q] = i16::max(
                max_scores[q], 
                i16::max(temp_m[q], i16::max(temp_i[q], temp_d[q]))
            );
        }
    }

    // Convert scaled integers back to log probabilities
    let mut result = vec![0.0; queries.len()];
    for (i, &score) in max_scores.iter().enumerate() {
        result[i] = score as f64 * scale_ratio;
    }

    result
}

// #[target_feature(enable = "avx2")]
// pub unsafe fn simd_metric_affine(
//     queries: &[&[u8]],
//     target: &[u8],
//     query_length: usize,
//     target_length: usize,
//     error_rates: &ErrorRatesAffine
// ) -> Vec<f64> {
//     assert!(query_length <= MAX_LEN, "Query length exceeds MAX_LEN");
//     assert!(queries.len() <= 16, "Number of queries exceeds 16");

//     // Find smallest absolute value to normalize by
//     let min_abs = error_rates.match_probability.abs().min(error_rates.substitution.abs())
//                   .min(error_rates.ins_open.abs()).min(error_rates.ins_extend.abs())
//                   .min(error_rates.del_open.abs()).min(error_rates.del_extend.abs());

//     // Calculate mathematically equivalent integer scores
//     let factor = 1.0;
//     let match_score_scaled = (error_rates.match_probability / min_abs * factor).round() as i16;
//     let sub_score_scaled = (error_rates.substitution / min_abs * factor).round() as i16;
//     let ins_open_score_scaled = (error_rates.ins_open / min_abs * factor).round() as i16;
//     let ins_extend_score_scaled = (error_rates.ins_extend / min_abs * factor).round() as i16;
//     let del_open_score_scaled = (error_rates.del_open / min_abs * factor).round() as i16;
//     let del_extend_score_scaled = (error_rates.del_extend / min_abs * factor).round() as i16;

//     // Store original scale factor for converting back at the end
//     let scale_ratio = min_abs / factor;

//     // Initialize DP rows for match, deletion, and insertion
//     let mut dp_M = [[unsafe { _mm256_set1_epi16(0) }; MAX_LEN + 1]; 2];
//     let mut dp_I = [[unsafe { _mm256_set1_epi16(i16::MIN) }; MAX_LEN + 1]; 2];
//     let mut dp_D = [[unsafe { _mm256_set1_epi16(i16::MIN) }; MAX_LEN + 1]; 2];

//     // Initialize first column with gap penalties
//     for i in 1..=query_length {
//         let curr = i & 1;
//         dp_M[curr][0] = unsafe { _mm256_set1_epi16(i16::MIN) };
//         dp_I[curr][0] = unsafe { _mm256_set1_epi16(i16::MIN) };
//         dp_D[curr][0] = unsafe { _mm256_set1_epi16(del_open_score_scaled + (i as i16 - 1) * del_extend_score_scaled) };
//     }

//     // Main DP loop
//     for i in 1..=query_length {
//         let prev = (i - 1) & 1;
//         let curr = i & 1;

//         for j in 1..=target_length {
//             let t = target[j - 1];

//             // Prepare match/mismatch scores for each query
//             let mut sub_vec = [sub_score_scaled; 16];
//             for q in 0..queries.len() {
//                 let q_base = queries[q][i - 1];
//                 let is_match = (q_base == t) as i16;
//                 sub_vec[q] = is_match * match_score_scaled + (1 - is_match) * sub_score_scaled;
//             }
//             let sub_simd = unsafe { _mm256_loadu_si256(sub_vec.as_ptr() as *const __m256i) };

//             // Update M matrix (Match/mismatch)
//             let diag_m = unsafe { _mm256_adds_epi16(dp_M[prev][j - 1], sub_simd) };
//             let diag_i = unsafe { _mm256_adds_epi16(dp_I[prev][j - 1], sub_simd) };
//             let diag_d = unsafe { _mm256_adds_epi16(dp_D[prev][j - 1], sub_simd) };
//             dp_M[curr][j] = unsafe { _mm256_max_epi16(diag_m, _mm256_max_epi16(diag_i, diag_d)) };

//             // Update I matrix (Insertion - gap in target)
//             let m_to_i = unsafe { _mm256_adds_epi16(dp_M[curr][j - 1], _mm256_set1_epi16(ins_open_score_scaled)) };
//             let i_to_i = unsafe { _mm256_adds_epi16(dp_I[curr][j - 1], _mm256_set1_epi16(ins_extend_score_scaled)) };
//             let d_to_i = unsafe { _mm256_adds_epi16(dp_D[curr][j - 1], _mm256_set1_epi16(ins_open_score_scaled)) };
//             dp_I[curr][j] = unsafe { _mm256_max_epi16(m_to_i, _mm256_max_epi16(i_to_i, d_to_i)) };

//             // Update D matrix (Deletion - gap in query)
//             let m_to_d = unsafe { _mm256_adds_epi16(dp_M[prev][j], _mm256_set1_epi16(del_open_score_scaled)) };
//             let i_to_d = unsafe { _mm256_adds_epi16(dp_I[prev][j], _mm256_set1_epi16(del_open_score_scaled)) };
//             let d_to_d = unsafe { _mm256_adds_epi16(dp_D[prev][j], _mm256_set1_epi16(del_extend_score_scaled)) };
//             dp_D[curr][j] = unsafe { _mm256_max_epi16(m_to_d, _mm256_max_epi16(i_to_d, d_to_d)) };
//         }
//     }

//     // Extract final values - for semi-global alignment, we only consider the final row
//     let final_row = query_length & 1;
//     let mut max_scores = vec![i16::MIN; queries.len()];

//     // Check all positions in the final row
//     for j in 1..=target_length {
//         // Check all three matrices at each position
//         for matrix_idx in 0..3 {
//             let scores = match matrix_idx {
//                 0 => dp_M[final_row][j],
//                 1 => dp_I[final_row][j],
//                 _ => dp_D[final_row][j],
//             };
            
//             let mut temp = [i16::MIN; 16];
//             unsafe { _mm256_storeu_si256(temp.as_mut_ptr() as *mut __m256i, scores) };

//             for q in 0..queries.len() {
//                 max_scores[q] = i16::max(max_scores[q], temp[q]);
//             }
//         }
//     }

//     // Convert scaled integers back to log probabilities
//     let mut result = vec![f64::MIN; queries.len()];
//     for (i, &score) in max_scores.iter().enumerate() {
//         result[i] = score as f64 * scale_ratio;
//     }

//     result
// }



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