/*
   Basic operations on sequences
*/

pub fn longest_common_prefix(seqs: &[&[u8]]) -> Option<Vec<u8>> {
    // Assuming not empty
    let first = &seqs[0];
    let mut prefix_length = first.len();

    for s in seqs.iter().skip(1) {
        let common_length = first.iter()
            .zip(s.iter())
            .take_while(|(a, b)| a == b)
            .count();
        
        prefix_length = prefix_length.min(common_length);

        if prefix_length == 0 {
            return None;
        }
    }

    Some(first[..prefix_length].to_vec())
}


pub fn longest_common_suffix(seqs: &[&[u8]]) -> Option<Vec<u8>> {
    // Or rev all seq and call prefix
    let first = &seqs[0];
    let mut suffix_length = first.len();

    for s in seqs.iter().skip(1) {
        let common_length = first.iter().rev()
            .zip(s.iter().rev())
            .take_while(|(a, b)| a == b)
            .count();
        
        suffix_length = suffix_length.min(common_length);

        if suffix_length == 0 {
            return None;
        }
    }

    Some(first[first.len()-suffix_length..].to_vec())
}


/*
    Operations on scores
*/


// Maybe add this later when flank scores appear to be relevant still



// pub fn remove_contained(ranges: Vec<(usize, usize, i32)>) -> Vec<(usize, usize, i32)> {
//     let mut sorted_ranges = ranges;
//     sorted_ranges.sort_by_key(|&(start, _, _)| start);

//     let mut result = Vec::new();
    
// }



// todo! more efficient way to do this, why was this such a hassle?
pub fn find_prominent_minima(points: &[(usize, i32)], min_prominence: f32) -> Vec<usize> {
    let mut minima = Vec::new();
    let n = points.len();
    let window_size = 10;

    if n < 2 {
        return minima;
    }

    for i in 0..n {
        // Check if current point is minimum in its window
        let mut is_minimum = true;
        
        // Look at left window
        let left_start = i.saturating_sub(window_size);
        for j in left_start..i {
            if points[j].0 + window_size >= points[i].0 && points[j].1 <= points[i].1 {
                is_minimum = false;
                break;
            }
        }

        // Look at right window
        if is_minimum {
            let right_end = (i + window_size + 1).min(n);
            for j in (i + 1)..right_end {
                if points[j].0 <= points[i].0 + window_size && points[j].1 < points[i].1 {
                    is_minimum = false;
                    break;
                }
            }
        }

        if !is_minimum {
            continue;
        }

        // Compute prominence using the window
        let mut left_max = points[i].1;
        let mut j = i;
        while j > 0 && points[j - 1].0 >= points[i].0 - window_size {
            j -= 1;
            left_max = left_max.max(points[j].1);
        }
        let left_prominence = (left_max - points[i].1) as f32;

        let mut right_max = points[i].1;
        j = i;
        while j < n - 1 && points[j + 1].0 <= points[i].0 + window_size {
            j += 1;
            right_max = right_max.max(points[j].1);
        }
        let right_prominence = (right_max - points[i].1) as f32;

        // Use the smaller prominence value
        let prominence = left_prominence.min(right_prominence);

        if prominence >= min_prominence {
            minima.push(points[i].0);
        }
    }

    minima
}


pub fn rel_dist_to_end(pos: isize, read_len: usize) -> isize {
    // If already negative, for a match starting before the read we can return 1
    if pos < 0 {
        return 1;
    }

    if pos <= (read_len / 2) as isize {
        if pos == 0 {
            1  // Left end
        } else {
            pos  // Distance from left end (already isize)
        }
    } else if pos == read_len as isize {
        -1  // Right end
    } else {
        -(read_len as isize - pos)  // Distance from right end
    }
}



#[cfg(test)]
mod test {

    use super::*;


    // Seqeunce tests 

    #[test]
    fn test_prefix_simple() {
        let seq1 = b"ATCGATCG".as_slice();
        let seq2 = b"ATCTTATCG".as_slice();
        let seq3 = b"ATCGGGG".as_slice();
        let seqs = vec![seq1, seq2, seq3];

        let prefix = longest_common_prefix(&seqs);
        assert_eq!(prefix, Some(b"ATC".to_vec()));
    }

    #[test]
    fn test_prefix_shorter_seq() {
        let seq1 = b"ATCGATCG".as_slice();
        let seq2 = b"ATCTTATCG".as_slice();
        let seq3 = b"AT".as_slice();
        let seqs = vec![seq1, seq2, seq3];

        let prefix = longest_common_prefix(&seqs);
        assert_eq!(prefix, Some(b"AT".to_vec()));
    }

    #[test]
    fn test_suffix_simple() {
        let seq1 = b"ATCGATCG".as_slice();
        let seq2 = b"ATCTTATCG".as_slice();
        let seq3 = b"ATCGGGG".as_slice();
        let seqs = vec![seq1, seq2, seq3];

        let suffix = longest_common_suffix(&seqs);
        assert_eq!(suffix, Some(b"G".to_vec()));
    }

    #[test]
    fn test_suffix_absent() {
        let seq1 = b"ATCGATCG".as_slice();
        let seq2 = b"ATCTTATCG".as_slice();
        let seq3 = b"ATCGGGA".as_slice();
        let seqs = vec![seq1, seq2, seq3];

        let suffix = longest_common_suffix(&seqs);
        assert_eq!(suffix, None);
    }

     // score tests 

   
    // #[test]
    // fn test_find_prominent_minima_left_edge() {
    //     let scores: Vec<(usize, i32)> = vec![(0, 1), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10)];
    //     let minima = find_prominent_minima(&scores, 1.0);
    //     assert_eq!(minima, vec![0]);
    // }

    // #[test]
    // fn test_find_prominent_minima_right_edge() {
    //     let scores: Vec<(usize, i32)> = vec![(0, 10), (1, 9), (2, 8), (3, 7), (4, 6), (5, 5), (6, 4), (7, 3), (8, 3), (9, 1)];
    //     let minima = find_prominent_minima(&scores, 1.0);
    //     assert_eq!(minima, vec![9]);
    // }

    // #[test]
    // fn test_find_prominent_minima_middle_vec() {
    //     let scores: Vec<(usize, i32)> = vec![(0, 10), (1, 9), (2, 8), (3, 7), (4, 6), (5, 7), (6, 8), (7, 9), (8, 10)];
    //     let minima = find_prominent_minima(&scores, 1.0);
    //     assert_eq!(minima, vec![4]);
    // }

    // #[test]
    // fn test_find_prominent_double_minima() {
    //     let scores: Vec<(usize, i32)> = vec![(0, 10), (1, 9), (2, 7), (3, 3), (4, 5), (5, 8), (6, 6), (7, 4), (8, 3)];
    //     let minima = find_prominent_minima(&scores, 1.0);
    //     assert_eq!(minima, vec![3, 8]);
    // }

    // #[test]
    // fn test_prominent_minima_no_minima() {
    //     let scores: Vec<(usize, i32)> = vec![(0, 10), (1,10), (2, 11), (3, 12)];
    //     let minima = find_prominent_minima(&scores, 1.0);
    //     assert_eq!(minima.len(), 1);
    // }
    
    // #[test]
    // fn test_prominent_minima_not_adjacent() {
    //     let scores: Vec<(usize, i32)> = vec![(0, 10), (1,11), (2, 8), (3, 4), (100, 8), (101, 6), (102, 4), (103, 8)];
    //     let minima = find_prominent_minima(&scores, 2.0);
    //     assert_eq!(minima, vec![3, 102]);
    // }
    
    // #[test]
    // fn test_prominent_minima_respecting_gap() {
    //     // Especially important for this use case to not have flase positive flanks
    //     let scores: Vec<(usize, i32)> = vec![(0, 10), (4, 8), (100, 6), (200, 8), (300, 10)];
    //     let minima = find_prominent_minima(&scores, 1.0);
    //     assert_eq!(minima.len(), 0);
    // }

    // #[test]
    // fn test_prominent_minima_plateau() {
    //     let scores: Vec<(usize, i32)> = vec![(0, 10), (1, 4), (2, 4), (3, 4), (4, 10)];
    //     let minima = find_prominent_minima(&scores, 1.0);
    //     println!("Platea Minima: {:?}", minima);

    // }

    // #[test]
    // fn test_remove_overlap_simple() {
    //     let ranges = vec![(1, 4, 10), (2, 6, 10), (8, 10, 10)];
    //     let result = remove_overlap(ranges);
    //     // Same edit distance for (1,4) and (2,6), we keep the shorter
    //     assert_eq!(result, vec![(1, 4, 10), (8, 10, 10)]);
    // }

    // #[test]
    // fn test_remove_overlap_complex() {
    //     // Not a realistic test example
    //     let ranges = vec![(1, 5, 10), (2, 3, 10), (4, 8, 10), (7, 10, 10), (12, 15, 10)];
    //     //                            lengths:    4       2       4       3       3
    //     //                                  1---5                  (4)
    //     //                                    23                   (1)
    //     //                                      4---8              (4)
    //     //                                         7--10           (3)
    //     //                                                12--15   (3)
    //     // Solution                          1--5         12--15 
    //     let result = remove_overlap(ranges);
    //     assert_eq!(result, vec![(2, 3, 10), (12, 15, 10)]);
    // }

    // #[test]
    // fn test_remove_overlap_diff_edits() {
    //     // Not a realistic test example
    //     let ranges = vec![(1, 5, 10), (2, 3, 100), (4, 8, 10), (7, 10, 10), (12, 15, 10)];
    //     //                            lengths:    4       2       4       3       3
    //     //                                  1---5                  (4)
    //     //                                    23                   (1) > 100 edits now
    //     //                                      4---8              (4)
    //     //                                         7--10           (3) < pick this, 2nd shortest
    //     //                                                12--15   (3)
    //     // Solution                          1--5         12--15 
    //     let result = remove_overlap(ranges);
    //     assert_eq!(result, vec![(7, 10, 10), (12, 15, 10)]);
    // }


}