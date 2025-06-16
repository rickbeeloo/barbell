use crate::types::*;
use pa_types::Pos;

/// Returns the reverse complement of a DNA sequence
/// Based on https://www.bioinformatics.org/sms/iupac.html
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    // Then complement the reversed sequence
    seq.iter()
        .rev()
        .map(|&base| match base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S', // S is self-complementary (G or C)
            b'W' => b'W', // W is self-complementary (A or T)
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            b'N' => b'N',
            b'*' => b'*', // Add explicit handling for mask character
            _ => base,
        })
        .collect::<Vec<u8>>()
}

/// Main external structs used in aligner
#[derive(Clone, Debug)]
pub struct FlankGroup {
    pub flank_seq: FlankSeq,
    pub orientation: Orientation,
    pub match_type: MatchType,
}

impl FlankGroup {
    pub fn new(
        sequences: Vec<Vec<u8>>,
        sequence_ids: Vec<String>,
        orientation: Orientation,
        match_type: MatchType,
    ) -> Self {
        let flank_seq = FlankSeq::init_from_sequences(&sequences, sequence_ids).unwrap();
        Self {
            flank_seq,
            orientation,
            match_type,
        }
    }

    pub fn reverse_complement(&self) -> Self {
        let flank_seq_rc = self.flank_seq.reverse_complement();
        let match_type_rc = match self.match_type {
            MatchType::Fbarcode => MatchType::Rbarcode,
            MatchType::Rbarcode => MatchType::Fbarcode,
            _ => MatchType::Flank,
        };
        let orientation_rc = match self.orientation {
            Orientation::Forward => Orientation::ReverseComplement,
            Orientation::ReverseComplement => Orientation::Forward,
        };
        Self {
            flank_seq: flank_seq_rc,
            orientation: orientation_rc,
            match_type: match_type_rc,
        }
    }
}

#[derive(Clone, Debug)]
pub struct FlankSeq {
    pub seq: Vec<u8>,
    pub mask_region: Option<(usize, usize)>,
    pub mask_queries: Vec<Vec<u8>>,
    pub mask_ids: Vec<String>,
}

impl FlankSeq {
    pub fn empty() -> Self {
        Self {
            seq: Vec::new(),
            mask_region: None,
            mask_queries: Vec::new(),
            mask_ids: Vec::new(),
        }
    }

    pub fn unmasked_len(&self) -> usize {
        let masked = self.mask_region.map_or(0, |(start, end)| end - start);
        self.seq.len() - masked
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    pub fn as_string(&self) -> String {
        String::from_utf8_lossy(&self.seq).to_string()
    }

    pub fn mask_len(&self) -> usize {
        self.mask_region.map_or(0, |(start, end)| end - start)
    }

    pub fn mask_covered(
        &self,
        traceback_path: &[pa_types::Pos],
        min_frac: f32,
    ) -> (f32, bool, Option<(usize, usize)>) {
        let (mask_start, mask_end) = match self.mask_region {
            Some(region) => region,
            None => return (1.0, true, None),
        };

        // Keep track of min and max r_pos that fall within the mask region
        let mut min_r_pos = usize::MAX;
        let mut max_r_pos = 0;

        // Count positions and track r_pos range
        let positions_in_mask = traceback_path
            .iter()
            .filter(|Pos(q_pos, r_pos)| {
                // Q and R is swapped in sassy
                if *q_pos as usize >= mask_start && (*q_pos as usize) <= mask_end {
                    // less than or euqal, as inclusive range
                    min_r_pos = min_r_pos.min(*r_pos as usize);
                    max_r_pos = max_r_pos.max(*r_pos as usize);
                    true
                } else {
                    false
                }
            })
            .count();

        let mask_length = mask_end - mask_start + 1; // +1 as inclusive range
        let frac = (positions_in_mask as f32) / (mask_length as f32);

        let r_pos_range = if positions_in_mask > 0 {
            Some((min_r_pos, max_r_pos))
        } else {
            None
        };

        (frac, frac >= min_frac, r_pos_range)
    }

    pub fn reverse_complement(&self) -> Self {
        let total_len = self.seq.len();

        let mask_region = self.mask_region.map(|(start, end)| {
            // Mirror the positions from the end of the sequence
            (total_len - end, total_len - start)
        });

        Self {
            seq: reverse_complement(&self.seq),
            mask_region,
            mask_queries: self
                .mask_queries
                .iter()
                .map(|bc| reverse_complement(bc))
                .collect(),
            mask_ids: self.mask_ids.iter().map(|id| id.clone()).collect(),
        }
    }

    pub fn set_mask_queries(&mut self, mask_queries: Vec<Vec<u8>>) {
        self.mask_queries = mask_queries;
    }

    pub fn set_mask_ids(&mut self, mask_ids: Vec<String>) {
        self.mask_ids = mask_ids;
    }

    pub fn init_from_sequences(sequences: &[Vec<u8>], sequence_ids: Vec<String>) -> Option<Self> {
        // We extract the shared prefixes and suffixes, input sequences should be of equal length
        let (prefix, suffix, full_len) = get_flanks(sequences, 1);

        //  println!("Prefix: {:?}", String::from_utf8_lossy(prefix.as_ref().unwrap()));
        //  println!("Suffix: {:?}", String::from_utf8_lossy(suffix.as_ref().unwrap()));

        // Initialize
        let flank_queries = match (prefix, suffix) {
            // [prefix]--MASK--[suffix]
            (Some(p), Some(s)) => {
                // It's possible there is only one sequence, then prefix  == suffix, and no mask left
                // in that case, we just use only prefix without mask
                // More exceptiosn probably?
                if p.len() == s.len() && p.len() == full_len {
                    Some(FlankSeq::with_flanks(Some(&p), None, 0))
                } else {
                    let mask_size = full_len - p.len() - s.len();
                    Some(FlankSeq::with_flanks(Some(&p), Some(&s), mask_size))
                }
            }

            // [prefix]--MASK
            (Some(p), None) => {
                let mask_size = full_len - p.len();
                Some(FlankSeq::with_flanks(Some(&p), None, mask_size))
            }

            // MASK--[suffix]
            (None, Some(s)) => {
                let mask_size = full_len - s.len();
                Some(FlankSeq::with_flanks(None, Some(&s), mask_size))
            }

            // No flanks found - should this be enforced?
            (None, None) => panic!("No shared flanks for your queries"),
        };

        // We now have the proper initialization of the flanking sequences
        // the queries for the mask will now be the slices from the original queries corresponding to the masked region
        // which could be different from the barcodes
        if let Some(mut flank_wo_mask_queries) = flank_queries {
            // Get the mask region from the flank
            let mask_region = flank_wo_mask_queries.mask_region.unwrap();

            // Get the mask queries from the original queries
            let mask_queries = sequences
                .iter()
                .map(
                    |s| s[mask_region.0..=mask_region.1].to_vec(), // note we always use inclusive range on the sequences
                )
                .collect();

            // Set the mask queries
            flank_wo_mask_queries.set_mask_queries(mask_queries);
            flank_wo_mask_queries.set_mask_ids(sequence_ids.to_vec());

            // Return the flank with mask queries
            return Some(flank_wo_mask_queries);
        }

        None
    }

    /// Creates a new FlankSeq with optional prefix and suffix sequences,
    /// and a mask region of specified size between them
    fn with_flanks(prefix: Option<&[u8]>, suffix: Option<&[u8]>, mask_size: usize) -> Self {
        let mut seq = Vec::new();
        let mut mask_start = 0;

        // Add prefix if present
        if let Some(p) = prefix {
            seq.extend_from_slice(p);
            mask_start = seq.len();
        }

        // Note mask size can exceed the length of the barcodes when barcodes would
        // share a common prefix or suffix.
        // Add mask region
        seq.extend(vec![b'N'; mask_size]);
        let mask_end = seq.len() - 1;

        // Add suffix if present
        if let Some(s) = suffix {
            seq.extend_from_slice(s);
        }

        let mask_region: Option<(usize, usize)> = if mask_size > 0 {
            Some((mask_start, mask_end))
        } else {
            None
        };

        Self {
            seq: seq.clone(),
            mask_region,
            mask_queries: Vec::new(), // Empty initialize we slice them out later
            mask_ids: Vec::new(),
        }
    }
}

pub fn longest_common_prefix(seqs: &[&[u8]]) -> Option<Vec<u8>> {
    let first = &seqs[0];
    let mut prefix_length = first.len();

    for s in seqs.iter().skip(1) {
        let common_length = first
            .iter()
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
    // Convert sequences to reversed vectors
    let reversed: Vec<Vec<u8>> = seqs
        .iter()
        .map(|s| s.iter().rev().copied().collect())
        .collect();
    let reversed_refs: Vec<&[u8]> = reversed.iter().map(|v| v.as_slice()).collect();
    longest_common_prefix(&reversed_refs).map(|mut prefix| {
        prefix.reverse();
        prefix
    })
}

pub fn get_flanks(
    seqs: &[Vec<u8>],
    min_length: usize,
) -> (Option<Vec<u8>>, Option<Vec<u8>>, usize) {
    // This does not make sense if the sequences are not equally long
    assert!(
        !seqs.iter().any(|s| s.len() != seqs[0].len()),
        "Sequences are not equally long"
    );

    let prefix = match longest_common_prefix(&seqs.iter().map(|s| s.as_slice()).collect::<Vec<_>>())
    {
        Some(p) if p.len() >= min_length => Some(p),
        _ => None,
    };

    let suffix = match longest_common_suffix(&seqs.iter().map(|s| s.as_slice()).collect::<Vec<_>>())
    {
        Some(s) if s.len() >= min_length => Some(s),
        _ => None,
    };

    (prefix, suffix, seqs[0].len())
}

#[cfg(test)]
mod test {

    use super::*;

    fn first_lowest_index(edits: &[i32]) -> usize {
        let min_edits = edits.iter().min().unwrap();
        edits.iter().position(|&e| e == *min_edits).unwrap()
    }

    const LONG_PREFIX: &str = "ACTAGCTACGATCATCTACGATCAGCTACGATCAGCTACGATCAGCTAC";
    const LONG_SUFFIX: &str = "TACGATCAGCTACGATCAGCTACGATCAGCTACGATCAGCTACGATCAGCTAC";
    const RANDOM_SEQ: &str =
        "ATCTGATCGACTTCTCTGCTGCGCGGCGCGCGAGCTTACTTCTGCGCGATCTGGTACGGACGATTCTCAGACTACGGCAGTCGACT";

    #[test]
    fn test_flank_reverse_complement() {
        let f = FlankSeq::init_from_sequences(
            &[b"TTTTAAAGGGG".to_vec(), b"TTTTCCCGGGG".to_vec()],
            vec!["seq1".to_string(), "seq2".to_string()],
        )
        .unwrap();
        let rc = f.reverse_complement();
        println!("Original: {}", f.as_string());
        println!("Reverse complement: {}", rc.as_string());
        assert_eq!(rc.as_string(), "CCCC***AAAA");
        assert_eq!(
            rc.mask_queries,
            vec![vec![b'T', b'T', b'T'], vec![b'G', b'G', b'G']]
        );
    }

    // #[test]
    // fn test_mask_covered_range() {
    //     //  We now where the mask, is, AAAA/TTTT
    //     let t = b"GGGGAAACCCC";
    //     let f = FlankSeq::init_from_sequences(
    //         &[b"GGGGAAACCCC".to_vec(), b"GGGGTTTCCCC".to_vec()],
    //         vec!["seq1".to_string(), "seq2".to_string()],
    //     )
    //     .unwrap();

    //     let s = search(&f.seq, t, 0.4);
    //     let first_lowest_index = first_lowest_index(&s.out);
    //     let traceback = s.trace(first_lowest_index);
    //     let (_, _, r_range) = f.mask_covered(&traceback.1[..], 1.0);

    //     // The mask should be positions 4,5,6 in the sequence (AAA)
    //     // And since this is an exact match, the r_range should match these positions
    //     let (start, end) = r_range.unwrap();
    //     assert_eq!(start, 4);
    //     assert_eq!(end, 6);

    //     // Verify these positions in the target sequence contain the mask
    //     assert_eq!(&t[start..=end], b"AAA");
    // }

    // #[test]
    // fn test_mask_partial_coverage() {
    //     //  We now where the mask, is, AAAA/TTTT
    //     let t = b"AAAAACCCC";
    //     let f = FlankSeq::init_from_sequences(
    //         &[b"GGGGAAAAAAAAACCCC".to_vec(), b"GGGGTTTTTTTTTCCCC".to_vec()],
    //         vec!["seq1".to_string(), "seq2".to_string()],
    //     )
    //     .unwrap();

    //     let s = search(&f.seq, t, 0.4);
    //     let first_lowest_index = first_lowest_index(&s.out);
    //     let traceback = s.trace(first_lowest_index);
    //     let (_, _, r_range) = f.mask_covered(&traceback.1[..], 1.0);
    //     let (start, end) = r_range.unwrap();
    //     assert_eq!(start, 0);
    //     assert_eq!(end, 4);
    //     assert_eq!(&t[start..=end], b"AAAAA");
    //     assert_eq!(end - start + 1, 5); // + 1 as inclusive range
    // }

    // #[test]
    // fn test_mask_covered_full_flank() {
    //     let t = b"GGGGAAACCCC";
    //     let f = FlankSeq::init_from_sequences(
    //         &[b"GGGGAAACCCC".to_vec(), b"GGGGTTTCCCC".to_vec()],
    //         vec!["seq1".to_string(), "seq2".to_string()],
    //     )
    //     .unwrap();

    //     let s = search(&f.seq, t, 0.4);
    //     let first_lowest_index = first_lowest_index(&s.out);
    //     let traceback = s.trace(first_lowest_index);
    //     let mask_covered = f.mask_covered(&traceback.1[..], 1.0);
    //     assert!(mask_covered.1);
    // }

    // #[test]
    // fn test_mask_covered_partial_flank() {
    //     let t = b"GGAAACCCC";
    //     let f = FlankSeq::init_from_sequences(
    //         &[b"GGGGAAACCCC".to_vec(), b"GGGGTTTCCCC".to_vec()],
    //         vec!["seq1".to_string(), "seq2".to_string()],
    //     )
    //     .unwrap();

    //     let s = search(&f.seq, t, 0.4);
    //     let first_lowest_index = first_lowest_index(&s.out);
    //     let traceback = s.trace(first_lowest_index);
    //     let mask_covered = f.mask_covered(&traceback.1[..], 1.0);
    //     assert!(mask_covered.1);
    // }

    // #[test]
    // fn test_mask_suffix_not_covered() {
    //     /*
    //        text      TTTTTTTTTTTTTTTPPPPPPPPPP
    //        pattern                  PPPPPPPPPP[MMMMMMMMMM]
    //     */
    //     let t = [RANDOM_SEQ.as_bytes(), LONG_PREFIX.as_bytes()].concat();

    //     let f = FlankSeq::init_from_sequences(
    //         &[
    //             {
    //                 let mut v = LONG_PREFIX.as_bytes().to_vec();
    //                 v.extend(vec![b'A'; 24]);
    //                 v
    //             },
    //             {
    //                 let mut v = LONG_PREFIX.as_bytes().to_vec();
    //                 v.extend(vec![b'T'; 24]);
    //                 v
    //             },
    //         ],
    //         vec!["seq1".to_string(), "seq2".to_string()],
    //     )
    //     .unwrap();

    //     let s = search(&f.seq, &t, 0.4);
    //     let first_lowest_index = first_lowest_index(&s.out);
    //     let traceback = s.trace(first_lowest_index);
    //     let mask_covered = f.mask_covered(&traceback.1[..], 1.0);
    //     assert!(!mask_covered.1);
    // }

    // #[test]
    // fn test_mask_prefix_not_covered() {
    //     /*
    //        text                   PPPPPPPPPPPPPPTTTTTTTTTTTTTTTTTTT
    //        pattern  [MMMMMMMMMMMM]PPPPPPPPPPPPPP
    //     */
    //     let t = [LONG_PREFIX.as_bytes(), RANDOM_SEQ.as_bytes()].concat();
    //     let f = FlankSeq::init_from_sequences(
    //         &[
    //             {
    //                 let mut v = vec![b'A'; 24];
    //                 v.extend(LONG_PREFIX.as_bytes());
    //                 v
    //             },
    //             {
    //                 let mut v = vec![b'T'; 24];
    //                 v.extend(LONG_PREFIX.as_bytes());
    //                 v
    //             },
    //         ],
    //         vec!["seq1".to_string(), "seq2".to_string()],
    //     )
    //     .unwrap();
    //     let s = search(&f.seq, &t, 0.4);
    //     let first_lowest_index = first_lowest_index(&s.out);
    //     let traceback = s.trace(first_lowest_index);
    //     let mask_covered = f.mask_covered(&traceback.1[..], 1.0);
    //     assert!(!mask_covered.1);
    // }

    // #[test]
    // fn test_mask_partially_covered() {
    //     /*
    //        text                  MMMMMMMMPPPPPPPPPPPPPPTTTTTTTTTTTTTTTTTTT
    //        pattern         [MMMMMMMMMMMM]PPPPPPP half covered mask
    //     */
    //     let t = [
    //         &vec![b'A'; 12],
    //         LONG_PREFIX.as_bytes(),
    //         RANDOM_SEQ.as_bytes(),
    //     ]
    //     .concat();

    //     let half_long_prefix = &LONG_PREFIX.as_bytes()[12..LONG_PREFIX.len()]; // Half mask

    //     let f = FlankSeq::init_from_sequences(
    //         &[
    //             {
    //                 let mut v = vec![b'A'; 24];
    //                 v.extend(LONG_PREFIX.as_bytes());
    //                 v
    //             },
    //             {
    //                 let mut v = vec![b'T'; 24];
    //                 v.extend(LONG_PREFIX.as_bytes());
    //                 v
    //             },
    //         ],
    //         vec!["seq1".to_string(), "seq2".to_string()],
    //     )
    //     .unwrap();

    //     let s = search(&f.seq, &t, 0.4);
    //     let first_lowest_index = first_lowest_index(&s.out);
    //     let traceback = s.trace(first_lowest_index);

    //     // Check if all is covered (returns 0.5 is covered, and false for full coverage)
    //     let mask_covered = f.mask_covered(&traceback.1[..], 1.0);
    //     assert!(!mask_covered.1);
    //     println!("Mask covered: {:?}", mask_covered);
    //     assert!(mask_covered.0 == 0.5);

    //     // Request with 0.5 coverage
    //     let mask_covered = f.mask_covered(&traceback.1[..], 0.5);
    //     assert!(mask_covered.1);
    //     assert!(mask_covered.0 == 0.5);
    // }

    #[test]
    fn test_mask_seq_case1() {
        /*
            text      GGGGAAACCCC
            pattern   GGGGTTTCCCC
                      ----   ---- > (AAA, TTT)
        */
        let f = FlankSeq::init_from_sequences(
            &[b"GGGGAAACCCC".to_vec(), b"GGGGTTTCCCC".to_vec()],
            vec!["seq1".to_string(), "seq2".to_string()],
        )
        .unwrap();

        let mask_queries = f.mask_queries;
        assert_eq!(
            mask_queries,
            vec![vec![b'A', b'A', b'A'], vec![b'T', b'T', b'T']]
        );
    }

    #[test]
    fn test_mask_seq_case2() {
        /*
            text      GGGGAACCCCC
            pattern   GGGGTTTCCCC
                      ----   ---- > (AAC, TTT)
        */
        let f = FlankSeq::init_from_sequences(
            &[b"GGGGAACCCCC".to_vec(), b"GGGGTTTCCCC".to_vec()],
            vec!["seq1".to_string(), "seq2".to_string()],
        )
        .unwrap();

        let mask_queries = f.mask_queries;
        assert_eq!(
            mask_queries,
            vec![vec![b'A', b'A', b'C'], vec![b'T', b'T', b'T']]
        );
    }

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
}
