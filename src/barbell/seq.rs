use pa_types::Pos;
use pa_bitpacking::search::*;
use crate::barbell::misc::*;
use plotly::{Plot, Scatter, Layout};
use crate::barbell::reader::{Query, QueryGroup};
use crate::barbell::strategy::*;

/// Returns the reverse complement of a DNA sequence
/// Based on https://www.bioinformatics.org/sms/iupac.html
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    println!("Input to reverse_complement: {}", String::from_utf8_lossy(seq));
    // Create a reversed vector first
    let reversed: Vec<u8> = seq.iter().rev().cloned().collect();
    println!("After reverse: {}", String::from_utf8_lossy(&reversed));
    
    // Then complement the reversed sequence
    let result = reversed.iter()
        .map(|&base| match base {
            b'A' => b'T', b'T' => b'A',
            b'C' => b'G', b'G' => b'C',
            b'R' => b'Y', b'Y' => b'R',
            b'S' => b'S', // S is self-complementary (G or C)
            b'W' => b'W', // W is self-complementary (A or T)
            b'K' => b'M', b'M' => b'K',
            b'B' => b'V', b'V' => b'B',
            b'D' => b'H', b'H' => b'D',
            b'N' => b'N',
            b'*' => b'*', // Add explicit handling for mask character
            _ => base,
        })
        .collect::<Vec<u8>>();
    println!("Output from reverse_complement: {}", String::from_utf8_lossy(&result));
    result
}



#[derive(Clone, Debug)]
pub enum MatchType {
    Flank, 
    Barcode,
}

#[derive(Clone, Debug)]
pub struct Match {
    pub some_id: String,
    pub match_type: MatchType,
    pub start: usize,
    pub end: usize,
    pub edits: i32,
    pub rel_dist_to_end: isize,
}


impl Match {
    pub fn new(some_id: String, match_type: MatchType, start: usize, end: usize, edits: i32, rel_dist_to_end: isize) -> Self {
        Self { some_id, match_type, start, end, edits, rel_dist_to_end }
    }
}




#[derive(Clone, Debug)]
pub struct Read<'a, S: Strategy> {
    seq: &'a [u8], // Reference to seq io buffer
    flank_locations: Vec<Match>,
    barcodes: Vec<BarcodeMatch>,
    strategy: S,
}

impl<'a, S: Strategy> Read<'a, S> {
    
    pub fn from_seq(seq: &'a [u8]) -> Self {
        Self { 
            seq, 
            flank_locations: Vec::new(),
            barcodes: Vec::new(),
            strategy: S::default(), // Note: Strategy trait needs to implement Default
        }
    }

    pub fn with_strategy(mut self, strategy: S) -> Self {
        self.strategy = strategy;
        self
    }

    pub fn annotate(&mut self, query_groups: &[QueryGroup]) -> Vec<Match> {  
        self.strategy.annotate(query_groups, self.seq)
    }

    pub fn final_assignment(&mut self, matches: &[Match]) -> Vec<Match> {
        self.strategy.final_assignment(matches)
    }

    
}


#[derive(Clone, Debug)]
pub struct BarcodeMatch {
    pub barcode_idx: usize,
    pub edits: i32,
    pub pos: usize,
}

impl BarcodeMatch {
    pub fn new(barcode_idx: usize, edits: i32, pos: usize) -> Self {
        Self { barcode_idx, edits, pos }
    }
}

#[derive(Clone, Debug)]
pub struct FlankMatch {
    pub start: usize,
    pub end: usize,
    pub edits: i32,
    pub barcodes: Option<Vec<BarcodeMatch>>, // Could be empty if no barcode passess filtering
}

impl FlankMatch {
    pub fn new(start: usize, end: usize, edits: i32) -> Self {
        Self { start, end, edits, barcodes: None }
    }
}


#[derive(Clone, Debug)]
pub struct FlankSeq {
    pub seq: Vec<u8>,
    mask_region: Option<(usize, usize)>,
}


impl FlankSeq {
    pub fn empty() -> Self {
        Self { 
            seq: Vec::new(),
            mask_region: None,
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

    pub fn mask_covered(&self, traceback_path: &[pa_types::Pos], min_frac: f32) -> (f32, bool) {

        // println!("Mask region: {:?}", self.mask_region);
        // If there is no mask, there is just one barcode, and hence no competition
        // we can just say it matched 
        let (mask_start, mask_end) = match self.mask_region {
            Some(region) => region,
            None => return (100.0, true),
        };

        // We use the traceback to check if at least min_frac of the mask region is covered
        let positions_in_mask = traceback_path
            .iter()
            .filter(|Pos(_, q_pos)| *q_pos as usize >= mask_start && (*q_pos as usize) < mask_end)
            .count();

        let mask_length = mask_end - mask_start;
        let frac = (positions_in_mask as f32) / (mask_length as f32);
        (frac, frac >= min_frac)
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
        }
    }
    
    /// Creates a new FlankSeq with optional prefix and suffix sequences,
    /// and a mask region of specified size between them
    pub fn with_flanks(prefix: Option<&[u8]>, suffix: Option<&[u8]>, mask_size: usize) -> Self {
        let mut seq = Vec::new();
        let mut mask_start = 0;
        
        // Add prefix if present
        if let Some(p) = prefix {
            seq.extend_from_slice(p);
            mask_start = seq.len();
        }
        
        // Add mask region
        seq.extend(vec![b'*'; mask_size]);
        let mask_end = seq.len();
        
        // Add suffix if present
        if let Some(s) = suffix {
            seq.extend_from_slice(s);
        }

        let mask_region = if mask_size > 0 {
            Some((mask_start, mask_end))
        } else {
            None
        };

        Self {
            seq: seq.clone(),
            mask_region,
        }
    }
}




mod test {  

    use super::*;
    use pa_bitpacking::search;
    

    fn first_lowest_index(edits: &[i32]) -> usize {
        let min_edits = edits.iter().min().unwrap();
        edits.iter()
            .position(|&e| e == *min_edits)
            .unwrap()
    }

    const LONG_PREFIX: &str = "ACTAGCTACGATCATCTACGATCAGCTACGATCAGCTACGATCAGCTAC";
    const LONG_SUFFIX: &str = "TACGATCAGCTACGATCAGCTACGATCAGCTACGATCAGCTACGATCAGCTAC";
    const RANDOM_SEQ: &str = "ATCTGATCGACTTCTCTGCTGCGCGGCGCGCGAGCTTACTTCTGCGCGATCTGGTACGGACGATTCTCAGACTACGGCAGTCGACT";

    #[test]
    fn test_mask_covered() {
        let t =     b"GGGGAAACCCC";    
        let f = FlankSeq::with_flanks(Some(b"GGGG"), Some(b"CCCC"), 3);
        let s = search(&f.seq, t, 0.4);
        let first_lowest_index = first_lowest_index(&s.out);
        let traceback = s.trace(first_lowest_index);
        let mask_covered = f.mask_covered(&traceback.1[..], 1.0);
        assert!(mask_covered.1);
    }

    #[test]
    fn test_flank_reverse_complement() {
        let f = FlankSeq::with_flanks(Some(b"GGGG"), Some(b"AAAA"), 3);
        let rc = f.reverse_complement();
        println!("Original: {}", f.as_string());
        println!("Reverse complement: {}", rc.as_string());
        assert_eq!(rc.as_string(), "TTTT***CCCC");
    }

    #[test]
    fn test_mask_suffix_not_covered() {
        /*
            text      TTTTTTTTTTTTTTTPPPPPPPPPP
            pattern                  PPPPPPPPPP[MMMMMMMMMM]
         */
        let t =  [RANDOM_SEQ.as_bytes(), LONG_PREFIX.as_bytes()].concat(); 
        let f = FlankSeq::with_flanks(Some(LONG_PREFIX.as_bytes()), None, 24);
        
        let s = search(&f.seq, &t, 0.4);
        let first_lowest_index = first_lowest_index(&s.out);
        let traceback = s.trace(first_lowest_index);
        let mask_covered = f.mask_covered(&traceback.1[..], 1.0);
        assert!(!mask_covered.1);
    }

    #[test]
    fn test_mask_prefix_not_covered() {
        /* 
            text                   PPPPPPPPPPPPPPTTTTTTTTTTTTTTTTTTT
            pattern  [MMMMMMMMMMMM]PPPPPPPPPPPPPP
         */
        let t =  [LONG_PREFIX.as_bytes(), RANDOM_SEQ.as_bytes()].concat();
        let f = FlankSeq::with_flanks(None, Some(LONG_PREFIX.as_bytes()), 24);
        let s = search(&f.seq, &t, 0.4);
        let first_lowest_index = first_lowest_index(&s.out);
        let traceback = s.trace(first_lowest_index);
        let mask_covered = f.mask_covered(&traceback.1[..], 1.0);
        assert!(!mask_covered.1);
    }

    #[test]
    fn test_mask_partially_covered() {
         /* 
            text                   PPPPPPPPPPPPPPTTTTTTTTTTTTTTTTTTT
            pattern         [MMMMMMMMMMMM]PPPPPPP half covered mask
         */
        let t =  [LONG_PREFIX.as_bytes(), RANDOM_SEQ.as_bytes()].concat();

        let half_long_prefix = &LONG_PREFIX.as_bytes()[12..LONG_PREFIX.len()]; // Half mask
        let f = FlankSeq::with_flanks(None, Some(half_long_prefix), 24);

        let s = search(&f.seq, &t, 0.4);
        let first_lowest_index = first_lowest_index(&s.out);
        let traceback = s.trace(first_lowest_index);
        
        // Check if all is covered (returns 0.5 is covered, and false for full coverage)
        let mask_covered = f.mask_covered(&traceback.1[..], 1.0);
        assert!(!mask_covered.1);
        assert!(mask_covered.0 == 0.5);

        // Request with 0.5 coverage
        let mask_covered = f.mask_covered(&traceback.1[..], 0.5);
        assert!(mask_covered.1);
        assert!(mask_covered.0 == 0.5);


    }


}

