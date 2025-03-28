use seq_io::fasta::{Reader,Record};
use crate::barbell::{misc::*, seq::*};
use colored::Colorize;
use crate::barbell::pattern_assign::*;


#[derive(Debug)]
pub struct Query {
    pub id: String,
    pub seq: Vec<u8>,
}

// Like this we don't later have to unwrap/match every time to see if we have something present
// we just generelize to searchig QueryGroup, seeded by flank, in the target
#[derive(Debug)]
pub struct QueryGroup {
    pub flank: Option<FlankSeq>,
    pub queries: Vec<Query>,
    pub query_ids: Vec<String>,
    pub match_type: MatchType, // As group identifier in output
    pub orientation: Orientation, // whether this is a reverse complement group, then we still now based on prefix char which were the same 
}

impl QueryGroup {
    pub fn new(flank: Option<FlankSeq>, queries: Vec<Query>, query_ids: Vec<String>, match_type: MatchType, orientation: Orientation) -> Self {
        Self { flank, queries, query_ids, match_type, orientation }
    }
}

fn read_fasta(path: &str) -> (Vec<Vec<u8>>, Vec<String>) {
    let mut reader = Reader::from_path(path).unwrap();
    let mut queries = Vec::new();
    let mut ids = Vec::new();
    
    while let Some(Ok(record)) = reader.records().next() {
        queries.push(record.seq().to_vec());
        ids.push(record.id().expect("Invalid seq id").to_string());
    }
    
    (queries, ids)
}

pub fn get_flanks(seqs: &[Vec<u8>], min_length: usize) -> (Option<Vec<u8>>, Option<Vec<u8>>, usize) {
    // This does not make sense if the sequences are not equally long
    assert!(
        !seqs.iter().any(|s| s.len() != seqs[0].len()),
        "Sequences are not equally long"
    );

    let prefix = match longest_common_prefix(&seqs.iter().map(|s| s.as_slice()).collect::<Vec<_>>()) {
        Some(p) if p.len() >= min_length => Some(p),
        _ => None,
    };
    
    let suffix = match longest_common_suffix(&seqs.iter().map(|s| s.as_slice()).collect::<Vec<_>>()) {
        Some(s) if s.len() >= min_length => Some(s),
        _ => None,
    };
    
    (prefix, suffix, seqs[0].len())
}

pub fn get_flank_seq(seqs: &[Vec<u8>], min_length: usize) -> Option<FlankSeq> {
    // Get shared prefix and suffix with minimum length requirement
    let (prefix, suffix, full_len) = get_flanks(seqs, min_length);

    match (prefix, suffix) {

        // [prefix]--MASK--[suffix]
        (Some(p), Some(s)) => {

            // It's possible there is only one sequence, then prefix  == suffix, and no mask left
            // in that case, we just use only prefix without mask 
            // More exceptiosn probably?
            if p.len() == s.len() && p.len() == full_len {
                Some(FlankSeq::with_flanks(Some(&p), None, 0))
            } else {
                let mask_size = full_len - p.len() - s.len();
                // println!("Mask size: {}", mask_size);
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
        (None, None) => None,
    }
}


pub fn read_queries(fasta_files: Vec<&str>, min_flank_length: Option<usize>) -> Vec<QueryGroup> {
    let mut query_groups = Vec::new();
    let min_flank_length = min_flank_length.unwrap_or(1);
    
    
    // For now we just support up to 2 fasta files, as single barcode and dual barcode
    // would it make sense to support more files?
    if fasta_files.len() > 2 {
        panic!("Only up to 2 fasta files are supported; if more needed please raise an issue");
    }

    // We just always call the first file forward, and the second file reverse barcodes
    let fasta_groups = vec![MatchType::Fbarcode, MatchType::Rbarcode];
    
    for (i, fasta_file) in fasta_files.iter().enumerate() {
        
        // Extract queries and get prefix/suffix flanks
        let (queries, ids) = read_fasta(fasta_file);
        let flank = get_flank_seq(&queries, min_flank_length);

        // Get match type 
        let match_type = fasta_groups[i].clone();

        // Collect sequences in forward (original) orientation  ----->
        let forward_queries = queries.iter().zip(&ids).map(|(seq, id)| Query {
            id: id.clone(),
            seq: seq.clone(),
        }).collect();
        query_groups.push(QueryGroup::new(flank.clone(), forward_queries, ids.clone(), match_type.clone(), Orientation::Forward));

        // Reverse complement group, keep the same id but reverse complement the sequence
        // same Fasta origin file so we keep the match type 
        let reverse_queries = queries.iter().zip(&ids).map(|(seq, id)| Query {
            id: id.clone(),
            seq: reverse_complement(seq),
        }).collect();
        let rc_flanks = flank.map(|f| f.reverse_complement());
        query_groups.push(QueryGroup::new(rc_flanks, reverse_queries, ids.clone(), match_type.clone(), Orientation::ReverseComplement));
    }

    println!("\n{}", "Query Groups".bold().underline());
    println!("{}:", "Found".bold());

    // Print statistics for each query group
    for g in &query_groups {
        let orientation = if g.orientation == Orientation::ReverseComplement { "RC".red() } else { "FW".green() };
        println!("  • Sequences: {} ({})", g.queries.len().to_string().bold(), orientation);
        
        if let Some(flank) = &g.flank {
            println!("    Pattern: {}", flank.as_string().dimmed());
        }
        println!();
    }

    query_groups
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_get_flanks() {
        let seqs = vec![
            b"AAAGGGGGTT".to_vec(), 
            b"AAGGGGGTTT".to_vec()
        //    __      __  <- AA and TT
        ];
        let (prefix, suffix, len) = get_flanks(&seqs, 1);
        assert_eq!(prefix, Some(b"AA".to_vec()));
        assert_eq!(suffix, Some(b"TT".to_vec()));
        assert_eq!(len, 10); // Length of first sequence
    }

    #[test]
    fn test_read_and_flanks() {
        let groups  = read_queries(vec!["data/test/left.fasta", "data/test/right.fasta"], None);
        
        let left_flank = groups[0].flank.as_ref().unwrap();
        let right_flank = groups[2].flank.as_ref().unwrap();

        let left_seq = left_flank.seq.clone();
        let right_seq = right_flank.seq.clone();

        // Left fasta
        /*
            >F_1
            AAAGGGTTTT
            >F_2
            AAGGGTTTTT
            __     ___
              NNNNN
            AA     TTT
        */

        assert_eq!(left_seq, b"AA****TTTT");

        // Right fasta
        /*
            >R_1
            AAAGGGGGTT
            >R_2
            AAGGGGGTTT
            __      __
              NNNNNN
            AA      TT
        */
        assert_eq!(right_seq, b"AA******TT");
        
    }

}
