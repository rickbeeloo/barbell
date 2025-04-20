use seq_io::fasta::{Reader,Record};
use colored::Colorize;
use crate::types::{MatchType, Orientation};
use crate::annotate::flank::FlankGroup;

#[derive(Debug, Clone)]
pub struct Query {
    pub id: String,
    pub seq: Vec<u8>,
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


pub fn read_queries(fasta_files: Vec<&str>) -> Vec<FlankGroup> {
    let mut flank_queries = Vec::new();
    
    // For now we just support up to 2 fasta files, as single barcode and dual barcode
    // would it make sense to support more files?
    if fasta_files.len() > 2 {
        panic!("Only up to 2 fasta files are supported; if more needed please raise an issue");
    }

    // We just always call the first file forward, and the second file reverse barcodes
    let fasta_groups = [MatchType::Fbarcode, MatchType::Rbarcode];
    
    for (i, fasta_file) in fasta_files.iter().enumerate() {
        
        // Read the sequences from user file
        let (queries, ids) = read_fasta(fasta_file);

        // Transform in Flank group, forward
        let match_type = fasta_groups[i].clone();
        let forward_group: FlankGroup = FlankGroup::new(queries, ids, Orientation::Forward, match_type);
        flank_queries.push(forward_group.clone());

        // Also add the reverse complement orientation
        let rc_group = forward_group.reverse_complement();
        flank_queries.push(rc_group);

    }

    println!("\n{}", "Query Groups".bold().underline());
    println!("{}:", "Found".bold());

    // Print statistics for each query group
    for g in &flank_queries {
        let orientation = if g.orientation == Orientation::ReverseComplement { "RC".red() } else { "FW".green() };
        println!("  • Sequences: {} ({})", g.flank_seq.mask_queries.len().to_string().bold(), orientation);
        println!("    Pattern: {}", g.flank_seq.as_string().dimmed());
        
        println!();
    }

    flank_queries
}