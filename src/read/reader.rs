use crate::annotate::flank::FlankGroup;
use crate::types::{MatchType, Orientation};
use colored::Colorize;
use seq_io::fasta::{Reader, Record};

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

pub fn group_names_to_match_types(group_names: Vec<String>) -> Vec<MatchType> {
    group_names
        .iter()
        .map(|name| {
            if name == "Fbar" {
                MatchType::Fbarcode
            } else if name == "Rbar" {
                MatchType::Rbarcode
            } else {
                panic!("Invalid group name: {} use Fbar or Rbar", name);
            }
        })
        .collect()
}

pub fn read_queries(fasta_files: Vec<&str>, group_names: Vec<String>) -> Vec<FlankGroup> {
    let mut flank_queries = Vec::new();

    if fasta_files.len() != group_names.len() {
        panic!("Number of fasta files and group names must be the same");
    }
    let group_types = group_names_to_match_types(group_names);

    for (i, fasta_file) in fasta_files.iter().enumerate() {
        // Read the sequences from user file
        let (queries, ids) = read_fasta(fasta_file);

        // Transform in Flank group, forward
        let match_type = group_types[i].clone();
        let forward_group: FlankGroup =
            FlankGroup::new(queries, ids, Orientation::Forward, match_type);
        flank_queries.push(forward_group.clone());

        // Also add the reverse complement orientation
        let rc_group = forward_group.reverse_complement();
        flank_queries.push(rc_group);
    }

    println!("\n{}", "Query Groups".bold().underline());
    println!("{}:", "Found".bold());

    // Print statistics for each query group
    for g in &flank_queries {
        let orientation = if g.orientation == Orientation::ReverseComplement {
            "RC".red()
        } else {
            "FW".green()
        };
        println!(
            "  • Sequences: {} ({})",
            g.flank_seq.mask_queries.len().to_string().bold(),
            orientation
        );
        println!("    Pattern: {}", g.flank_seq.as_string().dimmed());

        println!();
    }

    flank_queries
}
