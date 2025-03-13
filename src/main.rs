pub mod barbell;

use barbell::reader::*;
use barbell::demux::*;
use std::time::Instant;
use std::io::Write;

use crate::barbell::seq::*;
use seq_io::fastq::{Reader as FastqReader, Record}; // both Fasta and Fastq import same name struct
use crate::barbell::strategy::*;



pub fn process_reads(read_file: &str, queries: &[QueryGroup])  {
    let mut reader = FastqReader::from_path(read_file).expect("Failed to open read file");
    let mut count = 0;
    let mut chunk_count = 0;
    
    // // Skip the first 50 sequences
    // while let Some(record) = reader.next() {
    //     count += 1;
    //     if count >= 50 {
    //         break;
    //     }
    // }
    
    // Process remaining sequences in chunks of 50
    let mut done = 0;
    let mut tag_found = 0;
    let output_file = std::fs::File::create("output_200k_testing4.txt").expect("Failed to create output file");
    let mut writer = std::io::BufWriter::new(output_file);

    writeln!(writer, "read\tlabel\tstart\tend\tedits\tdist.to.end\tread.len").expect("Failed to write header");
    
    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        let seq = record.seq();
        let id = record.id().expect("Invalid fastq id");

        if id != "35d261dc-664e-4c8a-8614-ed1bd401043b" {
            continue;
        } else {
            println!("Processing read: {}", id);
        }

        
        let read = Read::from_seq(seq);
        let mut read = read.with_strategy(SimpleStrategy::default());
        
        let annotations = read.annotate(queries);
        let final_tags = read.final_assignment(&annotations);

        if final_tags.len() > 0 {
            tag_found += 1;
        }
        
        // Write matches to output file
        for tag in final_tags.iter() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                id, tag.some_id, tag.start, tag.end, tag.edits, tag.rel_dist_to_end, seq.len()
            ).expect("Failed to write to output file");
        }

        // If no tag is found lets write a line still but with no tag   
        if final_tags.len() == 0 {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                id, "no_tag", 0, 0, 0, 0, seq.len()
            ).expect("Failed to write to output file");
        }
        
        done += 1;
        if done % 10_000 == 0 {
            println!("Tags found {} out of {} reads", tag_found, done);
        }
    }
    
    println!("Done! Tags found {} out of {}", tag_found, done);
}


fn main() {
    println!("Starting Barbell");

     let groups =  read_queries(vec!["data/rapid_barcodes.fasta"], None);
     process_reads("/home/solprof/PhD/barbell-sg/data/pass_combined.fastq", &groups);



    // let demux = Demuxer::new(groups);   

    // demux.process_reads("data/long_reads.fastq");

}
