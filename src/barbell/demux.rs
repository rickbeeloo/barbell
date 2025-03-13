// use crate::barbell::reader::*;
// use crate::barbell::seq::*;
// use seq_io::fastq::{Reader as FastqReader, Record}; // both Fasta and Fastq import same name struct
// use crate::barbell::strategy::*;
// use std::io::Write;

// pub struct Demuxer{
//     query_groups: Vec<QueryGroup>,
//     q: f64, 
//     strategy: ScoreDiffStrategy,
// }

// impl Demuxer {
//     pub fn new(query_groups: Vec<QueryGroup>) -> Self {
//         Self { 
//             query_groups,
//             q: 0.5,           
//             strategy: ScoreDiffStrategy::default(),
//         }
//     }

//     // At least this is verbose for init as lib, but maybe not consumer pattern later
//     pub fn with_q(mut self, q: f64) -> Self {
//         self.q = q;
//         self
//     }

//     pub fn with_strategy(mut self, strategy: ScoreDiffStrategy) -> Self {
//         self.strategy = strategy;
//         self
//     }


//     fn rel_dist_to_end(pos: isize, read_len: usize) -> isize {
//         // If already negative, for a match starting before the read we can return 1
//         if pos < 0 {
//             return 1;
//         }
    
//         if pos <= (read_len / 2) as isize {
//             if pos == 0 {
//                 1  // Left end
//             } else {
//                 pos  // Distance from left end (already isize)
//             }
//         } else if pos == read_len as isize {
//             -1  // Right end
//         } else {
//             -(read_len as isize - pos)  // Distance from right end
//         }
//     }


//     pub fn process_reads(&self, read_file: &str)  {
//         let mut reader = FastqReader::from_path(read_file).expect("Failed to open read file");
//         let mut count = 0;
//         let mut chunk_count = 0;
        
//         // Skip the first 50 sequences
//         while let Some(record) = reader.next() {
//             count += 1;
//             if count >= 50 {
//                 break;
//             }
//         }
        
//         // Process remaining sequences in chunks of 50
//         let mut done = 0;
//         let mut assigned = 0;
//         let output_file = std::fs::File::create("output.txt").expect("Failed to create output file");
//         let mut writer = std::io::BufWriter::new(output_file);
        
//         while let Some(record) = reader.next() {
//             let record = record.expect("Error reading record");
//             let seq = record.seq();
//             let id = record.id().expect("Invalid fastq id");

//             if id != "f600409d-1224-4fa9-8052-372f35cab81f" {
//                 continue;
//             }

//             let annotation = self.strategy.annotate(&self.query_groups, seq);
//             let assignments = self.strategy.decide(&annotation);
                     
//             // for asn in assignments.iter() {
//             //     println!("{:?}", asn);
//             // }

//             if assignments.len() == 2 {
//                 assigned += 1;
//                 for assignment in assignments {
//                     // Start relative to read end:
//                     writeln!(
//                         writer,
//                         "{}\t{}\t{}\t{}",
//                         id,
//                         assignment.query_id,
//                         Demuxer::rel_dist_to_end(assignment.region_start as isize, seq.len()),
//                         assignment.dist
//                     ).expect("Failed to write to output file");
//                 }
//             }
            
//             done += 1;
//             if done % 1000 == 0 {
//                 println!("Processed {} reads, assigned {} reads", done, assigned);
//             }
//         }
        
//         writer.flush().expect("Failed to flush writer");
//     }
// }


