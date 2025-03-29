use crate::barbell::filter::AnnotationLine;
use crate::barbell::merge_sort::merge_sort_files;
use seq_io::fastq::{Reader,Record};
use csv;
use std::io::{BufWriter, Write};
use std::fs::File;




pub fn get_cuts(annotations: &[AnnotationLine]) {
    // Easier would be specifying the cuts in the filter, but that's a bit less clean cause then the filtered output file 
    // is different from the input file. For now lets try to get the biggest fragment

}


pub fn process_read_and_anno(seq: &[u8], qual: &[u8], annotations: &[AnnotationLine]) -> (Vec<u8>, Vec<u8>) {
    let cuts = annotations.first().unwrap().cut_positions.clone().unwrap_or_default();
    
    match cuts.len() {
        0 => (seq.to_vec(), qual.to_vec()),
        1 => {
            let cut_pos = cuts[0];
            let seq_len = seq.len();
            
            // Determine if cut is closer to left or right end
            if cut_pos <= seq_len / 2 {
                // Cut is closer to left end, keep right part
                (seq[cut_pos..].to_vec(), qual[cut_pos..].to_vec())
            } else {
                // Cut is closer to right end, keep left part
                (seq[..cut_pos].to_vec(), qual[..cut_pos].to_vec())
            }
        },
        2 => {
            let (start, end) = (cuts[0], cuts[1]);
            if end <= start {
                eprintln!("Warning: Invalid cut positions (end <= start): {} <= {}", end, start);
                return (seq.to_vec(), qual.to_vec());
            }
            (seq[start..end].to_vec(), qual[start..end].to_vec())
        },
        _ => {
            eprintln!("Warning: More than 2 cut positions specified: {:?}", cuts);
            (seq.to_vec(), qual.to_vec())
        }
    }
}

pub fn trim_matches(filtered_match_file: &str, read_fastq_file: &str, output_file: &str) {

    let mut reader = Reader::from_path(read_fastq_file).unwrap();
    let writer = BufWriter::new(File::create(format!("{}.trimmed.fastq", output_file))
        .expect("Failed to create output file"));
    let mut writer = writer;
    let mut total_reads = 0;
    let mut mapped_reads = 0;
    let mut trimmed_reads = 0;
    
    let mut matches_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(filtered_match_file)
        .expect("Failed to open matches file");

    // Group annotations by read ID while preserving order
    let mut current_read_id: Option<Vec<u8>> = None;
    let mut current_annotations: Vec<AnnotationLine> = Vec::new();

    for result in matches_reader.deserialize() {
        let anno: AnnotationLine = result.expect("Failed to parse annotation line");
        let read_id = anno.read.clone().into_bytes();

        // If we encounter a new read ID
        if current_read_id.as_ref() != Some(&read_id) {
            // Process previous group if it exists
            if let Some(prev_read_id) = current_read_id.take() {
                // Search through FASTQ records until we find matching read
                while let Some(record) = reader.next() {
                    let record = record.expect("Error reading record");
                    total_reads += 1;
                    
                    if record.id().unwrap().as_bytes() == prev_read_id.as_slice() {
                        mapped_reads += 1;
                        let (trimmed_seq, trimmed_qual) = process_read_and_anno(
                            record.seq(),
                            record.qual(),
                            &current_annotations
                        );
                        
                        // Only write if we actually trimmed something
                        if trimmed_seq.len() != record.seq().len() {
                            trimmed_reads += 1;
                            // Write FASTQ format: @header\nsequence\n+\nquality
                            writeln!(writer, "@{}", String::from_utf8_lossy(record.id().unwrap().as_bytes()))
                                .expect("Failed to write header");
                            writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_seq))
                                .expect("Failed to write sequence");
                            writeln!(writer, "+").expect("Failed to write separator");
                            writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_qual))
                                .expect("Failed to write quality");
                        }
                        break;
                    }
                }
            }
            
            // Start new group
            current_read_id = Some(read_id);
            current_annotations.clear();
        }
        
        current_annotations.push(anno);
    }

    // Process the last group
    if let Some(last_read_id) = current_read_id {
        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            total_reads += 1;
            
            if record.id().unwrap().as_bytes() == last_read_id.as_slice() {
                mapped_reads += 1;
                let (trimmed_seq, trimmed_qual) = process_read_and_anno(
                    record.seq(),
                    record.qual(),
                    &current_annotations
                );
                
                // Only write if we actually trimmed something
                if trimmed_seq.len() != record.seq().len() {
                    trimmed_reads += 1;
                    // Write FASTQ format: @header\nsequence\n+\nquality
                    writeln!(writer, "@{}", String::from_utf8_lossy(record.id().unwrap().as_bytes()))
                        .expect("Failed to write header");
                    writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_seq))
                        .expect("Failed to write sequence");
                    writeln!(writer, "+").expect("Failed to write separator");
                    writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_qual))
                        .expect("Failed to write quality");
                }
                break;
            }
        }
    }

    // Ensure all data is written
    writer.flush().expect("Failed to flush output");

    println!("Total reads processed: {}", total_reads);
    println!("Reads mapped to annotations: {}", mapped_reads);
    println!("Reads trimmed: {}", trimmed_reads);
    println!("Mapping rate: {:.2}%", (mapped_reads as f64 / total_reads as f64) * 100.0);
    println!("Trimming rate: {:.2}%", (trimmed_reads as f64 / mapped_reads as f64) * 100.0);
}