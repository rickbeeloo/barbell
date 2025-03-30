use std::io;
use std::fs::File;
use std::path::Path;
use csv::{ReaderBuilder, WriterBuilder};
use tempfile::tempdir;
use super::filter::AnnotationLine;

const CHUNK_SIZE: usize = 10000;

pub fn merge_sort_files(input_file: &str) -> io::Result<()> {
    let temp_dir = tempdir()?;
    let mut chunk_files = Vec::new();
    
    // Create temporary output file path by appending "_sorted" to the input file
    let temp_output = format!("{}_sorted", input_file);

    // Step 1: Split into sorted chunks
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(input_file)?;

    // Add debug print for input file
    println!("Reading from input file: {}", input_file);
    
    // Read and store the header
    let mut headers = reader.headers()?.clone();
    println!("Headers read: {:?}", headers);
    
    if !headers.iter().any(|h| h == "cuts") {
        headers.push_field("cuts");
    }

    let mut chunk_count = 0;
    let mut total_records = 0;
    loop {
        let mut chunk: Vec<AnnotationLine> = reader
            .deserialize()
            .take(CHUNK_SIZE)
            .collect::<Result<Vec<_>, _>>()?;
 
        if chunk.is_empty() {
            break;
        }

        total_records += chunk.len();
        println!("Processing chunk {} with {} records", chunk_count, chunk.len());

        chunk.sort_by_key(|r| (r.record_set_idx, r.record_idx));

        let chunk_path = temp_dir.path().join(format!("chunk_{}.tsv", chunk_count));
        let mut chunk_writer = WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path(&chunk_path)?;

        for record in &chunk {
            chunk_writer.write_record(&record.to_record())?;
        }
        chunk_writer.flush()?;
        
        chunk_files.push(chunk_path);
        chunk_count += 1;
    }

    println!("Total records processed: {}", total_records);
    println!("Number of chunks created: {}", chunk_count);

    // Step 2: Merge sorted chunks with header into temporary file
    merge_chunks(&chunk_files, &temp_output, headers)?;

    // Verify temp_output exists and has content
    if let Ok(metadata) = std::fs::metadata(&temp_output) {
        println!("Temporary output file size: {} bytes", metadata.len());
    }

    // Step 3: Only remove original file if temporary file exists and has content
    if std::path::Path::new(&temp_output).exists() {
        std::fs::remove_file(input_file)?;
        println!("Original file removed successfully");
        std::fs::copy(&temp_output, input_file)?;
        println!("New sorted file copied successfully");
        
        // Verify the final file
        if let Ok(metadata) = std::fs::metadata(input_file) {
            println!("Final sorted file size: {} bytes", metadata.len());
        }
    } else {
        return Err(io::Error::other(
            "Temporary output file was not created successfully"
        ));
    }

    Ok(())
}

fn merge_chunks<P: AsRef<Path>>(chunk_files: &[P], output_file: &str, headers: csv::StringRecord) -> io::Result<()> {
    let mut readers: Vec<csv::Reader<File>> = chunk_files
        .iter()
        .map(|path| {
            ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .from_path(path)
                .unwrap()
        })
        .collect();

    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output_file)?;

    // Write the header to the output file
    writer.write_record(&headers)?;

    let mut current_records: Vec<Option<AnnotationLine>> = vec![None; readers.len()];
    
    // Initialize with first records
    for (i, reader) in readers.iter_mut().enumerate() {
        current_records[i] = reader.deserialize().next().and_then(|r| r.ok());
    }

    // Merge chunks
    while current_records.iter().any(|r| r.is_some()) {
        let min_idx = current_records
            .iter()
            .enumerate()
            .filter(|(_, r)| r.is_some())
            .min_by(|(_, a), (_, b)| {
                let a = a.as_ref().unwrap();
                let b = b.as_ref().unwrap();
                (a.record_set_idx, a.record_idx).cmp(&(b.record_set_idx, b.record_idx))
            })
            .map(|(idx, _)| idx);

        if let Some(idx) = min_idx {
            if let Some(record) = &current_records[idx] {
                writer.write_record(&record.to_record())?;
            }

            current_records[idx] = readers[idx].deserialize().next().and_then(|r| r.ok());
        }
    }

    writer.flush()?;
    Ok(())
}

