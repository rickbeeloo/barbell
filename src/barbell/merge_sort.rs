use std::fs::{File, remove_file, copy};
use std::io;
use std::path::{Path, PathBuf};

use csv::{ReaderBuilder, WriterBuilder};
use tempfile::tempdir;

use super::filter::AnnotationLine;

const CHUNK_SIZE: usize = 1000;

pub fn merge_sort_files(input_file: &str) {
    let temp_dir = tempdir().expect("Failed to create temp directory");
    let mut chunk_files: Vec<PathBuf> = Vec::new();

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(input_file)
        .expect("Failed to open input file");

    let headers = reader.headers().expect("Failed to read headers").clone();
    let mut chunk_count = 0;

    loop {
        let chunk: Vec<AnnotationLine> = reader
            .deserialize()
            .take(CHUNK_SIZE)
            .collect::<Result<_, _>>()
            .expect("Failed to deserialize chunk");

        if chunk.is_empty() {
            break;
        }

        let mut sorted_chunk = chunk;
        sorted_chunk.sort_by_key(|r| (r.record_set_idx, r.record_idx));

        let chunk_path = temp_dir.path().join(format!("chunk_{}.tsv", chunk_count));
        let mut writer = WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(&chunk_path)
            .expect("Failed to create chunk file");

        writer.write_record(&headers).expect("Failed to write header");
        for record in &sorted_chunk {
            writer.write_record(&record.to_record()).expect("Failed to write record");
        }

        writer.flush().expect("Failed to flush writer");

        chunk_files.push(chunk_path);
        chunk_count += 1;
    }

    let temp_output = format!("{}_sorted", input_file);
    merge_chunks(&chunk_files, &temp_output, headers).expect("Failed to merge chunks");

    // Replace original with sorted
    remove_file(input_file).expect("Failed to remove original file");
    copy(&temp_output, input_file).expect("Failed to copy sorted file back");
    remove_file(temp_output).expect("Failed to remove temporary sorted file");
}

fn merge_chunks<P: AsRef<Path>>(
    chunk_files: &[P],
    output_file: &str,
    headers: csv::StringRecord,
) -> io::Result<()> {
    assert!(!chunk_files.is_empty(), "No chunks to merge");

    let mut readers: Vec<csv::Reader<File>> = chunk_files
        .iter()
        .map(|path| {
            ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(true)
                .from_path(path)
                .expect("Failed to open chunk file")
        })
        .collect();

    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output_file)
        .expect("Failed to create output file");

    writer.write_record(&headers)?;

    let mut current_records: Vec<Option<AnnotationLine>> = readers
        .iter_mut()
        .map(|r| r.deserialize::<AnnotationLine>().next().transpose())
        .collect::<Result<_, _>>()
        .expect("Failed to read initial records from chunks");

    while current_records.iter().any(|r| r.is_some()) {
        let min_idx = current_records
            .iter()
            .enumerate()
            .filter_map(|(i, r)| r.as_ref().map(|_| i))
            .min_by_key(|&i| {
                let record = current_records[i].as_ref().unwrap();
                (record.record_set_idx, record.record_idx)
            })
            .expect("Failed to find next record");

        if let Some(record) = current_records[min_idx].take() {
            writer.write_record(&record.to_record())?;
            current_records[min_idx] = readers[min_idx]
                .deserialize::<AnnotationLine>()
                .next()
                .transpose()
                .expect("Failed to deserialize next record");
        }
    }

    writer.flush()?;
    Ok(())
}