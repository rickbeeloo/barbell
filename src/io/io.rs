use flate2::read::MultiGzDecoder;
use seq_io::fastq::Reader;
use std::fs::File;
use std::io::Read;

/// Open a FASTQ file, transparently decompressing gzip if the path ends in `.gz`
pub fn open_fastq(path: &str) -> Reader<Box<dyn Read + Send>> {
    let reader: Box<dyn Read + Send> = if path.to_ascii_lowercase().ends_with(".gz") {
        let file = File::open(path)
            .unwrap_or_else(|e| panic!("Failed to open gzip FASTQ file '{path}': {e}"));
        Box::new(MultiGzDecoder::new(file))
    } else {
        let file = File::open(path)
            .unwrap_or_else(|e| panic!("Failed to open FASTQ file '{path}': {e}"));
        Box::new(file)
    };
    Reader::new(reader)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use seq_io::fastq::Record;
    use std::io::Write;
    use tempfile::NamedTempFile;

    const FASTQ_CONTENT: &[u8] = b"@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTTAAAA\n+\nIIIIIIII\n";

    #[test]
    fn test_open_fastq_plain() {
        let mut tmp = NamedTempFile::with_suffix(".fastq").unwrap();
        tmp.write_all(FASTQ_CONTENT).unwrap();
        tmp.flush().unwrap();

        let mut reader = open_fastq(tmp.path().to_str().unwrap());
        let mut ids: Vec<String> = Vec::new();
        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            ids.push(record.id().unwrap().to_string());
        }
        assert_eq!(ids, vec!["read1", "read2"]);
    }

    #[test]
    fn test_open_fastq_gzip() {
        let mut tmp = NamedTempFile::with_suffix(".fastq.gz").unwrap();
        let mut encoder = GzEncoder::new(&mut tmp, Compression::default());
        encoder.write_all(FASTQ_CONTENT).unwrap();
        encoder.finish().unwrap();
        let path = tmp.path().to_str().unwrap().to_string();

        let mut reader = open_fastq(&path);
        let mut ids: Vec<String> = Vec::new();
        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            ids.push(record.id().unwrap().to_string());
        }
        assert_eq!(ids, vec!["read1", "read2"]);
    }
}
