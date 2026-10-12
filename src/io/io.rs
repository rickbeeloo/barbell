use paraseq::BoxedReader;
use paraseq::fastx::{Collection, CollectionType};
use std::path::PathBuf;

/// Split a FASTQ record header into the read ID and optional description
pub fn split_fastq_header(header: &[u8]) -> anyhow::Result<(&str, &str)> {
    let header = std::str::from_utf8(header)
        .map_err(|err| anyhow::anyhow!("Input FASTQ read id is not valid UTF-8: {err}"))?;
    if let Some(split_idx) = header.find(char::is_whitespace) {
        let read_id = &header[..split_idx];
        let desc = header[split_idx..].trim_start();
        Ok((read_id, desc))
    } else {
        Ok((header, ""))
    }
}

/// Validate FASTQ input paths passed by the CLI or public API.
pub fn validate_fastq_paths(paths: &[PathBuf]) -> anyhow::Result<()> {
    if paths.is_empty() {
        anyhow::bail!("No FASTQ input files provided");
    }

    Ok(())
}

/// Open one or more FASTQ files as a paraseq single-read collection.
pub fn open_fastq_collection(paths: &[PathBuf]) -> anyhow::Result<Collection<BoxedReader>> {
    validate_fastq_paths(paths)?;
    Collection::from_paths(&paths, CollectionType::Single)
        .map_err(|e| anyhow::anyhow!("Failed to open FASTQ input: {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_validate_fastq_plain() {
        let mut tmp = NamedTempFile::with_suffix(".fastq").unwrap();
        tmp.write_all(b"@read1\nACGT\n+\nIIII\n").unwrap();
        tmp.flush().unwrap();

        let paths = vec![tmp.path().to_path_buf()];
        validate_fastq_paths(&paths).unwrap();
    }

    #[test]
    fn test_validate_fastq_gzip() {
        let mut tmp = NamedTempFile::with_suffix(".fastq.gz").unwrap();
        tmp.write_all(b"not actually compressed").unwrap();

        let paths = vec![tmp.path().to_path_buf()];
        validate_fastq_paths(&paths).unwrap();
    }

    #[test]
    fn test_validate_fastq_empty_errors() {
        assert!(validate_fastq_paths(&[]).is_err());
    }

    #[test]
    fn test_split_fastq_header_with_description() {
        let (read_id, desc) = split_fastq_header(b"read1 some description").unwrap();
        assert_eq!(read_id, "read1");
        assert_eq!(desc, "some description");
    }

    #[test]
    fn test_split_fastq_header_without_description() {
        let (read_id, desc) = split_fastq_header(b"read1").unwrap();
        assert_eq!(read_id, "read1");
        assert_eq!(desc, "");
    }
}
