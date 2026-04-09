use flate2::read::MultiGzDecoder;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_bam::Record as BamRecord;
use seq_io::fastq::{Record as SeqIoRecord, Reader};
use std::fs::File;
use std::io::Read as IoRead;

pub trait Record {
    fn id(&self) -> &str;
    fn seq(&self) -> &[u8];
    fn qual(&self) -> Option<&[u8]>;

    fn len(&self) -> usize {
        self.seq().len()
    }

    fn is_empty(&self) -> bool {
        self.seq().is_empty()
    }
}

pub trait FormatReader {
    type Rec: Record;

    fn next_record(&mut self) -> Option<anyhow::Result<Self::Rec>>;
}

/// Open a FASTQ file, decompressing gzip if the path ends in `.gz`
pub fn open_fastq(path: &str) -> Reader<Box<dyn IoRead + Send>> {
    let reader: Box<dyn IoRead + Send> = if path.to_ascii_lowercase().ends_with(".gz") {
        let file = File::open(path)
            .unwrap_or_else(|e| panic!("Failed to open gzip FASTQ file '{path}': {e}"));
        Box::new(MultiGzDecoder::new(file))
    } else {
        let file =
            File::open(path).unwrap_or_else(|e| panic!("Failed to open FASTQ file '{path}': {e}"));
        Box::new(file)
    };
    Reader::new(reader)
}

/// The supported input source types.
pub enum ReadSource {
    Fastq(Reader<Box<dyn IoRead + Send>>),
    Bam(bam::io::Reader<bgzf::io::Reader<File>>),
}

/// A generic read record used by annotation and trimming.
pub struct ReadRecord {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub original_bam_record: Option<BamRecord>,
}

impl Record for ReadRecord {
    fn id(&self) -> &str {
        &self.id
    }

    fn seq(&self) -> &[u8] {
        &self.seq
    }

    fn qual(&self) -> Option<&[u8]> {
        Some(&self.qual)
    }
}

impl ReadSource {
    pub fn next_record(&mut self) -> Option<anyhow::Result<ReadRecord>> {
        match self {
            ReadSource::Fastq(reader) => {
                let record = reader.next()?;
                match record {
                    Ok(record) => {
                        let (id, _) = record.id_desc().unwrap();
                        Some(Ok(ReadRecord {
                            id: id.to_string(),
                            seq: record.seq().to_vec(),
                            qual: record.qual().to_vec(),
                            original_bam_record: None,
                        }))
                    }
                    Err(err) => Some(Err(anyhow::anyhow!("FASTQ read failed: {err}"))),
                }
            }
            ReadSource::Bam(reader) => {
                let mut record = BamRecord::default();
                match reader.read_record(&mut record) {
                    Ok(0) => None,
                    Ok(_) => Some(Ok(ReadRecord {
                        id: record
                            .name()
                            .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
                            .unwrap_or_default(),
                        seq: record.sequence().iter().collect(),
                        qual: record.quality_scores().as_ref().to_vec(),
                        original_bam_record: Some(record),
                    })),
                    Err(err) => Some(Err(anyhow::anyhow!("BAM read failed: {err}"))),
                }
            }
        }
    }
}

impl FormatReader for ReadSource {
    type Rec = ReadRecord;

    fn next_record(&mut self) -> Option<anyhow::Result<Self::Rec>> {
        ReadSource::next_record(self)
    }
}

/// Open a generic read source: FASTQ, FASTQ.GZ, BAM or SAM.
pub fn open_reads(path: &str) -> anyhow::Result<ReadSource> {
    if path_is_bam(path) {
        let file = File::open(path)
            .map_err(|err| anyhow::anyhow!("Failed to open BAM file '{path}': {err}"))?;
        let reader = bam::io::Reader::new(file);
        Ok(ReadSource::Bam(reader))
    } else {
        Ok(ReadSource::Fastq(open_fastq(path)))
    }
}

fn path_is_bam(path: &str) -> bool {
    let path = path.to_ascii_lowercase();
    path.ends_with(".bam")
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use noodles_sam as sam;
    use noodles_sam::alignment::io::Write as _;
    use std::io::Write;
    use tempfile::{tempdir, NamedTempFile};

    const FASTQ_CONTENT: &[u8] = b"@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTTAAAA\n+\nIIIIIIII\n";

    #[test]
    fn test_open_fastq_plain() {
        let mut tmp = NamedTempFile::with_suffix(".fastq").unwrap();
        tmp.write_all(FASTQ_CONTENT).unwrap();
        tmp.flush().unwrap();

        let mut reader = open_fastq(tmp.path().to_str().unwrap());
        let mut ids = Vec::new();
        let mut seqs = Vec::new();
        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            ids.push(record.id().unwrap().to_string());
            seqs.push(String::from_utf8_lossy(record.seq()).to_string());
        }
        assert_eq!(ids, vec!["read1", "read2"]);
        assert_eq!(seqs, vec!["ACGTACGT", "TTTTAAAA"]);
    }

    #[test]
    fn test_open_fastq_gzip() {
        let mut tmp = NamedTempFile::with_suffix(".fastq.gz").unwrap();
        let mut encoder = GzEncoder::new(&mut tmp, Compression::default());
        encoder.write_all(FASTQ_CONTENT).unwrap();
        encoder.finish().unwrap();
        let path = tmp.path().to_str().unwrap().to_string();

        let mut reader = open_fastq(&path);
        let mut ids = Vec::new();
        let mut seqs = Vec::new();
        while let Some(record) = reader.next() {
            let record = record.expect("Error reading record");
            ids.push(record.id().unwrap().to_string());
            seqs.push(String::from_utf8_lossy(record.seq()).to_string());
        }
        assert_eq!(ids, vec!["read1", "read2"]);
        assert_eq!(seqs, vec!["ACGTACGT", "TTTTAAAA"]);
    }

    #[test]
    fn test_path_is_bam() {
        assert!(path_is_bam("reads.bam"));
        assert!(path_is_bam("READS.BAM"));
        assert!(!path_is_bam("reads.sam"));
        assert!(!path_is_bam("reads.fastq"));
    }

    #[test]
    fn test_open_reads_bam() {
        let tmpdir = tempdir().unwrap();
        let path = tmpdir.path().join("test.bam");
        let file = File::create(&path).unwrap();
        let mut writer = bam::io::Writer::new(file);

        let header = sam::Header::default();
        writer.write_header(&header).unwrap();

        let record = sam::alignment::RecordBuf::builder()
            .set_name("read1".to_string())
            .set_flags(sam::alignment::record::Flags::UNMAPPED)
            .set_sequence(b"ACGT".to_vec().into())
            .set_quality_scores(b"IIII".to_vec().into())
            .build();

        writer.write_alignment_record(&header, &record).unwrap();
        writer.try_finish().unwrap();
        drop(writer);

        let mut source = open_reads(path.to_str().unwrap()).unwrap();
        if let ReadSource::Bam(ref mut bam_reader) = source {
            let _header = bam_reader
                .read_header()
                .expect("Failed to read BAM header from generated BAM file");
        } else {
            panic!("expected BAM read source");
        }

        let rec = source.next_record().unwrap().unwrap();
        assert_eq!(rec.id, "read1");
        assert_eq!(rec.seq, b"ACGT");
        assert_eq!(rec.qual, b"IIII");
        assert!(rec.original_bam_record.is_some());
    }
}
