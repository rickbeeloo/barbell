//! Adapter trimming for Oxford Nanopore ligation-based (kit14) libraries.
//!
//! Reads from kit14 libraries (e.g. `SQK-LSK114`, `SQK-ULK114`, and the other
//! "114" kits) carry the same ligation adapter at both ends of the insert:
//!
//! ```text
//!   5' [adapter] ------------------ insert ------------------ [adapter] 3'
//!                \                                            /
//!                 both are the same double-stranded AMX114 adapter
//! ```
//!
//! Depending on which strand was captured, the first bases of the read are either
//! the adapter itself or its reverse complement:
//!
//! * forward adapter  : `CCTGTACTTCGTTCAGTTACGTATTGC`
//! * reverse adapter  : `GCAATACGTAACTGAACGAAGTACAGG`
//!
//! (These are the sequences Dorado uses to trim `SQK-LSK114`, `SQK-ULK114`,
//! `SQK-PCS114`, `SQK-RAD114`, `SQK-16S114-24`, `SQK-NBD114-*`, `SQK-PCB114-24`,
//! `SQK-RBK114-*`, `SQK-RPB114-24`, `SQK-HTB114-96`, ... see
//! `dorado/demux/adapter_primer_kits.cpp`.)
//!
//! This module finds (partial) adapter matches at the very start/end of a read and
//! reports the interval of the read that should be kept, in the same spirit as
//! [Porechop](https://github.com/rrwick/Porechop). The read is *not* reverse
//! complemented: trimming is always relative to the left (5') end of the read as
//! it is stored in the FASTQ.
//!
//! Adapters in the middle of a read are deliberately ignored, because those cannot
//! be resolved by plain end trimming (that requires splitting/filtering, which is
//! out of scope here).

use anyhow::anyhow;
use paraseq::Record;
use rayon::prelude::*;
use sassy::profiles::Iupac;
use sassy::{EncodedPatterns, Match, Searcher};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::config::AdapterTrimConfig;
use crate::io::io::validate_fastq_paths;

/// The forward (5') AMX114 ligation adapter used by all kit14 sequencing kits.
pub const ADAPTER_SEQUENCE: &[u8] = b"CCTGTACTTCGTTCAGTTACGTATTGC";

/// An adapter is only considered to sit at a read end when its outer edge is
/// within this many bases of the very start/end of the read. This tolerates the
/// few noisy bases that basecalling can leave in front of the adapter (Porechop
/// uses a similar slack via its `extra_end_trim`).
pub const ADAPTER_END_SLACK: usize = 10;

/// Overhang cost used for the adapter search. With `alpha = 0.5` a partial
/// alignment that only covers half of the adapter costs as much as a single
/// mismatch, which is what we want: a confidently matched adapter *sub*-sequence
/// at the read end is still evidence of an adapter.
const ADAPTER_SEARCH_ALPHA: f32 = 0.5;

/// Maximum number of edits allowed while searching for an adapter hit.
/// Together with the overhang cost this still allows partial hits (e.g. 12 bases
/// with 3 mismatches, or 8 perfectly matching bases), while keeping the search
/// selective.
const ADAPTER_MAX_EDITS: usize = 12;

/// A hit needs at least this many matching bases ...
///
/// This was tuned on simulated kit14 reads: with 14 matching bases and 75%
/// identity ~98.5% of reads carrying an adapter are trimmed while random reads
/// are only touched in <0.1% of the cases. Lowering it picks up a few more
/// (shorter) partial adapters at the cost of more false positives on random
/// sequence.
pub const ADAPTER_MIN_MATCH_LEN: usize = 14;

/// ... and at least this fraction of the alignment needs to match (gaps count as
/// mismatches).
pub const ADAPTER_MIN_IDENTITY: f64 = 0.75;

/// A trim that would leave less than this many bases is rejected (the read is
/// returned unchanged), to avoid emptying reads with many spurious matches.
const ADAPTER_MIN_INSERT_LEN: usize = 1;

/// A chained adapter has to be anchored within this many bases of the previous
/// cut. It allows a few bases in between while still rejecting an independent
/// adapter that merely happens to be close by.
const ADAPTER_CHAIN_OVERLAP: usize = 10;

/// Safety bound on the number of adapters trimmed in tandem at one read end.
const MAX_ADAPTERS_PER_END: usize = 8;

/// Which read end an adapter match was found on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdapterSide {
    /// The adapter (or its reverse complement) at the 5' end of the read.
    Start,
    /// The adapter (or its reverse complement) at the 3' end of the read.
    End,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdapterTrim {
    /// First base that is kept.
    pub keep_start: usize,
    /// One past the last base that is kept.
    pub keep_end: usize,
    /// Adapters that were found (left to right). Empty means nothing to trim.
    pub sides: [bool; 2],
}

impl AdapterTrim {
    pub fn nothing(read_len: usize) -> Self {
        Self {
            keep_start: 0,
            keep_end: read_len,
            sides: [false, false],
        }
    }

    pub fn is_empty(&self) -> bool {
        !self.sides[0] && !self.sides[1]
    }
}

/// Searches (partial) ligation adapters at the ends of Nanopore reads.
///
/// The internal search buffers make this cheap to reuse across reads, but it is
/// not `Sync`; instantiate one per thread.
pub struct AdapterTrimmer {
    searcher: Searcher<Iupac>,
    encoded: EncodedPatterns<Iupac>,
    min_identity: f64,
    min_match_len: usize,
    end_slack: usize,
}

impl Default for AdapterTrimmer {
    fn default() -> Self {
        Self::new(AdapterTrimConfig::default())
    }
}

impl AdapterTrimmer {
    pub fn new(config: AdapterTrimConfig) -> Self {
        let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(ADAPTER_SEARCH_ALPHA);
        // Pattern 0 is the adapter as found at the 5' end of a read, pattern 1 is
        // its reverse complement (which is what the 3' end looks like when the
        // other strand was captured).
        let fwd = ADAPTER_SEQUENCE.to_vec();
        let rc = reverse_complement(ADAPTER_SEQUENCE);
        let encoded = searcher.encode_patterns(&[fwd, rc]);
        Self {
            searcher,
            encoded,
            min_identity: config.min_identity,
            min_match_len: config.min_match_len,
            end_slack: config.end_slack,
        }
    }

    /// Find the adapters at the ends of `read` and return the interval that should
    /// be kept. `None` when there is nothing to trim.
    ///
    /// Only *terminal* adapters are considered, i.e. adapters that reach the very
    /// start/end of the read (allowing a few noisy bases). Adapters found in the
    /// middle of a read are ignored: those cannot be resolved by end trimming.
    ///
    /// Ligation artefacts can leave several adapters in tandem at a read end, so
    /// trimming continues as long as another adapter is anchored right at the
    /// current cut position.
    pub fn find_trim(&mut self, read: &[u8]) -> Option<AdapterTrim> {
        if read.is_empty() {
            return None;
        }

        let read_len = read.len();

        // Collect the hits first so that the (mutable) searcher borrow ends before
        // we inspect them.
        let hits: Vec<Match> = self
            .searcher
            .search_all_encoded_patterns(&self.encoded, read, ADAPTER_MAX_EDITS)
            .to_vec();

        let mut trim = AdapterTrim::nothing(read_len);
        if let Some((_, keep_start)) = self.resolve_side(&hits, read_len, AdapterSide::Start) {
            trim.keep_start = keep_start;
            trim.sides[0] = true;
        }
        if let Some((_, keep_end)) = self.resolve_side(&hits, read_len, AdapterSide::End) {
            trim.keep_end = keep_end;
            trim.sides[1] = true;
        }

        if trim.is_empty() {
            return None;
        }

        // Never return an empty interval; leaving the read as-is is better than
        // dropping it.
        if trim.keep_end <= trim.keep_start
            || trim.keep_end - trim.keep_start < ADAPTER_MIN_INSERT_LEN
        {
            return None;
        }
        Some(trim)
    }

    /// Walk inwards from one end of the read over consecutive adapters and return
    /// `(score, position_to_cut_at)`.
    ///
    /// The first adapter has to reach the very end of the read; each following one
    /// has to be anchored at the position of the previous cut (tandem adapters).
    fn resolve_side(
        &self,
        hits: &[Match],
        read_len: usize,
        side: AdapterSide,
    ) -> Option<(f64, usize)> {
        // Start at the very beginning/end of the read and move inwards.
        let mut pos = match side {
            AdapterSide::Start => 0,
            AdapterSide::End => read_len,
        };
        let mut result: Option<(f64, usize)> = None;

        // Each round has to advance `pos`, so this terminates; the bound is only
        // there to keep a pathological read from looping for a long time.
        for _ in 0..MAX_ADAPTERS_PER_END {
            let mut best: Option<(f64, usize)> = None;

            for m in hits {
                let Some((matches, identity)) =
                    score_match(m, self.min_match_len, self.min_identity)
                else {
                    continue;
                };
                // Weight matching bases by how clean the alignment is, so that an
                // alignment consisting mostly of insertions does not outrank a
                // clean (possibly shorter) one.
                let score = matches as f64 * identity;

                let (anchored, new_pos) = match side {
                    AdapterSide::Start => (
                        // Starts at/after the current cut (a few bases of junk in
                        // front of the adapter are allowed) ...
                        m.text_start <= pos + self.end_slack
                            // ... but anchored at it, not a re-alignment of the
                            // adapter we already walked past.
                            && m.text_start + ADAPTER_CHAIN_OVERLAP >= pos
                            // and it has to make progress.
                            && m.text_end > pos,
                        m.text_end,
                    ),
                    AdapterSide::End => (
                        m.text_end + self.end_slack >= pos
                            && m.text_end <= pos + ADAPTER_CHAIN_OVERLAP
                            && m.text_start < pos,
                        m.text_start,
                    ),
                };
                if !anchored {
                    continue;
                }

                let better = match best {
                    // On equal score walk over the outermost adapter first.
                    Some((best_score, best_pos)) => match side {
                        AdapterSide::Start => {
                            score > best_score || (score == best_score && new_pos > best_pos)
                        }
                        AdapterSide::End => {
                            score > best_score || (score == best_score && new_pos < best_pos)
                        }
                    },
                    None => true,
                };
                if better {
                    best = Some((score, new_pos));
                }
            }

            match best {
                Some((score, new_pos)) => {
                    pos = new_pos;
                    result = Some((score, new_pos));
                }
                None => break,
            }
        }

        result
    }

    /// Trim the adapters off a sequence (and quality string, if given).
    ///
    /// Returns `None` when there was nothing to trim.
    pub fn trim<'a>(
        &mut self,
        seq: &'a [u8],
        qual: Option<&'a [u8]>,
    ) -> Option<(&'a [u8], Option<&'a [u8]>)> {
        let trim = self.find_trim(seq)?;
        let seq = &seq[trim.keep_start..trim.keep_end];
        let qual = qual.map(|q| &q[trim.keep_start..trim.keep_end.min(q.len())]);
        Some((seq, qual))
    }
}

/// Number of matching bases and the identity of the alignment.
///
/// The identity is `matches / alignment_columns`, i.e. gaps count as mismatches.
/// This mirrors Porechop's "aligned region percent identity" and prevents
/// alignments full of insertions/deletions from being mistaken for an adapter.
fn score_match(m: &Match, min_match_len: usize, min_identity: f64) -> Option<(usize, f64)> {
    let mut matches = 0usize;
    let mut columns = 0usize;
    for op in &m.cigar.ops {
        let cnt = op.cnt as usize;
        columns += cnt;
        if op.op == pa_types::CigarOp::Match {
            matches += cnt;
        }
    }

    if columns == 0 || matches < min_match_len {
        return None;
    }
    let identity = matches as f64 / columns as f64;
    if identity < min_identity {
        return None;
    }
    Some((matches, identity))
}

/// Statistics of an adapter trimming run.
#[derive(Debug, Default, Clone, Copy)]
pub struct AdapterTrimStats {
    pub total: usize,
    pub trimmed_start: usize,
    pub trimmed_end: usize,
    pub trimmed_both: usize,
}

/// Default number of threads used for adapter trimming.
pub const DEFAULT_ADAPTER_THREADS: usize = 10;

/// Number of records that are read and trimmed as one parallel batch.
///
/// Small enough to keep the batch (which is copied out of the parser) bounded in
/// memory even for very long reads, large enough to keep the threads busy.
const ADAPTER_BATCH_RECORDS: usize = 512;

/// Size of the read-ahead buffer used in front of a decoder.
///
/// Uncompressed input is read straight from a `BufReader`; for gzip input this
/// decouples the parser from the decoder so that the decoder's background
/// workers are not woken up for (and stalled by) the parser's small reads.
const ADAPTER_READ_BUFFER: usize = 1 << 20;

/// One input file to trim, together with the file it is written to.
struct AdapterJob {
    input: PathBuf,
    output: PathBuf,
    gzip: bool,
}

/// A single read while it is being trimmed.
///
/// Holding the record data (instead of a reference into the parser's batch)
/// allows the reads of one batch to be trimmed by several threads at once, while
/// they are still written out in input order afterwards.
struct PendingRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    /// Which read ends were trimmed, for the statistics.
    sides: [bool; 2],
}

impl PendingRecord {
    fn new(id: &[u8], seq: &[u8], qual: Option<&[u8]>) -> Self {
        Self {
            id: id.to_vec(),
            seq: seq.to_vec(),
            qual: match qual {
                Some(qual) => qual.to_vec(),
                None => vec![b'!'; seq.len()],
            },
            sides: [false, false],
        }
    }

    /// Trim the ligation adapter off both read ends (if present).
    fn trim(&mut self, trimmer: &mut AdapterTrimmer) {
        let Some(trim) = trimmer.find_trim(&self.seq) else {
            return;
        };
        if trim.is_empty() {
            return;
        }

        self.seq = self.seq[trim.keep_start..trim.keep_end].to_vec();
        self.qual = self.qual[trim.keep_start..trim.keep_end].to_vec();
        self.sides = trim.sides;
    }
}

/// Standalone adapter trimming of FASTQ file(s).
///
/// Reads are written to `output` (a file, or a folder when `output` has no
/// extension) and keep their original order and orientation, also when
/// `config.threads` is greater than one: the reads of a batch are trimmed in
/// parallel but are always written out in input order afterwards.
pub fn trim_adapters_in_fastq(
    input: &[PathBuf],
    output: &str,
    config: &AdapterTrimConfig,
) -> anyhow::Result<AdapterTrimStats> {
    validate_fastq_paths(input)?;

    // Plan the jobs (and create the output folders) up front: this keeps the
    // trimming itself free of filesystem bookkeeping, and any I/O errors are
    // reported before any work is done.
    let jobs = plan_jobs(input, output, config)?;

    let threads = config.threads.max(1);
    // The searcher is not `Sync` and keeps internal buffers, so every worker gets
    // its own instance (see `AdapterTrimmer`).
    let pool = if threads > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .map_err(|err| {
                    anyhow!("Failed to create a thread pool of {threads} threads: {err}")
                })?,
        )
    } else {
        None
    };

    let mut stats = AdapterTrimStats::default();
    match pool.as_ref() {
        Some(pool) => {
            // Several input files (e.g. one per barcode) are trimmed at the same
            // time; each file is written by exactly one thread, so a file's reads
            // stay in input order.
            let results: Vec<anyhow::Result<AdapterTrimStats>> = pool.install(|| {
                jobs.par_iter()
                    .map(|job| trim_adapter_file(job, threads, Some(pool), config))
                    .collect()
            });
            for result in results {
                stats.merge(result?);
            }
        }
        None => {
            for job in &jobs {
                stats.merge(trim_adapter_file(job, threads, None, config)?);
            }
        }
    }

    Ok(stats)
}

/// Work out the input -> output file mapping and make sure the output folders
/// exist. A single input file can be written straight to `output`; multiple files
/// need one output per input, so they go into a folder.
fn plan_jobs(
    input: &[PathBuf],
    output: &str,
    config: &AdapterTrimConfig,
) -> anyhow::Result<Vec<AdapterJob>> {
    let single_file_output = input.len() == 1 && Path::new(output).extension().is_some();

    let mut jobs = Vec::with_capacity(input.len());
    for read_path in input {
        let out_path = if single_file_output {
            PathBuf::from(output)
        } else {
            let stem = read_path
                .file_name()
                .map(|name| name.to_string_lossy().to_string())
                .unwrap_or_else(|| "reads.fastq".to_string());
            Path::new(output).join(format!("{stem}.adapter_trimmed.fastq"))
        };
        if let Some(parent) = out_path.parent()
            && !parent.as_os_str().is_empty()
        {
            std::fs::create_dir_all(parent)?;
        }

        let gzip = config.gzip || out_path.extension().map(|ext| ext == "gz").unwrap_or(false);
        jobs.push(AdapterJob {
            input: read_path.clone(),
            output: out_path,
            gzip,
        });
    }

    Ok(jobs)
}

/// Open a FASTQ file as a batch reader, decoding gzip input in parallel.
///
/// Regular files are the fast path:
///
/// * gzip is decoded by `rapidgzip`, whose marker/window based decoder splits the
///   stream into independent blocks that are inflated by its own workers.
///   Decoding a multi-hundred-MB FASTQ serially would otherwise pin a single core
///   before the trimming threads ever see a byte;
/// * uncompressed input is wrapped in a plain `BufReader`.
///
/// Non-regular inputs (a FIFO, `/dev/stdin`, a process substitution) cannot be
/// decoded in parallel because every parallel path needs positional reads, so they
/// keep going through paraseq's `niffler`-based reader, which sniffs the format on
/// a single handle.
fn open_batch_reader(
    path: &Path,
    batch_records: usize,
    threads: usize,
) -> anyhow::Result<paraseq::fastq::Reader<Box<dyn Read + Send>>> {
    if !std::fs::metadata(path).is_ok_and(|meta| meta.is_file()) {
        return Ok(paraseq::fastq::Reader::from_path_with_batch_size(
            path,
            batch_records,
        )?);
    }

    let source: Box<dyn Read + Send> = if is_gzip_file(path) {
        Box::new(
            rapidgzip_core::Decoder::builder()
                .decoder_threads(threads)
                .build()?
                .open(path)?,
        )
    } else {
        Box::new(BufReader::with_capacity(
            ADAPTER_READ_BUFFER,
            File::open(path)?,
        ))
    };
    Ok(paraseq::fastq::Reader::with_batch_size(
        source,
        batch_records,
    )?)
}

/// Whether the regular file `path` holds gzip data, by magic number or extension.
///
/// The magic number is probed on a throwaway handle so that the parser's own
/// stream is never consumed. This is why the caller has to keep the function away
/// from non-regular files: opening a second handle on a pipe would steal bytes
/// from it.
fn is_gzip_file(path: &Path) -> bool {
    let mut magic = [0u8; 2];
    if File::open(path)
        .and_then(|mut file| file.read(&mut magic))
        .is_ok_and(|n| n == 2)
    {
        return magic == [0x1f, 0x8b];
    }
    path.extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

/// Trim the adapters of one FASTQ file, writing the result to `job.output`.
fn trim_adapter_file(
    job: &AdapterJob,
    threads: usize,
    pool: Option<&rayon::ThreadPool>,
    config: &AdapterTrimConfig,
) -> anyhow::Result<AdapterTrimStats> {
    let mut stats = AdapterTrimStats::default();

    let mut reader = open_batch_reader(&job.input, ADAPTER_BATCH_RECORDS, threads)
        .map_err(|err| anyhow!("Failed to open FASTQ file '{}': {err}", job.input.display()))?;
    let mut record_set = reader.new_record_set();

    let file = File::create(&job.output).map_err(|err| {
        anyhow!(
            "Failed to create output file '{}': {err}",
            job.output.display()
        )
    })?;
    let mut writer: Box<dyn Write> = if job.gzip {
        Box::new(flate2::write::GzEncoder::new(
            file,
            flate2::Compression::default(),
        ))
    } else {
        Box::new(BufWriter::new(file))
    };

    while record_set
        .fill(&mut reader)
        .map_err(|err| anyhow!("Error reading FASTQ file '{}': {err}", job.input.display()))?
    {
        let mut batch: Vec<PendingRecord> = record_set
            .iter()
            .map(|record| {
                let record = record.map_err(|err| {
                    anyhow!(
                        "Error reading FASTQ record in '{}': {err}",
                        job.input.display()
                    )
                })?;
                Ok(PendingRecord::new(
                    record.id(),
                    record.seq().as_ref(),
                    record.qual(),
                ))
            })
            .collect::<anyhow::Result<_>>()?;

        // Trim the batch: in parallel when more than one thread was requested,
        // otherwise in a plain loop that can reuse a single searcher.
        match pool {
            Some(pool) if batch.len() > 1 => {
                // One chunk per thread; rayon hands the chunks to its workers.
                let chunk_size = batch.len().div_ceil(threads).max(1);
                pool.install(|| {
                    batch
                        .par_chunks_mut(chunk_size)
                        // `for_each_init` gives every rayon worker its own
                        // searcher, which is much cheaper than building one per
                        // read (`AdapterTrimmer::new` encodes the adapter).
                        .for_each_init(
                            || AdapterTrimmer::new(config.clone()),
                            |trimmer, chunk| {
                                for record in chunk {
                                    record.trim(trimmer);
                                }
                            },
                        );
                });
            }
            _ => {
                let mut trimmer = AdapterTrimmer::new(config.clone());
                for record in &mut batch {
                    record.trim(&mut trimmer);
                }
            }
        }

        for record in &batch {
            stats.total += 1;
            stats.count(&record.sides);

            writer.write_all(b"@")?;
            writer.write_all(&record.id)?;
            writer.write_all(b"\n")?;
            writer.write_all(&record.seq)?;
            writer.write_all(b"\n+\n")?;
            writer.write_all(&record.qual)?;
            writer.write_all(b"\n")?;
        }
    }

    writer.flush()?;

    Ok(stats)
}

impl AdapterTrimStats {
    fn count(&mut self, sides: &[bool; 2]) {
        match (sides[0], sides[1]) {
            (true, true) => self.trimmed_both += 1,
            (true, false) => self.trimmed_start += 1,
            (false, true) => self.trimmed_end += 1,
            (false, false) => {}
        }
    }

    fn merge(&mut self, other: Self) {
        self.total += other.total;
        self.trimmed_start += other.trimmed_start;
        self.trimmed_end += other.trimmed_end;
        self.trimmed_both += other.trimmed_both;
    }

    pub fn trimmed(&self) -> usize {
        self.trimmed_start + self.trimmed_end + self.trimmed_both
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&c| RC[c as usize]).collect()
}

const RC: [u8; 256] = {
    let mut rc = [0; 256];
    let mut i = 0;
    while i < 256 {
        rc[i] = i as u8;
        i += 1;
    }
    rc[b'A' as usize] = b'T';
    rc[b'C' as usize] = b'G';
    rc[b'T' as usize] = b'A';
    rc[b'G' as usize] = b'C';
    rc[b'a' as usize] = b't';
    rc[b'c' as usize] = b'g';
    rc[b't' as usize] = b'a';
    rc[b'g' as usize] = b'c';
    rc[b'N' as usize] = b'N';
    rc[b'n' as usize] = b'n';
    rc
};

#[cfg(test)]
mod tests {
    use super::*;

    fn trimmer() -> AdapterTrimmer {
        AdapterTrimmer::new(AdapterTrimConfig::default())
    }

    #[test]
    fn test_full_forward_adapter_at_start() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [ADAPTER_SEQUENCE, insert.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, ADAPTER_SEQUENCE.len());
        assert_eq!(trim.keep_end, read.len());
        assert!(trim.sides[0] && !trim.sides[1]);
    }

    #[test]
    fn test_full_reverse_complement_adapter_at_start() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let rc_adapter = reverse_complement(ADAPTER_SEQUENCE);
        let read = [rc_adapter.as_slice(), insert.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, ADAPTER_SEQUENCE.len());
        assert!(trim.sides[0] && !trim.sides[1]);
    }

    #[test]
    fn test_full_forward_adapter_at_end() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [insert.as_slice(), ADAPTER_SEQUENCE].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, 0);
        assert_eq!(trim.keep_end, insert.len());
        assert!(!trim.sides[0] && trim.sides[1]);
    }

    #[test]
    fn test_full_reverse_complement_adapter_at_end() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let rc_adapter = reverse_complement(ADAPTER_SEQUENCE);
        let read = [insert.as_slice(), rc_adapter.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_end, insert.len());
        assert!(!trim.sides[0] && trim.sides[1]);
    }

    #[test]
    fn test_adapters_at_both_ends() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [
            ADAPTER_SEQUENCE,
            insert.as_slice(),
            &reverse_complement(ADAPTER_SEQUENCE),
        ]
        .concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, ADAPTER_SEQUENCE.len());
        assert_eq!(trim.keep_end, ADAPTER_SEQUENCE.len() + insert.len());
        assert!(trim.sides[0] && trim.sides[1]);
    }

    #[test]
    fn test_partial_adapter_at_start() {
        let mut t = trimmer();
        // Only the last 20 bases of the adapter are present (basecalling can drop
        // the first bases of a read).
        let partial = &ADAPTER_SEQUENCE[7..];
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [partial, insert.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, partial.len());
        assert!(trim.sides[0]);
    }

    #[test]
    fn test_adapter_with_mismatches_is_trimmed() {
        let mut t = trimmer();
        let mut adapter = ADAPTER_SEQUENCE.to_vec();
        adapter[3] = b'A';
        adapter[15] = b'C';
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [adapter.as_slice(), insert.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, adapter.len());
    }

    #[test]
    fn test_no_adapter_is_untouched() {
        let mut t = trimmer();
        let read = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        assert!(t.find_trim(&read).is_none());
        assert!(t.trim(&read, None).is_none());
    }

    #[test]
    fn test_middle_adapter_is_not_trimmed() {
        let mut t = trimmer();
        let flank = vec![b'A'; 300];
        let read = [flank.as_slice(), ADAPTER_SEQUENCE, flank.as_slice()].concat();
        assert!(t.find_trim(&read).is_none());
    }

    #[test]
    fn test_trim_returns_slices() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [ADAPTER_SEQUENCE, insert.as_slice()].concat();
        let qual = vec![b'I'; read.len()];
        let (seq, qual) = t.trim(&read, Some(&qual)).unwrap();
        assert_eq!(seq, insert);
        assert_eq!(qual.unwrap().len(), insert.len());
    }

    #[test]
    fn test_trim_does_not_empty_read() {
        let mut t = trimmer();
        // A read that is (almost) nothing but adapter should not be emptied
        // entirely by an over-eager hit.
        let read = ADAPTER_SEQUENCE.to_vec();
        if let Some(trim) = t.find_trim(&read) {
            assert!(trim.keep_end > trim.keep_start);
        }
    }

    #[test]
    fn test_adapters_survive_internal_gaps() {
        let mut t = trimmer();
        // Adapter with a deletion and a mismatch in the middle.
        let mut adapter = ADAPTER_SEQUENCE.to_vec();
        adapter.remove(12);
        adapter[5] = b'A';
        let insert = b"ACGTACGTACGTACGTACGTACGTACGT";
        let read = [adapter.as_slice(), insert.as_slice()].concat();
        let trim = t.find_trim(&read).expect("should still trim");
        assert!(trim.keep_start >= ADAPTER_SEQUENCE.len() - 3);
    }

    #[test]
    fn test_terminal_junk_before_adapter_is_removed() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        // A few spurious bases in front of the adapter (common in ONT reads).
        let read = [b"GGGG".as_slice(), ADAPTER_SEQUENCE, insert.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, 4 + ADAPTER_SEQUENCE.len());
    }

    #[test]
    fn test_two_adapters_in_tandem_at_start() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        // Ligation artefacts can leave two adapters in tandem.
        let read = [ADAPTER_SEQUENCE, ADAPTER_SEQUENCE, insert.as_slice()].concat();
        let trim = t.find_trim(&read).unwrap();
        assert_eq!(trim.keep_start, 2 * ADAPTER_SEQUENCE.len());
    }

    #[test]
    fn test_single_read_with_adapter_at_end_is_not_start_trimmed() {
        let mut t = trimmer();
        let insert = b"ACGTACGTACGTACGTACGT";
        let read = [insert.as_slice(), ADAPTER_SEQUENCE].concat();
        let trim = t.find_trim(&read).unwrap();
        assert!(!trim.sides[0], "must not trim the start of this read");
    }

    #[test]
    fn test_real_lsk114_read() {
        // First 493bp of the read used in dorado's adapter trimming tests
        // (tests/data/adapter_trim/lsk110_single_read.fastq). The 5' end is a
        // partial AMX114 adapter.
        let read = b"TACTTCGTTCCAGTTACGTATTGCTTGCAGGTCCAGACGGTGATGGTGATCGGCGGCGTGCTCATGATCGCAGGCGTGCTCCTTGGTCCGGTTCGGTTGATCGATGAAAAGCAGGAATCGCAAAAAAACAAAAACCCCGCCCGAAGGCGGGGTTTTCCCGAACAGCACAGTCGCAGGAATTAACAACATCCCCCCTATCAGGAGAGACGGCCCGCCGGCGAGCTGATCGAGCGGCCTGTGACCTGATCTG";
        let mut t = trimmer();
        let trim = t.find_trim(read).expect("adapter should be found");
        assert_eq!(trim.keep_start, 24);
        assert_eq!(trim.keep_end, read.len());
        assert!(trim.sides[0] && !trim.sides[1]);
    }

    /// Assemble a FASTQ file with `n` reads; every third read has an adapter at
    /// the start, every fifth read a reverse complement adapter at the end.
    fn write_test_fastq(path: &Path, n: usize) {
        let mut out = Vec::new();
        for i in 0..n {
            let mut seq = Vec::new();
            if i % 3 == 0 {
                seq.extend_from_slice(ADAPTER_SEQUENCE);
            }
            // Only ACGT(N), which is what the searcher accepts.
            for j in 0..36 {
                seq.push(b"ACGTN"[(i + j) % 5]);
            }
            if i % 5 == 0 {
                seq.extend_from_slice(&reverse_complement(ADAPTER_SEQUENCE));
            }
            let qual = vec![b'I'; seq.len()];
            out.extend_from_slice(format!("@read{i}\n").as_bytes());
            out.extend_from_slice(&seq);
            out.extend_from_slice(b"\n+\n");
            out.extend_from_slice(&qual);
            out.extend_from_slice(b"\n");
        }
        std::fs::write(path, out).unwrap();
    }

    #[test]
    fn test_trim_adapters_in_fastq_multithreaded_matches_single_threaded() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("reads.fastq");
        write_test_fastq(&input, 1500);

        let single_out = dir.path().join("single.fastq");
        let multi_out = dir.path().join("multi.fastq");

        let single_stats = trim_adapters_in_fastq(
            std::slice::from_ref(&input),
            single_out.to_str().unwrap(),
            &AdapterTrimConfig {
                threads: 1,
                ..Default::default()
            },
        )
        .unwrap();
        let multi_stats = trim_adapters_in_fastq(
            std::slice::from_ref(&input),
            multi_out.to_str().unwrap(),
            &AdapterTrimConfig {
                threads: 8,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(
            std::fs::read(&single_out).unwrap(),
            std::fs::read(&multi_out).unwrap(),
            "read order and content must not depend on the number of threads"
        );
        assert_eq!(single_stats.total, multi_stats.total);
        assert_eq!(single_stats.trimmed(), multi_stats.trimmed());
        assert_eq!(single_stats.trimmed_start, multi_stats.trimmed_start);
        assert_eq!(single_stats.trimmed_end, multi_stats.trimmed_end);
        assert_eq!(single_stats.trimmed_both, multi_stats.trimmed_both);
        assert!(single_stats.trimmed() > 0);
    }

    #[test]
    fn test_trim_adapters_in_fastq_multiple_files_in_folder() {
        let dir = tempfile::tempdir().unwrap();
        let files: Vec<PathBuf> = ["a.fastq", "b.fastq", "c.fastq"]
            .iter()
            .map(|name| {
                let path = dir.path().join(name);
                write_test_fastq(&path, 200);
                path
            })
            .collect();

        let single_dir = dir.path().join("single");
        let multi_dir = dir.path().join("multi");

        let single_stats = trim_adapters_in_fastq(
            &files,
            single_dir.to_str().unwrap(),
            &AdapterTrimConfig {
                threads: 1,
                ..Default::default()
            },
        )
        .unwrap();
        let multi_stats = trim_adapters_in_fastq(
            &files,
            multi_dir.to_str().unwrap(),
            &AdapterTrimConfig {
                threads: 4,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(single_stats.total, 3 * 200);
        assert_eq!(multi_stats.total, single_stats.total);
        assert_eq!(multi_stats.trimmed(), single_stats.trimmed());
        for input in &files {
            let name = format!(
                "{}.adapter_trimmed.fastq",
                input.file_name().unwrap().to_string_lossy()
            );
            assert_eq!(
                std::fs::read(single_dir.join(&name)).unwrap(),
                std::fs::read(multi_dir.join(&name)).unwrap(),
                "{name} differs between 1 and 4 threads"
            );
        }
    }

    /// Gzip input is decoded by a different reader than plain input (rapidgzip
    /// instead of `niffler` via paraseq), so it has to produce the same reads and
    /// the same trim decisions.
    #[test]
    fn test_trim_adapters_in_fastq_gzip_input_matches_plain_input() {
        let dir = tempfile::tempdir().unwrap();
        let plain = dir.path().join("reads.fastq");
        let gz = dir.path().join("reads.fastq.gz");
        write_test_fastq(&plain, 500);
        {
            let file = File::create(&gz).unwrap();
            let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            std::io::copy(&mut File::open(&plain).unwrap(), &mut encoder).unwrap();
            encoder.finish().unwrap();
        }

        let plain_out = dir.path().join("plain_out.fastq");
        let gz_out = dir.path().join("gz_out.fastq");

        let plain_stats = trim_adapters_in_fastq(
            std::slice::from_ref(&plain),
            plain_out.to_str().unwrap(),
            &AdapterTrimConfig {
                threads: 1,
                ..Default::default()
            },
        )
        .unwrap();
        let gz_stats = trim_adapters_in_fastq(
            std::slice::from_ref(&gz),
            gz_out.to_str().unwrap(),
            &AdapterTrimConfig {
                threads: 4,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(gz_stats.total, plain_stats.total);
        assert_eq!(gz_stats.trimmed(), plain_stats.trimmed());
        assert_eq!(gz_stats.trimmed_start, plain_stats.trimmed_start);
        assert_eq!(gz_stats.trimmed_end, plain_stats.trimmed_end);
        assert_eq!(gz_stats.trimmed_both, plain_stats.trimmed_both);
        assert_eq!(
            std::fs::read(&plain_out).unwrap(),
            std::fs::read(&gz_out).unwrap(),
            "gzip input must trim to the same bytes as plain input"
        );
    }

    /// A `.fastq` file that happens to hold gzip data (what `cat reads.fastq.gz
    /// > reads.fastq` produces) is detected by magic number, not by name.
    #[test]
    fn test_gzip_input_without_gz_extension_is_detected() {
        let dir = tempfile::tempdir().unwrap();
        let plain = dir.path().join("reads.fastq");
        let misnamed = dir.path().join("reads_no_suffix.fastq");
        write_test_fastq(&plain, 50);
        {
            let file = File::create(&misnamed).unwrap();
            let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            std::io::copy(&mut File::open(&plain).unwrap(), &mut encoder).unwrap();
            encoder.finish().unwrap();
        }

        assert!(is_gzip_file(&misnamed));
        assert!(!is_gzip_file(&plain));

        let out = dir.path().join("out.fastq");
        let stats = trim_adapters_in_fastq(
            std::slice::from_ref(&misnamed),
            out.to_str().unwrap(),
            &AdapterTrimConfig {
                threads: 2,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(stats.total, 50);
        assert!(stats.trimmed() > 0);
    }
}
