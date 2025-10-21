use crate::simulations::mutate::*;
use needletail::{FastxReader, Sequence, parse_fastx_file};
use rand::Rng;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

const RAPID_LEFT: &str = "GCTTGGGTGTTTAACC";
const RAPID_RIGHT: &str = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";

const READ_MIN_LEN: usize = 600;
const READ_MAX_LEN: usize = 4000;

const MAX_TRIM: usize = 20;
const MIN_DOULBE_SPACE: usize = 10;

const BARCODE_PATH: &str = "data/rapid_bars.txt";
// const RC_FRACTION: f64 = 0.5;

const MAX_EDITS: usize = 6;

// Assuming barcode is 5' - barcode - 3'
fn get_adapter_region(barcode: &[u8]) -> Vec<u8> {
    let mut seq = Vec::new();
    seq.extend_from_slice(RAPID_LEFT.as_bytes());
    seq.extend_from_slice(barcode);
    seq.extend_from_slice(RAPID_RIGHT.as_bytes());
    seq
}

fn create_random_sequence(length: usize) -> Vec<u8> {
    // Random dna sequence
    let base_seq = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();
    let seq = (0..length)
        .map(|_| base_seq[rng.random_range(0..4)])
        .collect();
    seq
}

fn create_random_quality_score(length: usize) -> Vec<u8> {
    let mut rng = rand::rng();
    (0..length)
        .map(|_| {
            let phred = rng.random_range(0..40); // 0..=39
            (phred + 33) as u8 // convert to ASCII
        })
        .collect()
}

fn read_barcodes(barcode_fasta: &str) -> Vec<(String, Vec<u8>)> {
    let mut barcodes: Vec<(String, Vec<u8>)> = vec![];
    let mut reader = parse_fastx_file(barcode_fasta).expect("valid path/file");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let norm_seq = seqrec.normalize(false);
        let seq_id = String::from_utf8_lossy(seqrec.id())
            .to_string()
            .trim()
            .to_string();
        barcodes.push((seq_id, norm_seq.into_owned()));
    }
    barcodes
}

struct SeqElement {
    seq: Vec<u8>,
    pos: usize,
}

struct MockRead {
    pub id: String,
    elements: Vec<SeqElement>,
    seq: Vec<u8>,
}

impl MockRead {
    pub fn new(id: String, read_len: usize) -> Self {
        Self {
            id,
            elements: vec![],
            seq: create_random_sequence(read_len),
        }
    }

    pub fn new_rand_len(id: String) -> Self {
        let read_len = rand::rng().random_range(READ_MIN_LEN..READ_MAX_LEN);
        Self::new(id, read_len)
    }

    pub fn add_element(&mut self, seq: &[u8], pos: usize) {
        self.elements.push(SeqElement {
            seq: seq.to_vec(),
            pos,
        });
    }

    pub fn render(&mut self) -> Vec<u8> {
        for element in &self.elements {
            self.seq
                .splice(element.pos..element.pos, element.seq.to_vec());
        }
        self.seq.clone()
    }
}

struct ReadCollection {
    pub reads: Vec<(MockRead, Option<String>)>,
    fastq_writer: BufWriter<File>,
    truth_writer: BufWriter<File>,
}

impl ReadCollection {
    pub fn new(fastq_out: &str, truth_out: &str) -> Self {
        let fastq_writer = File::create(fastq_out).expect("failed to create file");
        let truth_writer = File::create(truth_out).expect("failed to create file");
        let fastq_writer = BufWriter::new(fastq_writer);
        let truth_writer = BufWriter::new(truth_writer);
        Self {
            reads: vec![],
            fastq_writer,
            truth_writer,
        }
    }

    pub fn add_read(&mut self, read: MockRead, truth: Option<String>) {
        self.reads.push((read, truth));
    }

    pub fn dump(&mut self, rc_frac: f64) {
        for (read, truth) in &mut self.reads {
            let mut seq = read.render();

            // If we have a RC_FRACTION chance, we reverse complement the sequence
            if rand::rng().random_range(0.0..1.0) < rc_frac {
                seq = reverse_complement(&seq);
            }

            // Also add edits, random amount between 0 and MAX_EDITS
            let seq = mutate_sequence(&seq, 0, MAX_EDITS);

            let quality_score = create_random_quality_score(seq.len());
            let record = format!(
                "@{}\n{}\n+\n{}\n",
                read.id,
                String::from_utf8_lossy(&seq),
                String::from_utf8_lossy(&quality_score)
            );
            self.fastq_writer
                .write_all(record.as_bytes())
                .expect("failed to write record");

            if let Some(truth) = truth {
                let line = format!("{}\t{truth}\n", read.id);
                self.truth_writer
                    .write_all(line.as_bytes())
                    .expect("failed to write truth record");
            }
        }
    }
}

fn group1(fasta_out: &str, truth_out: &str, n: usize, barcode_path: &str, rc_frac: f64) {
    let mut read_collection = ReadCollection::new(fasta_out, truth_out);
    for i in 0..n {
        let read = MockRead::new_rand_len(format!("seq_{i}"));
        read_collection.add_read(read, None);
    }
    read_collection.dump(rc_frac);
}

fn group2(fasta_out: &str, truth_out: &str, n: usize, barcode_path: &str, rc_frac: f64) {
    // Read all barcodes for rapdi barcoding kit
    let barcodes = read_barcodes(barcode_path);
    let n_barcodes = barcodes.len();
    assert!(n_barcodes > 0);

    // Create the read collection
    let mut read_collection = ReadCollection::new(fasta_out, truth_out);

    // Add n reandom reads to read collection, where we have <L><B><R> regions
    for i in 0..n {
        let mut read = MockRead::new_rand_len(format!("seq_{i}"));
        // Select random barcode
        let barcode_idx = rand::rng().random_range(0..n_barcodes);
        let (barcode_name, barcode_seq) = &barcodes[barcode_idx];
        let adapter_region = get_adapter_region(barcode_seq);
        // Add to read as PREFIX (index = 0)
        read.add_element(&adapter_region, 0);
        // We add it to the collection with "barcode name" as truth as this
        // is what the demux tools should detect
        read_collection.add_read(read, Some(barcode_name.clone()));
    }

    // Dump the read collection
    read_collection.dump(rc_frac);
}

///Group III = Group II + random trim left end
fn group3(fasta_out: &str, truth_out: &str, n: usize, barcode_path: &str, rc_frac: f64) {
    // Read all barcodes for rapdi barcoding kit
    let barcodes = read_barcodes(barcode_path);
    let n_barcodes = barcodes.len();
    assert!(n_barcodes > 0);

    // Create the read collection
    let mut read_collection = ReadCollection::new(fasta_out, truth_out);

    // Add n reandom reads to read collection, where we have <L><B><R> regions
    for i in 0..n {
        let mut read = MockRead::new_rand_len(format!("seq_{i}"));
        // Select random barcode
        let barcode_idx = rand::rng().random_range(0..n_barcodes);
        let (barcode_name, barcode_seq) = &barcodes[barcode_idx];
        let adapter_region = get_adapter_region(barcode_seq);
        // Do random trim on the adapter region ONLY on the left side in this case as we prefix it so
        // experimentally only expect trimming on "overhang"
        let trimmed_adapter_region = random_trim_side(&adapter_region, MAX_TRIM, true, false);
        assert!(trimmed_adapter_region.len() <= adapter_region.len());
        // Add to read as PREFIX (index = 0)
        read.add_element(&trimmed_adapter_region, 0);
        // We add it to the collection with "barcode name" as truth as this
        // is what the demux tools should detect
        read_collection.add_read(read, Some(barcode_name.clone()));
    }

    // Dump the read collection
    read_collection.dump(rc_frac);
}

fn group4(fasta_out: &str, truth_out: &str, n: usize, barcode_path: &str, rc_frac: f64) {
    // Read all barcodes for rapdi barcoding kit
    let barcodes = read_barcodes(barcode_path);
    let n_barcodes = barcodes.len();
    assert!(n_barcodes > 0);

    // Create the read collection
    let mut read_collection = ReadCollection::new(fasta_out, truth_out);

    // Add n reandom reads to read collection, where we have <L><B><R> regions
    for i in 0..n {
        let mut read = MockRead::new_rand_len(format!("seq_{i}"));

        // Create first barcode region
        let first_barcode_idx = rand::rng().random_range(0..n_barcodes);
        let (first_barcode_name, first_barcode_seq) = &barcodes[first_barcode_idx];
        let first_adapter_region = get_adapter_region(first_barcode_seq);

        // Second barcode region - where we enforce DIFFERENT barcode
        let barcodes_left = (0..n_barcodes)
            .filter(|idx| *idx != first_barcode_idx)
            .collect::<Vec<_>>();

        let second_barcode_idx = barcodes_left[rand::rng().random_range(0..barcodes_left.len())];

        let (second_barcode_name, second_barcode_seq) = &barcodes[second_barcode_idx];

        assert_ne!(first_barcode_name, second_barcode_name); // Sanity check they are different
        let second_adapter_region = get_adapter_region(second_barcode_seq);

        // We combine these two seperated by the MIN_DOULBE_SPACE
        let random_seq = create_random_sequence(MIN_DOULBE_SPACE);
        let mut combined_seq = Vec::new();
        combined_seq.extend_from_slice(&first_adapter_region);
        combined_seq.extend_from_slice(&random_seq);
        combined_seq.extend_from_slice(&second_adapter_region);
        assert_eq!(
            combined_seq.len(),
            first_adapter_region.len() + MIN_DOULBE_SPACE + second_adapter_region.len()
        );

        // Then we insert this as prefix to the read
        read.add_element(&combined_seq, 0);

        // As truth lets actullay use the combined labels
        let truth_label = format!("{first_barcode_name}_{second_barcode_name}_double_front");
        read_collection.add_read(read, Some(truth_label));
    }

    // Dump the read collection
    read_collection.dump(rc_frac);
}

fn group5(fasta_out: &str, truth_out: &str, n: usize, barcode_path: &str, rc_frac: f64) {
    // Read all barcodes for rapdi barcoding kit
    let barcodes = read_barcodes(barcode_path);
    let n_barcodes = barcodes.len();
    assert!(n_barcodes > 0);

    // Create the read collection
    let mut read_collection = ReadCollection::new(fasta_out, truth_out);

    // Add n reandom reads to read collection, where we have <L><B><R> regions
    for i in 0..n {
        let mut read = MockRead::new_rand_len(format!("seq_{i}"));

        // Create first barcode region
        let first_barcode_idx = rand::rng().random_range(0..n_barcodes);
        let (first_barcode_name, first_barcode_seq) = &barcodes[first_barcode_idx];
        let first_adapter_region = get_adapter_region(first_barcode_seq);

        // We insert this as prefix
        read.add_element(&first_adapter_region, 0);

        // Second barcode region - where we enforce DIFFERENT barcode
        let barcodes_left = (0..n_barcodes)
            .filter(|idx| *idx != first_barcode_idx)
            .collect::<Vec<usize>>();

        // Random index in the range of left barcode indices
        let second_barcode_idx = rand::rng().random_range(0..barcodes_left.len());

        let (second_barcode_name, second_barcode_seq) =
            &barcodes[barcodes_left[second_barcode_idx]];

        assert_ne!(first_barcode_name, second_barcode_name); // Sanity check they are different
        let second_adapter_region = get_adapter_region(second_barcode_seq);

        // We insert this in the middle of the read
        let read_len = read.seq.len();
        let middle_index = (read_len / 2) - first_adapter_region.len();

        read.add_element(&second_adapter_region, middle_index);

        // As truth lets actullay use the combined labels
        let truth_label = format!("{first_barcode_name}_{second_barcode_name}_mid_insert");
        read_collection.add_read(read, Some(truth_label));
    }

    // Dump the read collection
    read_collection.dump(rc_frac);
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut rev_comp = Vec::new();
    for base in seq.iter().rev() {
        rev_comp.push(match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'A',
        });
    }
    rev_comp
}

fn group6(fasta_out: &str, truth_out: &str, n: usize, barcode_path: &str, rc_frac: f64) {
    // Read all barcodes for rapdi barcoding kit
    let barcodes = read_barcodes(barcode_path);
    let n_barcodes = barcodes.len();
    assert!(n_barcodes > 0);

    // Create the read collection
    let mut read_collection = ReadCollection::new(fasta_out, truth_out);

    // Add n reandom reads to read collection, where we have <L><B><R> regions
    for i in 0..n {
        let mut read = MockRead::new_rand_len(format!("seq_{i}"));

        // Create first barcode region
        let first_barcode_idx = rand::rng().random_range(0..n_barcodes);
        let (first_barcode_name, first_barcode_seq) = &barcodes[first_barcode_idx];
        let first_adapter_region = get_adapter_region(first_barcode_seq);

        // We insert this as prefix
        read.add_element(&first_adapter_region, 0);

        // Second barcode region - where we enforce DIFFERENT barcode
        let barcodes_left = (0..n_barcodes)
            .filter(|idx| *idx != first_barcode_idx)
            .collect::<Vec<_>>();
        let second_barcode_idx = rand::rng().random_range(0..barcodes_left.len());
        let (second_barcode_name, second_barcode_seq) = &barcodes[second_barcode_idx];
        assert_ne!(first_barcode_name, second_barcode_name); // Sanity check they are different
        let second_adapter_region = get_adapter_region(second_barcode_seq);

        // We insert this at the end, as suffix, we also reverse complement it to match
        // expected of dual end
        let rev_comp_second_adapter = reverse_complement(&second_adapter_region);
        println!(
            "Forward: {}",
            String::from_utf8_lossy(&second_adapter_region)
        );
        println!(
            "Backward: {}",
            String::from_utf8_lossy(&rev_comp_second_adapter)
        );
        read.add_element(
            &rev_comp_second_adapter,
            read.seq.len() + first_adapter_region.len(),
        );

        // As truth lets actullay use the combined labels
        let truth_label = format!("{first_barcode_name}_{second_barcode_name}_front_back");
        read_collection.add_read(read, Some(truth_label));
    }

    // Dump the read collection
    read_collection.dump(rc_frac);
}

#[derive(Debug)]
pub enum TestGroup {
    GroupI,   // Random sequence, no barcode OR adapter at all
    GroupII,  // Valid barcode + adapter
    GroupIII, // Group II + random trim left end
    GroupIV,  // INVALID: GroupII * 2, double barcode
    GroupV,   // INVALID: Barcode at the frond + one at the middle
    GroupVI,  // INVALID:Barocde at both ends
}

// Implement string for enum
impl std::fmt::Display for TestGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

pub fn create_testdata(n: usize, sim_out_dir: &str, barcode_file: &str, rc_frac: f64) {
    // If sim_out_dir does not exist, create it
    if !Path::new(sim_out_dir).exists() {
        println!("Creating directory: {sim_out_dir}");
        std::fs::create_dir_all(sim_out_dir).expect("failed to create directory");
    }

    let all_groups = [
        TestGroup::GroupI,
        TestGroup::GroupII,
        TestGroup::GroupIII,
        TestGroup::GroupIV,
        TestGroup::GroupV,
    ];

    let group_fns = vec![group1, group2, group3, group4, group5, group6];
    for (group, group_fn) in all_groups.iter().zip(group_fns) {
        let fasta_out = format!("{sim_out_dir}/{group}.fastq");
        let truth_out = format!("{sim_out_dir}/{group}_truth.txt");
        group_fn(
            fasta_out.as_str(),
            truth_out.as_str(),
            n,
            barcode_file,
            rc_frac,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_create_random_sequence() {
    //     let seq = create_random_sequence(10);
    //     assert_eq!(seq.len(), 10);
    //     println!("Sequence: {}", String::from_utf8_lossy(&seq));
    // }

    // #[test]
    // fn test_get_adatper() {
    //     let barcode = b"xxxx";
    //     let region = get_adapter_region(barcode);
    //     let expected = b"GCTTGGGTGTTTAACCxxxxGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";
    //     assert_eq!(region, expected);
    // }

    // #[test]
    // fn test_read_barcodes() {
    //     let barcode_file = "data/rapid_bars.txt";
    //     let barcodes = read_barcodes(barcode_file);
    //     assert!(barcodes.len() > 1);
    //     let first_barcode = barcodes.first().unwrap();
    //     assert_eq!(first_barcode.0, "BC01");
    //     assert_eq!(first_barcode.1, "TACATGCTCCTGTTGTTAGGGAGG".as_bytes());
    // }

    // #[test]
    // fn test_mock_read() {
    //     let mut read = MockRead::new("seq_0".to_string(), 100);
    //     read.add_element(b"AAAA", 0);
    //     read.add_element(b"CCCC", 20);
    //     read.add_element(b"GGGG", 40);
    //     let seq = read.render();
    //     assert_eq!(&seq[0..4], "AAAA".as_bytes());
    //     assert_eq!(&seq[20..24], "CCCC".as_bytes());
    //     assert_eq!(&seq[40..44], "GGGG".as_bytes());
    //     println!("Mock read: {}", String::from_utf8_lossy(&seq));
    // }

    // #[test]
    // fn test_simulation_data() {
    //     let n = 100;
    //     let sim_out_dir = "simulated_data";
    //     let barcode_file = "data/rapid_bars.txt";
    //     create_testdata(n, sim_out_dir, barcode_file, 0.0);
    // }
}
