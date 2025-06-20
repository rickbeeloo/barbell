use crate::search::distribution::{TailSide, get_fp_threshold};
use needletail::{FastxReader, Sequence, parse_fastx_file};
use rand::Rng;
use rayon::prelude::*;
use sassy::profiles::{Iupac, Profile};
use sassy::search::Searcher;

#[derive(Clone, Debug, PartialEq, PartialOrd, Ord, Eq)]
pub enum BarcodeType {
    Fbar,
    Rbar,
    FFlank, // Not public, in case Fbar is not detected
    RFlank, // Not public, in case Rbar is not detected
}

impl BarcodeType {
    pub fn as_flank(&self) -> Self {
        match self {
            BarcodeType::Fbar => BarcodeType::FFlank,
            BarcodeType::Rbar => BarcodeType::RFlank,
            _ => panic!("Cannot convert {:?} to flank", self),
        }
    }
}

pub struct Barcode {
    pub seq: Vec<u8>,
    pub label: String,
    pub match_type: BarcodeType,
    pub k_cutoff: Option<usize>, // Should be 'tuned'
}

impl Barcode {
    pub fn new(seq: &[u8], label: &str, match_type: BarcodeType) -> Self {
        if !Iupac::valid_seq(seq) {
            panic!("Sequence contains character not supported by IUPAC");
        }
        Self {
            seq: seq.to_vec(),
            label: label.to_string(),
            match_type,
            k_cutoff: None,
        }
    }
}

/// Tune multiple sequences in parallel and return their k-cutoffs
pub fn tune_sequences_parallel(
    sequences: &[&[u8]],
    n_iter: usize,
    fp_target: f32,
    alpha: f32,
    n_threads: usize,
) -> Vec<usize> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .unwrap();

    pool.install(|| {
        sequences
            .par_iter()
            .map(|seq| tune_single_sequence(seq, n_iter, fp_target, alpha))
            .collect()
    })
}

fn tune_single_sequence(seq: &[u8], n_iter: usize, fp_target: f32, alpha: f32) -> usize {
    let min_len: usize = 0;
    let min_len = seq.len() * 2;
    let max_len = min_len + 1;
    // let max_len = seq.len() * 2;
    println!("Seq len: {}", seq.len());
    let mut costs = Vec::new();
    let mut searcher = Searcher::<Iupac>::new_rc_with_overhang(alpha);

    for _ in 0..n_iter {
        let random_seq = random_dna_seq(min_len, max_len);
        let matches = searcher.search(seq, &random_seq, seq.len());

        if matches.is_empty() {
            continue;
        }
        let lowest_cost = matches.iter().map(|m| m.cost).min().unwrap();
        costs.push(lowest_cost);
    }
    let cutoff = get_fp_threshold(costs, fp_target, TailSide::Left);
    println!("Cutoff: {}", cutoff);
    cutoff as usize
}

pub fn random_dna_seq(min_len: usize, max_len: usize) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let len = rng.gen_range(min_len..max_len);
    let mut seq = Vec::new();
    for _ in 0..len {
        let base = match rng.gen_range(0..4) {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };
        seq.push(base);
    }
    seq
}

pub struct BarcodeGroup {
    pub flank: Vec<u8>,
    pub bar_region: (usize, usize),
    pub barcodes: Vec<Barcode>,  // Think always assume IUPAC anyway
    pub k_cutoff: Option<usize>, // Should be 'tuned'
}

impl BarcodeGroup {
    pub fn new(
        query_seqs: Vec<Vec<u8>>,
        query_labels: Vec<String>,
        barcode_type: BarcodeType,
    ) -> Self {
        // Often a group shares shared flanks, we assume they do for now, and extract
        // their prefix and suffixes, if present
        let query_seqs_refs: Vec<&[u8]> = query_seqs.iter().map(|seq| seq.as_slice()).collect();
        let (prefix, suffix) = Self::get_flanks(&query_seqs_refs);
        let mut flank = Vec::new();

        let prefix_len = prefix.as_ref().unwrap_or(&Vec::new()).len();
        let suffix_len = suffix.as_ref().unwrap_or(&Vec::new()).len();
        let mask_size = query_seqs[0].len() - prefix_len - suffix_len;

        // Add prefix if present
        if let Some(p) = prefix {
            flank.extend_from_slice(&p);
        }
        if mask_size > 0 {
            flank.extend(vec![b'N'; mask_size]);
        }
        if let Some(s) = suffix {
            flank.extend_from_slice(&s);
        }

        // Slice out all the masked region sequences
        let mut barcodes = Vec::new();
        for (i, seq) in query_seqs.iter().enumerate() {
            let start = prefix_len;
            let end = start + mask_size;
            println!(
                "Barcode region: {:?}",
                String::from_utf8_lossy(&seq[start..end])
            );
            barcodes.push(Barcode::new(
                &seq[start..end],
                &query_labels[i],
                barcode_type.clone(),
            ));
        }

        Self {
            flank,
            bar_region: (prefix_len, prefix_len + mask_size - 1),
            barcodes,
            k_cutoff: None,
        }
    }

    pub fn new_from_fasta(fasta_file: &str) -> Self {
        let mut reader = parse_fastx_file(fasta_file).expect("valid path/file");
        let mut bar_seqs = Vec::new();
        let mut labels = Vec::new();

        // Collect all records first
        while let Some(record) = reader.next() {
            let seqrec = record.expect("invalid record");
            let norm_seq = seqrec.normalize(false);
            labels.push(std::str::from_utf8(seqrec.id()).unwrap().to_string());
            bar_seqs.push(norm_seq.into_owned());
        }
        Self::new(bar_seqs, labels, BarcodeType::Fbar)
    }

    /// Tune the group (flank) and all barcodes in parallel
    pub fn tune_group(&mut self, n_iter: usize, fp_target: f32, alpha: f32, n_threads: usize) {
        // Collect all sequences to tune: flank + all barcodes
        let mut sequences = vec![self.flank.as_slice()];
        for barcode in &self.barcodes {
            sequences.push(barcode.seq.as_slice());
        }

        // Tune all sequences in parallel
        let cutoffs = tune_sequences_parallel(&sequences, n_iter, fp_target, alpha, n_threads);

        // Set the cutoffs
        self.k_cutoff = Some(cutoffs[0]); // First cutoff is for the flank
        for (i, barcode) in self.barcodes.iter_mut().enumerate() {
            barcode.k_cutoff = Some(cutoffs[i + 1]); // Skip first (flank), rest are barcodes
        }
    }

    pub fn get_flanks(seqs: &[&[u8]]) -> (Option<Vec<u8>>, Option<Vec<u8>>) {
        // This does not make sense if the sequences are not equally long
        assert!(
            !seqs.iter().any(|s| s.len() != seqs[0].len()),
            "All sequences per group must be equally long"
        );

        let prefix = Self::longest_common_prefix(seqs);
        let suffix = Self::longest_common_suffix(seqs);

        (prefix, suffix)
    }

    pub fn longest_common_prefix(seqs: &[&[u8]]) -> Option<Vec<u8>> {
        if seqs.is_empty() {
            return None;
        }

        let first = seqs[0];
        let mut common_len = first.len();

        for seq in &seqs[1..] {
            common_len = common_len.min(
                first
                    .iter()
                    .zip(seq.iter())
                    .take_while(|(a, b)| a == b)
                    .count(),
            );
            if common_len == 0 {
                return None;
            }
        }

        Some(first[..common_len].to_vec())
    }

    pub fn longest_common_suffix(seqs: &[&[u8]]) -> Option<Vec<u8>> {
        if seqs.is_empty() {
            return None;
        }

        let first = seqs[0];
        let mut common_len = first.len();

        for seq in &seqs[1..] {
            common_len = common_len.min(
                first
                    .iter()
                    .rev()
                    .zip(seq.iter().rev())
                    .take_while(|(a, b)| a == b)
                    .count(),
            );
            if common_len == 0 {
                return None;
            }
        }

        Some(first[first.len() - common_len..].to_vec())
    }
}

mod tests {
    use super::*;

    #[test]
    fn test_longest_common_prefix() {
        let seqs = vec![
            b"ACGTAGAGAG".as_slice(),
            b"ACGTAGACTA".as_slice(),
            b"ACGAGCAGGA".as_slice(),
        ];
        let prefix = BarcodeGroup::longest_common_prefix(&seqs);
        assert_eq!(prefix, Some(b"ACG".to_vec()));
    }

    #[test]
    fn test_longest_common_suffix() {
        let seqs = vec![
            b"ACGTAGAGAGGGA".as_slice(),
            b"ACGTTAGACTAGA".as_slice(),
            b"ACGAGCAGGAGAA".as_slice(),
        ];
        let suffix = BarcodeGroup::longest_common_suffix(&seqs);
        assert_eq!(suffix, Some(b"A".to_vec()));
    }

    #[test]
    fn test_barcode_group() {
        let seqs = vec![b"AAATTTGGG".as_slice(), b"AAACCCGGG".as_slice()];
        let labels = vec!["s1".to_string(), "s2".to_string()];
        let barcode_type = BarcodeType::Fbar;
        let barcode_group = BarcodeGroup::new(
            vec![seqs[0].to_vec(), seqs[1].to_vec()],
            labels,
            barcode_type,
        );
        assert_eq!(barcode_group.flank, b"AAANNNGGG".to_vec());
        assert_eq!(barcode_group.bar_region, (3, 5));
        assert_eq!(barcode_group.barcodes.len(), 2);
        assert_eq!(barcode_group.barcodes[0].seq, b"TTT".to_vec());
        assert_eq!(barcode_group.barcodes[1].seq, b"CCC".to_vec());
    }

    #[test]
    #[should_panic]
    fn test_barcode_group_invalid_seq() {
        let seqs = vec![b"@@@@@@@@@".as_slice(), b"AAACCCGGG".as_slice()];
        let labels = vec!["s1".to_string(), "s2".to_string()];
        let barcode_type = BarcodeType::Fbar;
        let _ = BarcodeGroup::new(
            vec![seqs[0].to_vec(), seqs[1].to_vec()],
            labels,
            barcode_type,
        );
    }

    #[test]
    #[should_panic]
    fn test_barcode_group_unequal_length() {
        let seqs = vec![b"AAATTTGGG".as_slice(), b"AAAAAAACCCGGG".as_slice()];
        let labels = vec!["s1".to_string(), "s2".to_string()];
        let barcode_type = BarcodeType::Fbar;
        let _ = BarcodeGroup::new(
            vec![seqs[0].to_vec(), seqs[1].to_vec()],
            labels,
            barcode_type,
        );
    }

    #[test]
    fn test_barcode_tune_k() {
        let barcode = b"CACAAAGACACCGACAACTTTCTT";
        let cutoff = tune_single_sequence(barcode, 10000, 0.01, 0.5);
        // We always expect k to be less than half of the length of the barcode
        // as that is what we get for random sequences on average
        println!("Got cut off: {}", cutoff);
        assert!(cutoff < barcode.len() / 2);
    }

    #[test]
    fn test_barcode_group_tune_k() {
        let mut barcode_group = BarcodeGroup::new(
            vec![
                b"CACAAAGACACAAAAAAAAAAAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
                    .to_vec(),
                b"CACAAAGACACTTTTTTTTTTTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
                    .to_vec(),
            ],
            vec!["s1".to_string(), "s2".to_string()],
            BarcodeType::Fbar,
        );
        barcode_group.tune_group(1000, 0.01, 0.5, 4);
        // We always expect k to be less than half of the length of the barocde
        // as that is what we get for random sequences on average
        println!(
            "Got cut off: {}",
            barcode_group.k_cutoff.unwrap_or(usize::MAX)
        );
        assert!(barcode_group.k_cutoff.unwrap_or(usize::MAX) < barcode_group.flank.len() / 2);

        // Also check that barcodes were tuned
        for barcode in &barcode_group.barcodes {
            assert!(barcode.k_cutoff.is_some());
            println!(
                "Barcode {} cutoff: {}",
                barcode.label,
                barcode.k_cutoff.unwrap()
            );
        }
    }

    #[test]
    fn test_fasta_read() {
        let example_file = "examples/rapid_bars.fasta";
        let barcode_group = BarcodeGroup::new_from_fasta(example_file);
        let expected_flank = b"TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCTTGGGTGTTTAACCNNNNNNNNNNNNNNNNNNNNNNNNGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";
        let flank = barcode_group.flank;
        assert_eq!(flank, expected_flank);
        assert_eq!(barcode_group.bar_region, (52, 75));
        assert_eq!(barcode_group.barcodes.len(), 97);
        assert_eq!(
            barcode_group.barcodes[0].seq,
            b"AAGAAAGTTGTCGGTGTCTTTGTG".to_vec() // NB01 fwd
        );
    }
}
