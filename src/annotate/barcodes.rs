use colored::Colorize;
use needletail::{Sequence, parse_fastx_file};
use sassy::profiles::{Iupac, Profile};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, PartialOrd, Ord, Eq, Serialize, Deserialize)]
pub enum BarcodeType {
    Ftag,
    Rtag,
    Fflank, // Not public, in case Fbar is not detected
    Rflank, // Not public, in case Rbar is not detected
}

impl BarcodeType {
    pub fn as_flank(&self) -> Self {
        match self {
            BarcodeType::Ftag => BarcodeType::Fflank,
            BarcodeType::Rtag => BarcodeType::Rflank,
            _ => panic!("Cannot convert {self:?} to flank"),
        }
    }

    pub fn as_str(&self) -> &str {
        match self {
            BarcodeType::Ftag => "Ftag",
            BarcodeType::Rtag => "Rtag",
            BarcodeType::Fflank => "Fflank",
            BarcodeType::Rflank => "Rflank",
        }
    }
}

#[derive(Clone, Debug, PartialEq, PartialOrd, Ord, Eq, Serialize, Deserialize)]
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

#[derive(Clone, Debug, PartialEq, PartialOrd, Ord, Eq, Serialize, Deserialize)]
pub struct BarcodeGroup {
    pub flank: Vec<u8>,
    pub flank_prefix: Vec<u8>,
    pub flank_suffix: Vec<u8>,
    pub bar_region: (usize, usize),
    pub barcodes: Vec<Barcode>, // Think always assume IUPAC anyway
    pub k_cutoff: Option<usize>,
    pub barcode_type: BarcodeType,
    pub prefix_k: Option<usize>,
    pub suffix_k: Option<usize>,
}

impl BarcodeGroup {
    pub fn new(
        query_seqs: Vec<Vec<u8>>,
        query_labels: Vec<String>,
        barcode_type: BarcodeType,
    ) -> Self {
        // Often a group shares shared flanks, we assume they do for now, and extract
        // their prefix and suffixes, if present
        if query_seqs.len() == 1 {
            panic!(
                "For now we only support 'groups' if you only have one query, just add a 'random' second query to your fasta file with the same left and right flank, just chnge the barcode"
            )
        }
        let query_seqs_refs: Vec<&[u8]> = query_seqs.iter().map(|seq| seq.as_slice()).collect();
        let (prefix, suffix) = Self::get_flanks(&query_seqs_refs);
        let mut flank = Vec::new();

        let prefix_len = prefix.as_ref().unwrap_or(&Vec::new()).len();
        let suffix_len = suffix.as_ref().unwrap_or(&Vec::new()).len();
        let mask_size = query_seqs[0].len() - prefix_len - suffix_len;

        // Add prefix if present
        if let Some(p) = &prefix {
            flank.extend_from_slice(p);
        }
        if mask_size > 0 {
            flank.extend(vec![b'N'; mask_size]);
        }
        if let Some(s) = &suffix {
            flank.extend_from_slice(s);
        }

        // Slice out all the masked region sequences
        let mut barcodes = Vec::new();
        for (i, seq) in query_seqs.iter().enumerate() {
            //let start = prefix_len.saturating_sub(PADDING);
            //let end = (start + mask_size + PADDING).min(seq.len());
            barcodes.push(Barcode::new(
                seq, // trying the whole sequence
                &query_labels[i],
                barcode_type.clone(),
            ));
        }

        Self {
            flank,
            flank_prefix: prefix.clone().unwrap_or_default(),
            flank_suffix: suffix.clone().unwrap_or_default(),
            bar_region: (prefix_len, prefix_len + mask_size - 1), //inclsuve
            barcodes,
            k_cutoff: None,
            barcode_type,
            prefix_k: None,
            suffix_k: None,
        }
    }

    pub fn display(&self) {
        let (mask_start, mask_end) = self.bar_region;
        let left_flank = &self.flank[..mask_start];
        let right_flank = &self.flank[mask_end + 1..]; // As inclusive mask end
        // Create colored string to show flank composition
        // left flank & right_flank = blue, mask = cyan
        let left_flank_str = String::from_utf8_lossy(left_flank).blue();
        let right_flank_str = String::from_utf8_lossy(right_flank).blue();
        let mask_size = mask_end - mask_start + 1;
        let mask_str = "-".repeat(mask_size).to_string().bright_yellow();

        println!("{left_flank_str}{mask_str}{right_flank_str}");
        for barcode in self.barcodes.iter().take(5) {
            println!(
                "{}: {}",
                barcode.label.green(),
                String::from_utf8_lossy(&barcode.seq).bright_yellow()
            );
        }
        if self.barcodes.len() > 2 {
            println!("...+{} more", self.barcodes.len() - 2);
        }
    }

    pub fn new_from_fasta(fasta_file: &str, bar_type: BarcodeType) -> Self {
        let mut reader = parse_fastx_file(fasta_file).expect("Query file not found");
        let mut bar_seqs: Vec<Vec<u8>> = Vec::new();
        let mut labels = Vec::new();

        // Collect all records first
        while let Some(record) = reader.next() {
            let seqrec = record.expect("invalid record");
            let norm_seq = seqrec.normalize(false); // Makes all uppercase, important for set_perc threshold
            labels.push(std::str::from_utf8(seqrec.id()).unwrap().to_string());
            bar_seqs.push(norm_seq.into_owned());
        }
        Self::new(bar_seqs, labels, bar_type)
    }

    pub fn set_perc_threshold(&mut self, perc: f32) {
        // The flank has an N mask so we should not count these in the perc as
        // they always match. Also the user can supply N so we just count N instead then
        let n_count: usize = self.flank.iter().filter(|&c| *c == b'N').count();
        let informative_flank = self.flank.len() - n_count;
        self.k_cutoff = Some((informative_flank as f32 * perc).ceil() as usize);
        // Also set k for each barcode
        for barcode in self.barcodes.iter_mut() {
            barcode.k_cutoff = Some((barcode.seq.len() as f32 * perc) as usize);
        }
        // Do the same for the suffixes
        self.prefix_k = Some((self.flank_prefix.len() as f32 * perc).ceil() as usize);
        self.suffix_k = Some((self.flank_suffix.len() as f32 * perc).ceil() as usize);
    }

    pub fn set_flank_threshold(&mut self, flank_threshold: usize) {
        self.k_cutoff = Some(flank_threshold);
    }

    pub fn set_barcode_threshold(&mut self, barcode_threshold: usize) {
        for barcode in self.barcodes.iter_mut() {
            barcode.k_cutoff = Some(barcode_threshold);
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
    // use crate::annotate::tune::tune_single_sequence;

    #[test]
    fn test_display() {
        let seqs = [
            b"AAATTTGGG".as_slice(),
            b"AAACCCGGG".as_slice(),
            b"AAATATGGG".as_slice(),
        ];
        let labels = vec!["s1".to_string(), "s2".to_string(), "s3".to_string()];
        let barcode_type = BarcodeType::Ftag;
        let barcode_group = BarcodeGroup::new(
            vec![seqs[0].to_vec(), seqs[1].to_vec(), seqs[2].to_vec()],
            labels,
            barcode_type,
        );
        barcode_group.display();
    }

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
        let seqs = [b"AAATTTGGG".as_slice(), b"AAACCCGGG".as_slice()];
        let labels = vec!["s1".to_string(), "s2".to_string()];
        let barcode_type = BarcodeType::Ftag;
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
        let seqs = [b"@@@@@@@@@".as_slice(), b"AAACCCGGG".as_slice()];
        let labels = vec!["s1".to_string(), "s2".to_string()];
        let barcode_type = BarcodeType::Ftag;
        let _ = BarcodeGroup::new(
            vec![seqs[0].to_vec(), seqs[1].to_vec()],
            labels,
            barcode_type,
        );
    }

    #[test]
    #[should_panic]
    fn test_barcode_group_unequal_length() {
        let seqs = [b"AAATTTGGG".as_slice(), b"AAAAAAACCCGGG".as_slice()];
        let labels = vec!["s1".to_string(), "s2".to_string()];
        let barcode_type = BarcodeType::Ftag;
        let _ = BarcodeGroup::new(
            vec![seqs[0].to_vec(), seqs[1].to_vec()],
            labels,
            barcode_type,
        );
    }

    #[test]
    fn test_fasta_read() {
        let example_file = "examples/rapid_bars.fasta";
        let barcode_group = BarcodeGroup::new_from_fasta(example_file, BarcodeType::Ftag);
        let expected_flank = b"GCTTGGGTGTTTAACCNNNNNNNNNNNNNNNNNNNNNNNNGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";
        let flank = barcode_group.flank;
        assert_eq!(flank, expected_flank);
        assert_eq!(barcode_group.bar_region, (16, 39));
        assert_eq!(&flank[16..=39], b"NNNNNNNNNNNNNNNNNNNNNNNN");
        assert_eq!(barcode_group.barcodes.len(), 96);
        assert_eq!(
            barcode_group.barcodes[0].seq,
            b"AAGAAAGTTGTCGGTGTCTTTGTG".to_vec() // NB01 fwd
        );
    }
}
