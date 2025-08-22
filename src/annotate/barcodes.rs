use crate::kits::kits::*;
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
    pub pad_region: (usize, usize), // barcode + padding
    pub barcodes: Vec<Barcode>,
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
        let (pad_start, pad_end) = (prefix_len.saturating_sub(10), prefix_len + mask_size + 10);
        for (i, seq) in query_seqs.iter().enumerate() {
            let start = pad_start;
            let end: usize = pad_end.min(seq.len());

            //assert_eq!(*&seq[start..end].len(), 24 + 20);
            barcodes.push(Barcode::new(
                &seq[start..end], // trying the whole sequence
                &query_labels[i],
                barcode_type.clone(),
            ));
        }

        Self {
            flank,
            flank_prefix: prefix.clone().unwrap_or_default(),
            flank_suffix: suffix.clone().unwrap_or_default(),
            bar_region: (prefix_len, prefix_len + mask_size - 1), //inclsuve
            pad_region: (pad_start, pad_end),
            barcodes,
            k_cutoff: None,
            barcode_type,
            prefix_k: None,
            suffix_k: None,
        }
    }

    pub fn display(&self, n: usize) {
        // This will show with padding, if we get the bar region
        let (mask_start, mask_end) = self.bar_region;

        // Create colored string to show flank composition
        // left flank & right_flank = blue, mask = cyan
        let left_flank_str = String::from_utf8_lossy(self.flank_prefix.as_slice()).blue();
        let right_flank_str = String::from_utf8_lossy(self.flank_suffix.as_slice()).blue();
        let mask_size = mask_end - mask_start + 1;
        let mask_str = "-".repeat(mask_size).to_string().bright_yellow();

        println!("{left_flank_str}{mask_str}{right_flank_str}");

        let left_flank_len = self.flank_prefix.len();
        for barcode in self.barcodes.iter().take(n) {
            // We have "padding" which we remove before printing
            let (pad_start, _pad_end) = self.pad_region;
            let (bar_start, bar_end) = self.bar_region;
            let len = barcode.seq.len();
            let start_pos = bar_start.saturating_sub(pad_start).min(len);
            //fixme: better logic here although in practice it should never happen
            let mut end_pos = (bar_end + 1).saturating_sub(pad_start).min(len);
            if end_pos < start_pos {
                end_pos = start_pos;
            }

            let just_bar_slice = &barcode.seq[start_pos..end_pos];

            // Align the barcode sequence so it starts under the mask region
            let label_text = format!("{}: ", barcode.label);
            let visible_label_len = label_text.len();
            let pad_spaces = left_flank_len.saturating_sub(visible_label_len);
            let pad_str = if pad_spaces == 0 {
                " ".to_string()
            } else {
                " ".repeat(pad_spaces)
            };

            println!(
                "{}{}{}",
                label_text.green(),
                pad_str,
                String::from_utf8_lossy(just_bar_slice).bright_yellow()
            );
        }
        if self.barcodes.len() > 2 {
            println!("...+{} more", self.barcodes.len() - 2);
        }
    }

    pub fn new_from_kit(kit: &str) -> Vec<Self> {
        // Get all kit info
        let kit_config = get_kit_info(kit);

        // We have to build all the queries, either one or two groups
        let mut groups = Vec::new();

        // Check if we have left flank
        if let Some(left_prefix) = kit_config.left_prefix {
            let mut query_seqs = Vec::new();
            let mut query_labels = Vec::new();
            let prefix = left_prefix.to_string();
            let suffix = kit_config.left_suffix.unwrap_or_default().to_string();
            let barcode_list = get_barcodes(kit_config.barcodes.from, kit_config.barcodes.to);
            for barcode_name in barcode_list {
                let barcode_seq = lookup_barcode_seq(&barcode_name)
                    .expect("Barcode not found - odd - raise issue");
                let left_flank_query = format!("{prefix}{barcode_seq}{suffix}");
                query_seqs.push(left_flank_query.as_bytes().to_vec());
                query_labels.push(barcode_name);
            }
            groups.push(BarcodeGroup::new(
                query_seqs,
                query_labels,
                BarcodeType::Ftag,
            ));
        }

        // Check if we have right flank
        if let Some(right_prefix) = kit_config.right_prefix {
            let mut query_seqs = Vec::new();
            let mut query_labels = Vec::new();
            let prefix = right_prefix.to_string();
            let suffix = kit_config.right_suffix.unwrap_or_default().to_string();
            let barcode_list = get_barcodes(
                kit_config.barcodes2.unwrap().from,
                kit_config.barcodes2.unwrap().to,
            );
            for barcode_name in barcode_list {
                let barcode_seq = lookup_barcode_seq(&barcode_name)
                    .expect("Barcode not found - odd - raise issue");
                let right_flank_query = format!("{prefix}{barcode_seq}{suffix}");
                query_seqs.push(right_flank_query.as_bytes().to_vec());
                query_labels.push(barcode_name);
            }
            groups.push(BarcodeGroup::new(
                query_seqs,
                query_labels,
                BarcodeType::Rtag,
            ));
        }

        // Collect all records first
        groups
    }

    pub fn new_from_fasta(fasta_file: &str, bar_type: BarcodeType) -> Self {
        let mut reader = parse_fastx_file(fasta_file).expect("Query file not found");
        let mut bar_seqs: Vec<Vec<u8>> = Vec::new();
        let mut labels = Vec::new();

        // Collect all records first
        while let Some(record) = reader.next() {
            let seqrec = record.expect("invalid record");
            let norm_seq = seqrec.normalize(true); // Makes all uppercase, important for set_perc threshold
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

    pub fn get_effective_len(&self) -> usize {
        // Simply length of prefix + suffix
        self.flank_prefix.len() + self.flank_suffix.len()
    }
}

#[cfg(test)]
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
        barcode_group.display(5);
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
        // We add a 10bp padding on each side of the mask, but since we don't have 10 chars
        // here it just maxes out and takes full flanks
        assert_eq!(barcode_group.barcodes[0].seq, b"AAATTTGGG".to_vec());
        assert_eq!(barcode_group.barcodes[1].seq, b"AAACCCGGG".to_vec());
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
            &barcode_group.barcodes[0].seq[10..10 + 24],
            b"AAGAAAGTTGTCGGTGTCTTTGTG".to_vec() // NB01 fwd
        );
    }

    #[test]
    fn new_from_kit_rapid_barcodes() {
        let groups = BarcodeGroup::new_from_kit("SQK-NBD114-96");
        for group in groups {
            group.display(10);
        }
    }
}
