use pa_types::{Cost, Pos};
use sassy::profiles::{Iupac, Profile};
use sassy::search::*;

use crate::search::barcodes::{BarcodeGroup, BarcodeType, Tuneable};
use crate::search::interval::collapse_overlapping_matches;

pub struct Demuxer {
    overhang_searcher: Searcher<Iupac>,
    searcher: Searcher<Iupac>,
    queries: Vec<BarcodeGroup>,
}

#[derive(Debug, Clone)]
pub struct BarbellMatch {
    pub read_start_bar: usize,
    pub read_end_bar: usize,
    pub read_start_flank: usize,
    pub read_end_flank: usize,
    pub bar_start: usize,
    pub bar_end: usize,
    pub match_type: BarcodeType,
    pub flank_cost: Cost,
    pub barcode_cost: Cost,
    pub label: String,
    pub strand: Strand,
}

impl BarbellMatch {
    pub fn new(
        read_start_bar: usize,
        read_end_bar: usize,
        read_start_flank: usize,
        read_end_flank: usize,
        bar_start: usize,
        bar_end: usize,
        match_type: BarcodeType,
        flank_cost: Cost,
        barcode_cost: Cost,
        label: String,
        strand: Strand,
    ) -> Self {
        Self {
            read_start_bar,
            read_end_bar,
            read_start_flank,
            read_end_flank,
            bar_start,
            bar_end,
            match_type,
            flank_cost,
            barcode_cost,
            label,
            strand,
        }
    }
}

impl Demuxer {
    pub fn new(alpha: f32) -> Self {
        // Use rc searcher with alpha, where alpha is overhang cost
        let overhang_searcher = Searcher::<Iupac>::new_rc_with_overhang(alpha);
        let searcher = Searcher::<Iupac>::new_rc();
        Self {
            overhang_searcher,
            searcher,
            queries: Vec::new(),
        }
    }

    pub fn add_query_group(&mut self, query_group: BarcodeGroup) -> &mut Self {
        self.queries.push(query_group);
        self
    }

    fn slice_masked_region<'a>(
        &self,
        read: &'a [u8],
        flank: &BarcodeGroup,
        flank_match: &Match,
    ) -> (&'a [u8], (usize, usize)) {
        let (mask_start, mask_end) = flank.bar_region;
        //  (&read[mask_start..=mask_end], (mask_start, mask_end))
        let path = flank_match.to_path();

        // Collect unique positions and find min/max r_pos for positions in mask region
        let positions: std::collections::HashSet<_> = path.iter().collect();
        let (min_r_pos, max_r_pos) = positions
            .iter()
            .filter(|pos| pos.0 as usize >= mask_start && pos.0 as usize <= mask_end)
            .fold((usize::MAX, 0), |(min, max), pos| {
                let r_pos = pos.1 as usize;
                (min.min(r_pos), max.max(r_pos))
            });

        if min_r_pos == usize::MAX {
            return (&[], (0, 0));
        }

        let start = min_r_pos.saturating_sub(crate::WIGGLE_ROOM);
        let end = (max_r_pos + crate::WIGGLE_ROOM).min(read.len());

        println!(
            "Mask region: [{}-{}] c: {}{:?}",
            start,
            end,
            flank_match.cost,
            String::from_utf8_lossy(&read[start..=end])
        );
        (&read[start..=end], (start, end))
    }

    pub fn demux(&mut self, read: &[u8]) -> Vec<BarbellMatch> {
        let mut results: Vec<BarbellMatch> = Vec::new();

        // We first search for the top level of the barcodeGroups
        for barcode_group in self.queries.iter() {
            let flank = &barcode_group.flank;

            let flank_matches =
                self.overhang_searcher
                    .search(flank, &read, barcode_group.k_cutoff.unwrap_or(0));

            // If we found a flank, we slice out the masked region and search for the barcodes in the
            // masked region, we should be AWARE OF THE STRAND as sassy now returns both directions (Fwd, Rc)
            for flank_match in flank_matches {
                println!(
                    "Flank slice: {:?}",
                    String::from_utf8_lossy(
                        &read[flank_match.start.1 as usize..flank_match.end.1 as usize],
                    )
                );
                let (mask_region, (mask_start, _mask_end)) =
                    self.slice_masked_region(read, barcode_group, &flank_match);

                // If no mask match found we can just skip to the next one
                if mask_region.is_empty() {
                    continue;
                }

                // Then we compare each query against the masked region
                for barcode in barcode_group.barcodes.iter() {
                    // Look for barcode mathces in mask slice
                    let bms = self.overhang_searcher.search(
                        &barcode.seq,
                        &mask_region,
                        barcode.k_cutoff.unwrap_or(0),
                    );

                    // If no barcode matches within k we can just skip
                    if bms.is_empty() {
                        continue;
                    }

                    // We make sure they flank and barcode match on the same strand
                    for bm in bms.iter().filter(|bm| bm.strand == flank_match.strand) {
                        results.push(BarbellMatch::new(
                            //t - storing flank matches to more accurately filter out overlaps later on
                            bm.start.1 as usize + mask_start,
                            bm.end.1 as usize + mask_start,
                            flank_match.start.1 as usize,
                            flank_match.end.1 as usize,
                            //q
                            bm.start.0 as usize,
                            bm.end.0 as usize,
                            barcode.match_type.clone(),
                            flank_match.cost,
                            bm.cost,
                            barcode.label.clone(),
                            bm.strand,
                        ));
                    }
                }
            }
        }
        results
        //collapse_overlapping_matches(&results, 0.5)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::search::interval::collapse_overlapping_matches;

    /// Helper function to check if a value is within wiggle room of expected
    fn assert_within_wiggle_room(actual: usize, expected: usize, wiggle_room: usize) {
        let min = expected.saturating_sub(wiggle_room);
        let max = expected + wiggle_room;
        assert!(
            actual >= min && actual <= max,
            "Expected {} to be within [{}, {}] of expected {}",
            actual,
            min,
            max,
            expected
        );
    }

    #[test]
    fn test_slice_masked_region() {
        let mut demuxer = Demuxer::new(0.5);
        let read = b"GGGGGAAATTTGGGCCCCCCCCCCCCCCCCCCCCCC".to_vec();
        //                         AAATTTGGG <- match (0 edits/cost)
        let flank = BarcodeGroup::new(
            &[b"AAATTTGGG", b"AAAXXXGGG"],
            &["s1", "s2"],
            BarcodeType::Fbar,
        );
        demuxer.add_query_group(flank);

        let matches = demuxer.demux(&read);
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].label, "s1");
        assert_eq!(matches[0].read_start_bar, 8);
        assert_eq!(matches[0].read_end_bar, 11); // Exclusive index
        assert_eq!(matches[0].bar_start, 0);
        assert_eq!(matches[0].bar_end, 3); // This is inclusive, maybe we should stay consistent exclusive?
        assert_eq!(matches[0].match_type, BarcodeType::Fbar);
        assert_eq!(matches[0].barcode_cost, 0);
        assert_eq!(matches[0].strand, Strand::Fwd);
    }

    #[test]
    fn test_demux_rc() {
        let mut demuxer = Demuxer::new(0.5);
        let read = b"GGGGGAAATTTGGGCCCCCCCCCCCCCCCCCCCCCC".to_vec();
        //                         AAATTTGGG <- match (0 edits/cost)
        let flank = BarcodeGroup::new(
            &[
                Iupac::reverse_complement("AAATTTGGG".as_bytes()).as_slice(),
                Iupac::reverse_complement("AAACCCGGG".as_bytes()).as_slice(),
            ],
            &["s1", "s2"],
            BarcodeType::Fbar,
        );
        demuxer.add_query_group(flank);

        let matches = demuxer.demux(&read);

        assert_eq!(matches.len(), 2);

        // Same as before, but now rc
        assert_eq!(matches[0].label, "s1");
        assert_eq!(matches[0].read_start_bar, 8);
        assert_eq!(matches[0].read_end_bar, 11); // Exclusive index
        assert_eq!(matches[0].bar_start, 0);
        assert_eq!(matches[0].bar_end, 3); // This is inclusive, maybe we should stay consistent exclusive?
        assert_eq!(matches[0].match_type, BarcodeType::Fbar);
        assert_eq!(matches[0].barcode_cost, 0);
        assert_eq!(matches[0].strand, Strand::Rc); // Same as above EXCEPT strand
    }

    #[test]
    fn test_just_overhang() {
        let mut demuxer = Demuxer::new(0.5);
        let read = b"TTTTTTGGGXXXXXXXXXXXXXXXXXXXXXXXXXXXx".to_vec();
        //                    |||||| -> 4 * 0.5 = 2 cost
        //             AAATTTTTTTTTTGGG <- match (0 edits/cost)
        //            TTTTTT
        //         TTTTTTTTTT
        let mut flank = BarcodeGroup::new(
            &[b"AAATTTTTTTTTTGGG", b"AAACCCCCCCCCCGGG"],
            &["s1", "s2"],
            BarcodeType::Fbar,
        );
        // Tune flank seq
        flank.k_cutoff = Some(4);
        // Tuen each of the barcodes within flank
        for barcode in flank.barcodes.iter_mut() {
            barcode.k_cutoff = Some(2);
        }

        demuxer.add_query_group(flank);
        let matches = demuxer.demux(&read);
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].label, "s1");

        // This can vary a little bit though
        assert_within_wiggle_room(matches[0].read_start_bar, 0, crate::WIGGLE_ROOM);
        assert_within_wiggle_room(matches[0].read_end_bar, 7, crate::WIGGLE_ROOM);

        assert_within_wiggle_room(matches[0].bar_start, 3, crate::WIGGLE_ROOM);
        assert_within_wiggle_room(matches[0].bar_end, 10, crate::WIGGLE_ROOM); // This is inclusive, maybe we should stay consistent exclusive?
        assert_eq!(matches[0].match_type, BarcodeType::Fbar);
        assert_eq!(matches[0].barcode_cost, 2);
        assert_eq!(matches[0].strand, Strand::Fwd); // Same as above EXCEPT strand
    }

    #[test]
    fn match_in_middle() {
        let mut demuxer = Demuxer::new(0.5);
        let read = b"CCCCCCCCCCCCCAAATTTTTTTTTTTGGGCCCCCCCCCCCC".to_vec();
        //                                 AAATTTTTTTTTTTGGG <- match (0 edits/cost)
        let flank = BarcodeGroup::new(
            &[b"AAATTTTTTTTTTTGGG", b"AAACCCCCCCCCCCGGG"],
            &["s1", "s2"],
            BarcodeType::Fbar,
        );
        demuxer.add_query_group(flank);

        let res = demuxer.demux(&read);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].label, "s1");
        assert_eq!(res[0].read_start_bar, 16);
        assert_eq!(res[0].read_end_bar, 27); // Exclusive index
        assert_eq!(res[0].bar_start, 0);
        assert_eq!(res[0].bar_end, 11); // This is inclusive, maybe we should stay consistent exclusive?
        assert_eq!(res[0].match_type, BarcodeType::Fbar);
        assert_eq!(res[0].barcode_cost, 0);
        assert_eq!(res[0].strand, Strand::Fwd);
    }

    #[test]
    pub fn search_real_read() {
        println!("Test real read");
        let read: Vec<u8> = b"TGTTATATTTCCCTGTACTTCGTTCCAGTTATTTTTATGCAAAAAACCGGTGTTTAACCACCACTGCCATGTATCAAAGTACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCAACAGGAAAACTATTTTCTGCAGG".to_vec();

        // FWD group
        let mut fwd_barcode_group = BarcodeGroup::new(
            &[b"TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCTTGGGTGTTTAACCACCACTGCCATGTATCAAAGTACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA", b"TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCTTGGGTGTTTAACCTTCGGATTCTATCGTGTTTCCCTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"],
            &["bar22_fwd", "bar4_fwd"],
            BarcodeType::Fbar,
        );
        fwd_barcode_group.tune_k(1000, 0.001, 0.5);
        for bar in fwd_barcode_group.barcodes.iter_mut() {
            bar.tune_k(1_000, 0.001, 0.5);
        }
        let mut demuxer = Demuxer::new(0.5);
        demuxer.add_query_group(fwd_barcode_group);
        println!("FWD demux");
        let fwd_matches = demuxer.demux(&read);
        println!("\n");

        // RC group
        let mut rc_barcode_group = BarcodeGroup::new(
            &[Iupac::reverse_complement("TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCTTGGGTGTTTAACCACCACTGCCATGTATCAAAGTACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA".as_bytes()).as_slice(), Iupac::reverse_complement("TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCTTGGGTGTTTAACCTTCGGATTCTATCGTGTTTCCCTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA".as_bytes()).as_slice()],
            &["bar22_rc", "bar4_rc"],
            BarcodeType::Fbar,
        );
        rc_barcode_group.tune_k(1000, 0.001, 0.5);
        for bar in rc_barcode_group.barcodes.iter_mut() {
            bar.tune_k(1_000, 0.001, 0.5);
        }
        let mut demuxer_rc = Demuxer::new(0.5);
        demuxer_rc.add_query_group(rc_barcode_group);
        println!("RC demux");
        let rc_matches = demuxer_rc.demux(&read);
        println!("\n");

        println!("FWD matches");
        for m in fwd_matches.iter() {
            println!("m: {:?}", m);
        }
        println!("RC matches");
        for m in rc_matches.iter() {
            println!("m: {:?}", m);
        }

        let collapsed_fwd_matches = collapse_overlapping_matches(&fwd_matches, 0.5);
        let collapsed_rc_matches = collapse_overlapping_matches(&rc_matches, 0.5);

        assert_eq!(collapsed_fwd_matches.len(), 1);
        assert_eq!(collapsed_rc_matches.len(), 1);
        let fwd_first = collapsed_fwd_matches[0].clone();
        let rc_first = collapsed_rc_matches[0].clone();

        assert_eq!(fwd_first.label, "bar22_fwd");
        assert_eq!(rc_first.label, "bar22_rc");
        // read starts
        assert_eq!(fwd_first.read_start_bar, 59);
        assert_eq!(rc_first.read_start_bar, 59);
        // Read ends
        assert_eq!(fwd_first.read_end_bar, 83);
        assert_eq!(rc_first.read_end_bar, 83);
        // Cost should be 0 for both
        assert_eq!(fwd_first.barcode_cost, 0);
        assert_eq!(rc_first.barcode_cost, 0);
    }
}
