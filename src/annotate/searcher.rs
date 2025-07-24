use crate::annotate::barcodes::{BarcodeGroup, BarcodeType};
use crate::annotate::interval::collapse_overlapping_matches;
use crate::filter::pattern::Cut;
use colored::*;
use needletail::{FastxReader, Sequence, parse_fastx_file};
use pa_types::Cost;
use sassy::profiles::Iupac;
use sassy::profiles::Profile;
use sassy::{Match, Searcher, Strand};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::thread_local;

thread_local! {
    static OVERHANG_SEARCHER: std::cell::RefCell<Option<Searcher<Iupac>>> = std::cell::RefCell::new(None);
}

#[derive(Clone)]
pub struct Demuxer {
    alpha: f32,
    queries: Vec<BarcodeGroup>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BarbellMatch {
    pub read_id: String,
    pub read_len: usize,
    pub rel_dist_to_end: isize,

    // Where barcode starts in read
    pub read_start_bar: usize,
    pub read_end_bar: usize,

    // Where flank starts in read
    pub read_start_flank: usize,
    pub read_end_flank: usize,

    // Where in the barcode the match starts
    pub bar_start: usize,
    pub bar_end: usize,

    pub match_type: BarcodeType,
    pub flank_cost: Cost,
    pub barcode_cost: Cost,
    pub label: String,
    #[serde(
        serialize_with = "serialize_strand",
        deserialize_with = "deserialize_strand"
    )]
    pub strand: Strand,

    #[serde(
        serialize_with = "serialize_cuts",
        deserialize_with = "deserialize_cuts"
    )]
    pub cuts: Option<Vec<(Cut, usize)>>,
}

// Custom serialization for strand field
fn serialize_strand<S>(strand: &Strand, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    match strand {
        Strand::Fwd => serializer.serialize_str("Fwd"),
        Strand::Rc => serializer.serialize_str("Rc"),
    }
}

// Custom deserialization for strand field
fn deserialize_strand<'de, D>(deserializer: D) -> Result<Strand, D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = String::deserialize(deserializer)?;
    match s.as_str() {
        "Fwd" => Ok(Strand::Fwd),
        "Rc" => Ok(Strand::Rc),
        _ => Err(serde::de::Error::custom(format!("Invalid strand: {}", s))),
    }
}

// Custom serialization for cuts field
fn serialize_cuts<S>(cuts: &Option<Vec<(Cut, usize)>>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    match cuts {
        None => serializer.serialize_str(""),
        Some(cuts_vec) => {
            let cuts_str = cuts_vec
                .iter()
                .map(|(cut, pos)| format!("{}:{}", cut, pos))
                .collect::<Vec<_>>()
                .join(",");
            serializer.serialize_str(&cuts_str)
        }
    }
}

// Custom deserialization for cuts field
fn deserialize_cuts<'de, D>(deserializer: D) -> Result<Option<Vec<(Cut, usize)>>, D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = String::deserialize(deserializer)?;

    if s.is_empty() {
        return Ok(None);
    }

    let cuts_vec: Result<Vec<_>, _> = s
        .split(',')
        .map(|part| {
            let mut split = part.split(':');
            let cut_str = split
                .next()
                .ok_or_else(|| serde::de::Error::custom("Invalid cut format: missing cut part"))?;
            let pos_str = split.next().ok_or_else(|| {
                serde::de::Error::custom("Invalid cut format: missing position part")
            })?;

            let cut = Cut::from_string(cut_str).ok_or_else(|| {
                serde::de::Error::custom(format!("Invalid cut string: {}", cut_str))
            })?;
            let pos = pos_str
                .parse::<usize>()
                .map_err(|_| serde::de::Error::custom(format!("Invalid position: {}", pos_str)))?;

            Ok((cut, pos))
        })
        .collect();

    cuts_vec.map(Some)
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
        read_len: usize,
        read_id: String,
        rel_dist_to_end: isize,
        cuts: Option<Vec<(Cut, usize)>>,
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
            read_len,
            read_id,
            rel_dist_to_end,
            cuts,
        }
    }
}

fn rel_dist_to_end(pos: isize, read_len: usize) -> isize {
    // If already negative, for a match starting before the read we can return 1
    if pos < 0 {
        return 1;
    }

    if pos <= (read_len / 2) as isize {
        if pos == 0 {
            1 // Left end
        } else {
            pos // Distance from left end (already isize)
        }
    } else if pos == read_len as isize {
        -1 // Right end
    } else {
        -(read_len as isize - pos) // Distance from right end
    }
}

impl Demuxer {
    pub fn new(alpha: f32) -> Self {
        Self {
            alpha,
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
        let flank_len = flank.flank.len();
        let (mask_start, mask_end) = flank.bar_region;
        // println!("Matching orientation: {:?}", flank_match.strand);
        // println!(
        //     "Flank region: {}",
        //     String::from_utf8_lossy(&read[flank_match.text_start..flank_match.text_end])
        // );

        // println!("Mask start - mask end: {} - {}", mask_start, mask_end);
        // let (mask_start, mask_end) = if flank_match.strand == Strand::Rc {
        //     println!("Flank len: {}", flank_len);
        //     println!(
        //         "Adjusting mask start - mask end: {} - {}",
        //         flank_len - mask_start,
        //         flank_len - mask_end
        //     );
        //     (flank_len - mask_end, flank_len - mask_start)
        // } else {
        //     (mask_start, mask_end)
        // };

        //  (&read[mask_start..=mask_end], (mask_start, mask_end))

        // Reverse cigar
        let mut flank_match = flank_match.clone();
        flank_match.cigar.reverse();

        let path = flank_match.to_path();
        // for pos in path.iter() {
        //     println!("Q - R: {} - {}", pos.0, pos.1);
        // }

        //  println!("Mask start - mask end: {} - {}", mask_start, mask_end);

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
        let end = (max_r_pos + crate::WIGGLE_ROOM).min(read.len() - 1);

        // println!(
        //     "Mask region: [{}-{}] c: {}{:?}",
        //     start,
        //     end,
        //     flank_match.cost,
        //     String::from_utf8_lossy(&read[start..=end])
        // );
        (&read[start..=end], (start, end))
    }

    pub fn demux(&mut self, read_id: &str, read: &[u8]) -> Vec<BarbellMatch> {
        let start_time = std::time::Instant::now();
        let mut results: Vec<BarbellMatch> = Vec::new();

        // Initialize thread-local searcher if not already done
        OVERHANG_SEARCHER.with(|cell| {
            if cell.borrow().is_none() {
                *cell.borrow_mut() = Some(Searcher::<Iupac>::new_rc_with_overhang(self.alpha));
            }
        });

        // We first search for the top level of the barcodeGroups
        for (i, barcode_group) in self.queries.iter().enumerate() {
            //println!("Barcode group: {}", i);

            let flank = &barcode_group.flank;
            //  println!("Flank: {}", String::from_utf8_lossy(flank));
            let k = barcode_group.k_cutoff.unwrap_or(0);
            // println!("K: {}", k);
            // println!("Q: {}", String::from_utf8_lossy(flank));
            // println!("Read: {}", String::from_utf8_lossy(read));

            //  println!("Looking for flank: {}", String::from_utf8_lossy(flank));
            let flank_matches = OVERHANG_SEARCHER.with(|cell| {
                if let Some(ref mut searcher) = *cell.borrow_mut() {
                    searcher.search(flank, &read, k)
                } else {
                    unreachable!();
                    //vec![]
                }
            });
            // println!("\tFlank matches: {}", flank_matches.len());

            // If we found a flank, we slice out the masked region and search for the barcodes in the
            // masked region, we should be AWARE OF THE STRAND as sassy now returns both directions (Fwd, Rc)
            for flank_match in flank_matches {
                // println!(
                //     "[{}] Flank match: {}",
                //     read_id,
                //     String::from_utf8_lossy(
                //         &read[flank_match.start.1 as usize..flank_match.end.1 as usize]
                //     )
                // );
                // println!(
                //     "[{}] Flank match RC: {}",
                //     read_id,
                //     String::from_utf8_lossy(&Iupac::reverse_complement(
                //         &read[flank_match.start.1 as usize..flank_match.end.1 as usize]
                //     ))
                // );
                let (mask_region, (mask_start, _mask_end)) =
                    self.slice_masked_region(read, barcode_group, &flank_match);

                // If no mask match found we can just skip to the next one
                if mask_region.is_empty() {
                    continue;
                }

                // println!(
                //     "Mask region FWD: {:?}",
                //     String::from_utf8_lossy(mask_region)
                // );
                // println!(
                //     "Mask region RC: {:?}",
                //     String::from_utf8_lossy(&Iupac::reverse_complement(mask_region))
                // );
                // Then we compare each query against the masked region
                let mut barcode_found = false;
                // println!("Number of barcodes: {}", barcode_group.barcodes.len());
                for barcode in barcode_group.barcodes.iter() {
                    // println!(
                    //     "\tLooking for barcode: {} with seq: {}",
                    //     barcode.label,
                    //     String::from_utf8_lossy(&barcode.seq)
                    // );
                    // Look for barcode mathces in mask slice
                    let bms = OVERHANG_SEARCHER.with(|cell| {
                        let k = barcode.k_cutoff.unwrap_or(0);
                        //  println!("Using k: {}", k);
                        if let Some(ref mut searcher) = *cell.borrow_mut() {
                            // println!("\tLooking for barcode: {}", barcode.label);

                            searcher.search(&barcode.seq, &mask_region, k)
                        } else {
                            vec![]
                        }
                    });

                    // We make sure they flank and barcode match on the same strand
                    for bm in bms.iter().filter(|bm| bm.strand == flank_match.strand) {
                        // println!("\tBm cost: {}", bm.cost);

                        barcode_found = true;
                        let rel_dist = rel_dist_to_end(flank_match.text_start as isize, read.len());

                        results.push(BarbellMatch::new(
                            //t - storing flank matches to more accurately filter out overlaps later on
                            bm.text_start as usize + mask_start,
                            bm.text_end as usize + mask_start,
                            flank_match.text_start as usize,
                            flank_match.text_end as usize,
                            //q
                            bm.pattern_start as usize,
                            bm.pattern_end as usize,
                            barcode.match_type.clone(),
                            flank_match.cost,
                            bm.cost,
                            barcode.label.clone(),
                            bm.strand,
                            read.len(),
                            read_id.to_string(),
                            rel_dist,
                            None, // Only added after filtering is used
                        ));
                    }
                }
                if !barcode_found {
                    // println!("Came here!");
                    // We just add the flank as a match
                    results.push(BarbellMatch::new(
                        //t - storing flank matches to more accurately filter out overlaps later on
                        flank_match.text_start as usize + mask_start,
                        flank_match.text_end as usize + mask_start,
                        flank_match.text_start as usize,
                        flank_match.text_end as usize,
                        //q
                        0,
                        0,
                        barcode_group.barcodes[0].match_type.as_flank().clone(),
                        flank_match.cost,
                        barcode_group.barcodes[0].seq.len() as Cost,
                        "flank".to_string(),
                        flank_match.strand,
                        read.len(),
                        read_id.to_string(),
                        rel_dist_to_end(
                            (flank_match.text_start as usize + mask_start) as isize,
                            read.len(),
                        ),
                        None, // Only added after filtering is used
                    ));
                }
            }
        }
        let end_time = std::time::Instant::now();
        //  println!("Time taken: {:?}", end_time.duration_since(start_time));
        // results
        collapse_overlapping_matches(&results, 0.5)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotate::interval::collapse_overlapping_matches;

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
            vec![b"AAATTTGGG".to_vec(), b"AAAXXXGGG".to_vec()],
            vec!["s1".to_string(), "s2".to_string()],
            BarcodeType::Fbar,
        );
        demuxer.add_query_group(flank);

        let matches = demuxer.demux("test", &read);
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
        let read = b"GGGGGAAATTTTTTTTTTTGGGCCCCCCCCCCCCCCCCCCCCCC".to_vec();
        //                         AAATTTGGG <- match (0 edits/cost)
        let mut bar_group = BarcodeGroup::new(
            vec![
                Iupac::reverse_complement("AAATTTTTTTTTTTGGG".as_bytes()).to_vec(),
                Iupac::reverse_complement("AAACCCCCCCCCCCGGG".as_bytes()).to_vec(),
            ],
            vec!["s1".to_string(), "s2".to_string()],
            BarcodeType::Fbar,
        );
        //bar_group.tune_group_random_sequences(1000, 0.001, 0.5, 1);
        bar_group.tune_set_perc_threshold(0.1);

        demuxer.add_query_group(bar_group);

        let matches = demuxer.demux("test", &read);

        assert!(matches.len() < 3 && matches.len() > 0);

        // Same as before, but now rc
        assert_eq!(matches[0].label, "s1");
        assert_eq!(matches[0].read_start_bar, 8);
        assert_eq!(matches[0].read_end_bar, 19); // Exclusive index
        assert_eq!(matches[0].bar_start, 0);
        assert_eq!(matches[0].bar_end, 11); // This is inclusive, maybe we should stay consistent exclusive?
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
            vec![b"AAATTTTTTTTTTGGG".to_vec(), b"AAACCCCCCCCCCGGG".to_vec()],
            vec!["s1".to_string(), "s2".to_string()],
            BarcodeType::Fbar,
        );
        // Tune flank seq
        flank.k_cutoff = Some(4);
        // Tuen each of the barcodes within flank
        for barcode in flank.barcodes.iter_mut() {
            barcode.k_cutoff = Some(2);
        }

        demuxer.add_query_group(flank);
        let matches = demuxer.demux("test", &read);
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
            vec![b"AAATTTTTTTTTTTGGG".to_vec(), b"AAACCCCCCCCCCCGGG".to_vec()],
            vec!["s1".to_string(), "s2".to_string()],
            BarcodeType::Fbar,
        );
        demuxer.add_query_group(flank);

        let res = demuxer.demux("test", &read);
        assert_eq!(res.len(), 2);
        assert_eq!(res[1].label, "s1");
        assert_eq!(res[1].read_start_bar, 16);
        assert_eq!(res[1].read_end_bar, 27); // Exclusive index
        assert_eq!(res[1].bar_start, 0);
        assert_eq!(res[1].bar_end, 11); // This is inclusive, maybe we should stay consistent exclusive?
        assert_eq!(res[1].match_type, BarcodeType::Fbar);
        assert_eq!(res[1].barcode_cost, 0);
        assert_eq!(res[1].strand, Strand::Fwd);
    }

    #[test]
    fn test_multi_match() {
        let mut demuxer = Demuxer::new(0.5);
        let read = b"CCCCCCCCCCCCCAAATTTTTTTTTTTGGGCCCCCCCCCCCC".to_vec();
        //                                 AAATTTTTTTTTTTGGG <- match (0 edits/cost)
        let flank = BarcodeGroup::new(
            vec![b"AAATTTTTTTTTTTGGG".to_vec(), b"AAACCCCCCCCCCCGGG".to_vec()],
            vec!["s1".to_string(), "s2".to_string()],
            BarcodeType::Fbar,
        );
        demuxer.add_query_group(flank);
    }

    #[test]
    pub fn search_real_read() {
        println!("Test real read");
        let read: Vec<u8> = b"TGTTATATTTCCCTGTACTTCGTTCCAGTTATTTTTATGCAAAAAACCGGTGTTTAACCACCACTGCCATGTATCAAAGTACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCAACAGGAAAACTATTTTCTGCAGG".to_vec();

        // FWD group
        let mut fwd_barcode_group = BarcodeGroup::new(
            vec![
                b"TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCTTGGGTGTTTAACCACCACTGCCATGTATCAAAGTACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA".to_vec(),
                b"TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCTTGGGTGTTTAACCTTCGGATTCTATCGTGTTTCCCTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA".to_vec()
            ],
            vec!["bar22_fwd".to_string(), "bar4_fwd".to_string()],
            BarcodeType::Fbar,
        );
        fwd_barcode_group.tune_group_random_sequences(1000, 0.001, 0.5, 4);
        let mut demuxer = Demuxer::new(0.5);
        demuxer.add_query_group(fwd_barcode_group);
        println!("FWD demux");
        let fwd_matches = demuxer.demux("test", &read);
        println!("\n");

        // RC group
        let mut rc_barcode_group = BarcodeGroup::new(
            vec![
                Iupac::reverse_complement("TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCTTGGGTGTTTAACCACCACTGCCATGTATCAAAGTACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA".as_bytes()).to_vec(),
                Iupac::reverse_complement("TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCTTGGGTGTTTAACCTTCGGATTCTATCGTGTTTCCCTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA".as_bytes()).to_vec()
            ],
            vec!["bar22_rc".to_string(), "bar4_rc".to_string()],
            BarcodeType::Fbar,
        );
        rc_barcode_group.tune_group_random_sequences(1000, 0.001, 0.5, 4);
        let mut demuxer_rc = Demuxer::new(0.5);
        demuxer_rc.add_query_group(rc_barcode_group);
        println!("RC demux");
        let rc_matches = demuxer_rc.demux("test", &read);
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

    #[test]
    fn test_overhanging_real_read() {
        let read = b"ATGTTTTTTTTTTTTGCCGATATAACCGTTTCATATCGGAGGGAATGGAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCATCAGCTTCAACAGCATCAGCACAGTTTTGAACCTGATATGGATGAGTTTCGTGAACTGCGTAATTTTAACTTGGGTGATAGCTTACACGCCGTACGTGGAAGCAGGCCGCGCGTGGGCAAGGTTTATATATTAAAGTCTTTGAGCAGCATAATGATCAACCGAGTATGGAAATCCATTATGCAAACATGCAAGTACCGAGCATGAAGAGAAATTGAGCTTAATGATGGGGCTGATTGAGCAATGTGAGCAACTGCAATGTAGCTATGCGGTATTTTTGCCTCAAGCTCGATTAGCTGTAGGCACAGGTGCAAACCAGTTGCTTCAGGCTAAAAAACTTCTGGCGCAGGCTTAATTATGATGCATGCTGATGTAGAACCTCGATTTAATGACTTTGCTGTTATTACTGCTTGCACAAGTCTTGTTTAATCCCGTTCTGTTAACTGCAATTTTTATTCTCCTCTGTTTATTTGTCAGCTTTAAAGATGAAACAAAAACAGTATCCAA";
        let bar10 = b"TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCTTGGGTGTTTAACCACCACTGCCATGTATCAAAGTACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";
        let bar11 = b"TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCTTGGGTGTTTAACCCGTCAACTGACAGTGGTTCGTACTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";

        let mut demuxer = Demuxer::new(0.5);
        let mut fwd_barcode_group = BarcodeGroup::new(
            vec![bar10.to_vec(), bar11.to_vec()],
            vec!["bar22".to_string(), "bar16".to_string()],
            BarcodeType::Fbar,
        );

        let tune_start_time = std::time::Instant::now();
        fwd_barcode_group.tune_group_random_sequences(10_000, 0.001, 0.5, 4);
        let tune_end_time = std::time::Instant::now();
        println!(
            "Tune time: {:?}",
            tune_end_time.duration_since(tune_start_time)
        );

        demuxer.add_query_group(fwd_barcode_group);

        let matches = demuxer.demux("test", read.as_slice());
        for m in matches.iter() {
            println!("m: {:?}", m);
        }
    }

    #[test]
    #[ignore = "Full test only for debugging"]
    fn test_full_rapid_flow() {
        let example_file = "examples/rapid_bars.fasta";
        let mut barcode_group = BarcodeGroup::new_from_fasta(example_file, BarcodeType::Fbar);
        barcode_group.tune_group_random_sequences(100_000, 0.0001, 0.5, 16);
        let mut demuxer = Demuxer::new(0.5);
        demuxer.add_query_group(barcode_group);

        // Read file test
        let mut reader =
            parse_fastx_file("/home/solprof/Downloads/sub.fastq").expect("valid path/file");
        while let Some(record) = reader.next() {
            let seqrec = record.expect("invalid record");
            let norm_seq = seqrec.normalize(false);
            println!("Read id: {}", String::from_utf8_lossy(seqrec.id()));
            let result = demuxer.demux("test", &norm_seq.into_owned());
            for m in result.iter() {
                println!("\tm: {:?}", m);
            }
            println!("\n");
        }
    }

    #[test]
    fn test_serialization_deserialization() {
        use crate::filter::pattern::{Cut, CutDirection};
        use serde_json;

        // Test with no cuts
        let match1 = BarbellMatch::new(
            0,
            10,
            0,
            10,
            0,
            10,
            BarcodeType::Fbar,
            0,
            0,
            "test".to_string(),
            Strand::Fwd,
            100,
            "read1".to_string(),
            0,
            None,
        );

        let json1 = serde_json::to_string(&match1).unwrap();
        let deserialized1: BarbellMatch = serde_json::from_str(&json1).unwrap();
        assert_eq!(match1.cuts, deserialized1.cuts);

        // Test with cuts
        let cuts = Some(vec![
            (Cut::new(1, CutDirection::Before), 50),
            (Cut::new(2, CutDirection::After), 75),
        ]);
        let match2 = BarbellMatch::new(
            0,
            10,
            0,
            10,
            0,
            10,
            BarcodeType::Fbar,
            0,
            0,
            "test".to_string(),
            Strand::Fwd,
            100,
            "read2".to_string(),
            0,
            cuts.clone(),
        );

        let json2 = serde_json::to_string(&match2).unwrap();
        let deserialized2: BarbellMatch = serde_json::from_str(&json2).unwrap();
        assert_eq!(match2.cuts, deserialized2.cuts);

        // Verify the JSON format
        assert!(json1.contains("\"cuts\":\"\""));
        assert!(json2.contains("\"cuts\":\"Before(1):50,After(2):75\""));
    }

    #[test]
    fn overflow_rel_dist_bug() {
        let pos = 167_isize;
        let read_len = 173_usize;
        let rel_pos = rel_dist_to_end(pos, read_len);
        println!("rel_pos: {}", rel_pos);
    }

    #[test]
    fn example_rc_miss() {
        let read = b"ACCAGATTGCTGGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCTTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAGATTAAGCCATGCATGTCTAAGTATAAACAAATTCATACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATAGTTTATTTGATGGTTTCTTGCTACATGGATAACTGTGGTAATTCTAGAGCTAATACATGCTGAAAAGCCCCGACTTCTGGAAGGGGTGTATTTATTAGATAAAAAACCAATGACTTCGGTCTTCTTGGTGATTCATAATAACTTCTCGAATCGCATGGCCTCGCGCCGGCGATGCTTCATTCAAATATCTGCCCTATCAACTTTCGATGGTAGGATAGAGGCCTACCATGGTATCAACGGGTAACGGGAATTAGGGTTCGATTCCGGAGAGGGAGCCTAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCCGACACGGGGA";
        let read_rc = Iupac::reverse_complement(read);

        let mut demuxer = Demuxer::new(0.5);

        /*

                >5R
        AATGTACTTCGTTCAGTTACGTATTGCTGGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCTTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTACCTGGTTGATYCTG
        >6R
        AATGTACTTCGTTCAGTTACGTATTGCTGGTGCTGTTCTCGCAAAGGCAGAAAGTAGTCTTAACCTTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTACCTGGTTGATYCTG
         */

        let mut fwd_barcode_group = BarcodeGroup::new(
            vec![
                b"AATGTACTTCGTTCAGTTACGTATTGCTGGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCTTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTACCTGGTTGATYCTG".to_vec(),
                b"AATGTACTTCGTTCAGTTACGTATTGCTGGTGCTGTTCTCGCAAAGGCAGAAAGTAGTCTTAACCTTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTACCTGGTTGATYCTG".to_vec()],
            vec!["5R".to_string(), "6R".to_string()],
            BarcodeType::Rbar,
        );
        fwd_barcode_group.tune_group_random_sequences(100, 0.001, 0.5, 16);
        demuxer.add_query_group(fwd_barcode_group);
        println!("FWD demux");
        let fwd_matches = demuxer.demux("test", read.as_slice());
        for fw_m in fwd_matches.iter() {
            println!("fw_m: {:?}", fw_m);
        }
        println!("\n");
        println!("RC demux");
        let rc_matches = demuxer.demux("test", read_rc.as_slice());
        for rc_m in rc_matches.iter() {
            println!("rc_m: {:?}", rc_m);
        }
        println!("\n");
    }

    #[test]
    fn test_rel_end() {
        let read_len = 3537;
        let pos = 3537;
        let rel_dist = rel_dist_to_end(pos, read_len);
        println!("rel dist: {}", rel_dist);
    }
}
