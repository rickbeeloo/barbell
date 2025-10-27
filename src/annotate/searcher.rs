use crate::annotate::barcodes::{Barcode, BarcodeGroup, BarcodeType};
use crate::annotate::cigar_parse::*;
use crate::annotate::interval::collapse_overlapping_matches;
use crate::filter::pattern::Cut;
use cigar_lodhi_rs::*;
use pa_types::*;
use pa_types::{CigarOp, Cost};
use sassy::profiles::Iupac;
use sassy::{Match, Searcher, Strand};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
// remove later
use crate::annotate::edit_model::*;

pub struct Demuxer {
    queries: Vec<BarcodeGroup>,
    #[allow(dead_code)]
    verbose: bool,
    // Fractional thresholds 0..1
    min_score_frac: f64,
    min_score_diff_frac: f64,
    // Lodhi parameters used throughout
    lodhi: Lodhi,
    perfect_scores: Vec<f64>, // query idx > perfect score
    overhang_searcher: Searcher<Iupac>,
    regular_searcher: Searcher<Iupac>,
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
        _ => Err(serde::de::Error::custom(format!("Invalid strand: {s}"))),
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
                .map(|(cut, pos)| format!("{cut}:{pos}"))
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
                serde::de::Error::custom(format!("Invalid cut string: {cut_str}"))
            })?;
            let pos = pos_str
                .parse::<usize>()
                .map_err(|_| serde::de::Error::custom(format!("Invalid position: {pos_str}")))?;

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

/// Relative distance to read end
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
    pub fn new(alpha: f32, verbose: bool, min_score_frac: f64, min_score_diff_frac: f64) -> Self {
        Self {
            queries: Vec::new(),
            verbose,
            min_score_frac,
            min_score_diff_frac,
            perfect_scores: Vec::new(),
            lodhi: Lodhi::new(3, 0.5),
            overhang_searcher: Searcher::<Iupac>::new_rc_with_overhang(alpha),
            regular_searcher: Searcher::<Iupac>::new_rc(),
        }
    }

    /// Add query group to demuxer
    pub fn add_query_group(&mut self, query_group: BarcodeGroup) -> &mut Self {
        // Get perfect score for this group and store it
        let perfect_score = self.get_perfect_match_score(&query_group);
        self.queries.push(query_group);
        self.perfect_scores.push(perfect_score);
        self
    }

    /// Get perfect score for barcodegroup based on lodhi
    fn get_perfect_match_score(&mut self, barcode_group: &BarcodeGroup) -> f64 {
        // Todo: this does include padding, we might not always want that
        let l_bar = barcode_group.pad_region.1 - barcode_group.pad_region.0;
        let perfect = Cigar {
            ops: vec![CigarElem {
                op: CigarOp::Match,
                cnt: l_bar as i32,
            }],
        };
        self.lodhi.compute(&perfect)
    }

    //fixme: would beneift from some more clean up
    /// Demultiplex read
    pub fn demux(
        &mut self,
        read_id: &str,
        read: &[u8],
        search_lonely_bars: bool,
    ) -> Vec<BarbellMatch> {
        let mut results: Vec<BarbellMatch> = Vec::new();

        if search_lonely_bars {
            for (group_i, barcode_group) in self.queries.iter().enumerate() {
                for barcode_and_flank in barcode_group.barcodes.iter() {
                    let just_bar = &barcode_and_flank.seq
                        [crate::PADDING..barcode_and_flank.seq.len() - crate::PADDING];
                    let cut_off = get_edit_cut_off(just_bar.len());

                    // Looking for just barcode match
                    let just_bar_region_match =
                        self.regular_searcher.search(just_bar, &read, cut_off);

                    for flank_match in just_bar_region_match {
                        results.push(BarbellMatch::new(
                            flank_match.text_start,
                            flank_match.text_end,
                            flank_match.text_start,
                            flank_match.text_end,
                            flank_match.pattern_start,
                            flank_match.pattern_end,
                            barcode_group.barcode_type.as_bar(),
                            flank_match.cost,
                            flank_match.cost, // Keep same cost as barcode for now
                            barcode_and_flank.label.clone(),
                            flank_match.strand,
                            read.len(),
                            read_id.to_string(),
                            rel_dist_to_end(flank_match.text_start as isize, read.len()),
                            None,
                        ));
                    }
                }
            }
        }

        for (group_i, barcode_group) in self.queries.iter().enumerate() {
            let flank = &barcode_group.flank;
            let flank_k = barcode_group.k_cutoff.unwrap_or(0);
            let perfect_bar_score = self.perfect_scores[group_i];

            let flank_matches = self.overhang_searcher.search(flank, &read, flank_k);

            // If the search was a Fflank/Rflank we only have to search the flanks and not extract the barcode region
            if barcode_group.barcode_type == BarcodeType::Fadapter
                || barcode_group.barcode_type == BarcodeType::Radapter
            {
                for flank_match in flank_matches.iter() {
                    results.push(BarbellMatch::new(
                        flank_match.text_start,
                        flank_match.text_end,
                        flank_match.text_start,
                        flank_match.text_end,
                        0,
                        0,
                        barcode_group.barcode_type.clone(),
                        flank_match.cost,
                        flank_match.cost, // Keep same cost as barcode for now
                        "flank".to_string(),
                        flank_match.strand,
                        read.len(),
                        read_id.to_string(),
                        rel_dist_to_end(flank_match.text_start as isize, read.len()),
                        None,
                    ));
                }
            }

            for flank_match in flank_matches.iter() {
                // Get barcode region
                let (mut mask_start, mut mask_end) = barcode_group.bar_region;

                // If the match is in reverse complement we have to flip the locations based
                // on the length
                if flank_match.strand == Strand::Rc {
                    let flank_len = flank.len();
                    let (start, end) = barcode_group.bar_region;
                    mask_start = flank_len - end;
                    mask_end = flank_len - start;
                }

                // Extract read positions for barcode matching region
                let Some((barcode_region_start, barcode_region_end)) =
                    get_matching_region(flank_match, mask_start, mask_end)
                else {
                    continue; // no room for barcode?
                };

                // We have the barcode match region but to align against it we add some
                // padding on both sides of the barcode to anchor the alignment
                let barcode_region_start =
                    (flank_match.text_start + barcode_region_start).saturating_sub(crate::PADDING);
                let barcode_region_end =
                    (flank_match.text_start + barcode_region_end + crate::PADDING).min(read.len());

                let barcode_region = &read[barcode_region_start..barcode_region_end];
                let mut candidates: Vec<(Match, &Barcode)> = Vec::new();

                for barcode_and_flank in barcode_group.barcodes.iter() {
                    let fwd_best_hit = self
                        .regular_searcher
                        .search(
                            &barcode_and_flank.seq,
                            &barcode_region,
                            barcode_and_flank.seq.len(),
                        )
                        .into_iter()
                        .filter(|m| m.strand == flank_match.strand)
                        .min_by_key(|m| m.cost);

                    if let Some(fwd_hit) = fwd_best_hit {
                        candidates.push((fwd_hit, barcode_and_flank));
                    }
                }

                if candidates.is_empty() {
                    results.push(BarbellMatch::new(
                        flank_match.text_start,
                        flank_match.text_end,
                        flank_match.text_start,
                        flank_match.text_end,
                        0,
                        0,
                        barcode_group.barcodes[0].match_type.as_flank().clone(),
                        flank_match.cost,
                        barcode_group.barcodes[0].seq.len() as Cost,
                        "flank".to_string(),
                        flank_match.strand,
                        read.len(),
                        read_id.to_string(),
                        rel_dist_to_end(flank_match.text_start as isize, read.len()),
                        None,
                    ));

                    continue;
                }

                // Score using Lodhi
                let mut scored: Vec<(f64, f64, Match, &Barcode)> =
                    Vec::with_capacity(candidates.len());

                for (m, b) in candidates.into_iter() {
                    let s = self.lodhi.compute(&m.cigar);
                    let s_norm = if perfect_bar_score > 0.0 {
                        s / perfect_bar_score
                    } else {
                        0.0
                    };
                    scored.push((s_norm, s, m, b));
                }

                // sort by normalized score high to low
                scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

                /*
                    Important here is the structure:
                    ---------------------------------------- seq
                        |                |     pad_start, pad_end <- the query in the search
                         ===          ===      "padding"
                            |         |        bar_start, bar_end (here 'mask')


                So now for the query (including the padding) we want to extract the matching region
                exlcuding this padding to get the barcode match
                while we cannot ensure we have padding in our target, this is fixed in our barcode+flank query

                */

                let (pad_start, pad_end) = barcode_group.pad_region;
                let (bar_start, bar_end) = barcode_group.bar_region;
                let rel_bar_start = bar_start - pad_start;
                let rel_bar_end = pad_end - bar_end + (bar_end - bar_start);

                let bar_read_region = map_pat_to_text(
                    &scored[0].2,
                    rel_bar_start as i32,
                    rel_bar_end as i32,
                    scored[0].2.strand == Strand::Rc,
                );

                let ((bar_start, bar_end), (read_bar_start, read_bar_end)) =
                    bar_read_region.expect("No barcode match region found; unusual");

                // Apply fractional thresholds
                let top_norm = scored[0].0;
                let mut is_valid_barcode_match = top_norm >= self.min_score_frac;
                if scored.len() > 1 {
                    is_valid_barcode_match = is_valid_barcode_match
                        && (top_norm - scored[1].0) >= self.min_score_diff_frac;
                }

                if is_valid_barcode_match {
                    results.push(BarbellMatch::new(
                        barcode_region_start + read_bar_start,
                        barcode_region_start + read_bar_end,
                        flank_match.text_start,
                        flank_match.text_end,
                        barcode_region_start + bar_start,
                        barcode_region_start + bar_end,
                        scored[0].3.match_type.clone(),
                        flank_match.cost,
                        scored[0].2.cost as Cost,
                        scored[0].3.label.clone(),
                        scored[0].2.strand,
                        read.len(),
                        read_id.to_string(),
                        rel_dist_to_end(flank_match.text_start as isize, read.len()),
                        None,
                    ));
                } else {
                    results.push(BarbellMatch::new(
                        flank_match.text_start,
                        flank_match.text_end,
                        flank_match.text_start,
                        flank_match.text_end,
                        0,
                        0,
                        barcode_group.barcodes[0].match_type.as_flank().clone(),
                        flank_match.cost,
                        scored[0].2.cost as Cost,
                        "flank".to_string(),
                        flank_match.strand,
                        read.len(),
                        read_id.to_string(),
                        rel_dist_to_end(flank_match.text_start as isize, read.len()),
                        None,
                    ));
                }
            }
        }

        let results = collapse_overlapping_matches(&results, 0.99);
        results
    }
}

mod tests {
    use super::*;
    use crate::annotate::barcodes::*;
    use rand::prelude::*;
    use rand::seq::SliceRandom;
    use std::hint::black_box;

    fn edit(seq: &[u8], num_errors: usize) -> Vec<u8> {
        let mut rng = thread_rng();
        let mut edited = seq.to_vec();

        for _ in 0..num_errors {
            // Skip if empty
            if edited.is_empty() {
                break;
            }

            // Choose operation
            match rng.gen_range(0..3) {
                // Substitution
                0 => {
                    let bases = b"ATGC";
                    let pos = rng.gen_range(0..edited.len());
                    let new_base = *bases.choose(&mut rng).unwrap();
                    edited[pos] = new_base;
                }

                // Insertion
                1 => {
                    let bases = b"ATGC";
                    let new_base = *bases.choose(&mut rng).unwrap();
                    let pos = rng.gen_range(0..=edited.len()); // insertion *can* go at the end
                    edited.insert(pos, new_base);
                }

                // Deletion
                2 => {
                    if !edited.is_empty() {
                        let pos = rng.gen_range(0..edited.len());
                        edited.remove(pos);
                    }
                }

                _ => unreachable!(),
            }
        }

        edited
    }

    fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        let mut rng = thread_rng();
        let mut seq = seq.to_vec();
        seq.reverse();
        seq.iter_mut().for_each(|base| {
            *base = match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => unreachable!(),
            }
        });
        seq
    }

    #[test]
    #[ignore]
    fn fuzz_demuxer() {
        let mut rng = thread_rng();

        let mut demuxer = Demuxer::new(0.8, false, 0.2, 0.1);
        let query_groups: Vec<BarcodeGroup> = BarcodeGroup::new_from_kit("SQK-RBK114-96", false);
        for group in query_groups {
            demuxer.add_query_group(group);
        }

        let left_flank: &[u8] = b"GCTTGGGTGTTTAACC";
        let barcode: &[u8] = b"AAGAAAGTTGTCGGTGTCTTTGTG";
        let right_flank: &[u8] = b"GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";

        let parts = vec![left_flank, barcode, right_flank];
        let attempts = 1_000_000;

        for _ in 0..attempts {
            // Choose a random subset of parts (size 1–3)
            let num_parts = rng.gen_range(1..=3);
            let mut chosen_parts = parts
                .choose_multiple(&mut rng, num_parts)
                .cloned()
                .collect::<Vec<_>>();

            // Random order
            chosen_parts.shuffle(&mut rng);

            // Concatenate
            let combined: Vec<u8> = chosen_parts
                .iter()
                .flat_map(|p| p.iter().copied())
                .collect();

            // Apply random number of edits (0–80)
            let num_edits = rng.gen_range(0..=80);
            let mut mutated = edit(&combined, num_edits);

            // Lets see if we rc it
            let rc_bool = rng.gen_bool(0.5);
            if rc_bool {
                mutated = reverse_complement(&mutated);
            }

            // Slice only first 35
            let slice_end = 35.min(mutated.len());
            mutated = mutated[..slice_end].to_vec();

            // Run the demuxer and consume result
            let matches = demuxer.demux("fuzzed", &mutated, false);
            black_box(matches);
        }
    }
}
