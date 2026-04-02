use crate::annotate::barcodes::{BarcodeGroup, BarcodeType};
use crate::annotate::cigar_parse::*;
use crate::annotate::interval::collapse_overlapping_matches;
use crate::filter::pattern::Cut;
use cigar_lodhi_rs::*;
use pa_types::*;
use pa_types::{CigarOp, Cost};
use sassy::profiles::Iupac;
use sassy::{Match, Searcher, Strand};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

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
    // Buffers between demux calls
    results_buf: Vec<BarbellMatch>,
    best_per_pattern_buf: Vec<Option<Match>>,
    candidates_buf: Vec<(Match, usize)>, // (match, barcode index)
    scored_buf: Vec<(f64, f64, Match, usize)>, // (norm, raw, match, barcode index)
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
            results_buf: Vec::new(),
            best_per_pattern_buf: Vec::new(),
            candidates_buf: Vec::new(),
            scored_buf: Vec::new(),
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

    fn push_flank_only_result(
        results_buf: &mut Vec<BarbellMatch>,
        read_id: &str,
        read_len: usize,
        barcode_group: &BarcodeGroup,
        flank_match: &Match,
    ) {
        results_buf.push(BarbellMatch::new(
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
            read_len,
            read_id.to_string(),
            rel_dist_to_end(flank_match.text_start as isize, read_len),
            None,
        ));
    }

    fn collect_candidates_for_region(
        regular_searcher: &mut Searcher<Iupac>,
        best_per_pattern_buf: &mut Vec<Option<Match>>,
        candidates_buf: &mut Vec<(Match, usize)>,
        barcode_group: &BarcodeGroup,
        flank_match: &Match,
        barcode_region: &[u8],
        k_cutoff: usize,
    ) {
        let full_k_cutoff = barcode_group.barcodes[0].seq.len();

        // Keep only the best (lowest-cost) hit per encoded pattern index.
        best_per_pattern_buf.clear();
        best_per_pattern_buf.resize(barcode_group.barcodes.len(), None);

        for m in regular_searcher.search_encoded_patterns(
            barcode_group
                .encoded_barcodes
                .get_patterns(flank_match.strand),
            barcode_region,
            k_cutoff,
        ) {
            let idx: usize = m.pattern_idx;
            if idx >= best_per_pattern_buf.len() {
                continue;
            }

            let should_replace = match &best_per_pattern_buf[idx] {
                Some(best) => m.cost < best.cost,
                None => true,
            };
            if should_replace {
                best_per_pattern_buf[idx] = Some(m.clone());
            }
        }

        // Fallback: if only none or one pattern matched with strict cutoff, do a deeper pass.
        let matched_patterns = best_per_pattern_buf.iter().filter(|m| m.is_some()).count();

        if matched_patterns <= 1 && k_cutoff < full_k_cutoff {
            best_per_pattern_buf.fill(None);
            for m in regular_searcher.search_encoded_patterns(
                barcode_group
                    .encoded_barcodes
                    .get_patterns(flank_match.strand),
                barcode_region,
                full_k_cutoff,
            ) {
                let idx: usize = m.pattern_idx;
                if idx >= best_per_pattern_buf.len() {
                    continue;
                }

                let should_replace = match &best_per_pattern_buf[idx] {
                    Some(best) => m.cost < best.cost,
                    None => true,
                };
                if should_replace {
                    best_per_pattern_buf[idx] = Some(m.clone());
                }
            }
        }

        candidates_buf.clear();
        for (idx, maybe_match) in best_per_pattern_buf.iter_mut().enumerate() {
            if let Some(m) = maybe_match.take() {
                candidates_buf.push((m, idx));
            }
        }
    }

    fn score_and_push_result(
        lodhi: &mut Lodhi,
        scored_buf: &mut Vec<(f64, f64, Match, usize)>,
        candidates_buf: &mut Vec<(Match, usize)>,
        results_buf: &mut Vec<BarbellMatch>,
        read_id: &str,
        read_len: usize,
        barcode_group: &BarcodeGroup,
        flank_match: &Match,
        barcode_region_start: usize,
        perfect_bar_score: f64,
        min_score_frac: f64,
        min_score_diff_frac: f64,
    ) {
        if candidates_buf.is_empty() {
            Self::push_flank_only_result(
                results_buf,
                read_id,
                read_len,
                barcode_group,
                flank_match,
            );
            return;
        }

        // Score using Lodhi
        scored_buf.clear();
        for (m, barcode_idx) in candidates_buf.drain(..) {
            let s = lodhi.compute(&m.cigar);
            let s_norm = if perfect_bar_score > 0.0 {
                s / perfect_bar_score
            } else {
                0.0
            };
            scored_buf.push((s_norm, s, m, barcode_idx));
        }

        // sort by normalized score high to low
        scored_buf.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

        let (pad_start, _pad_end) = barcode_group.pad_region;
        let (bar_start, bar_end) = barcode_group.bar_region;
        let rel_bar_start = bar_start - pad_start;
        let rel_bar_end = bar_end - pad_start;

        let bar_read_region =
            map_pat_to_text_with_cost(&scored_buf[0].2, rel_bar_start as i32, rel_bar_end as i32);

        let ((bar_start, bar_end), (read_bar_start, read_bar_end), bar_cost) =
            bar_read_region.expect("No barcode match region found; unusual");

        // Apply fractional thresholds
        let top_norm = scored_buf[0].0;
        let mut is_valid_barcode_match = top_norm >= min_score_frac;
        if scored_buf.len() > 1 {
            is_valid_barcode_match =
                is_valid_barcode_match && (top_norm - scored_buf[1].0) >= min_score_diff_frac;
        }

        if is_valid_barcode_match {
            let top_barcode = &barcode_group.barcodes[scored_buf[0].3];
            results_buf.push(BarbellMatch::new(
                barcode_region_start + read_bar_start,
                barcode_region_start + read_bar_end,
                flank_match.text_start,
                flank_match.text_end,
                barcode_region_start + bar_start,
                barcode_region_start + bar_end,
                top_barcode.match_type.clone(),
                flank_match.cost,
                bar_cost as Cost,
                top_barcode.label.clone(),
                scored_buf[0].2.strand,
                read_len,
                read_id.to_string(),
                rel_dist_to_end(flank_match.text_start as isize, read_len),
                None,
            ));
        } else {
            Self::push_flank_only_result(
                results_buf,
                read_id,
                read_len,
                barcode_group,
                flank_match,
            );
        }
    }

    //fixme: would beneift from some more clean up
    /// Demultiplex read
    pub fn demux(&mut self, read_id: &str, read: &[u8]) -> Vec<BarbellMatch> {
        self.results_buf.clear();

        for (group_i, barcode_group) in self.queries.iter().enumerate() {
            let flank = &barcode_group.flank;
            let flank_k = barcode_group.k_cutoff.unwrap_or(0);
            let perfect_bar_score = self.perfect_scores[group_i];

            let flank_matches = self.overhang_searcher.search(flank, &read, flank_k);

            for flank_match in flank_matches.iter() {
                // Get barcode region
                let (mask_start, mask_end) = barcode_group.bar_region;

                // Extract read positions for barcode matching region
                let Some((barcode_region_start, barcode_region_end)) =
                    get_matching_region(flank_match, mask_start, mask_end)
                else {
                    continue; // no room for barcode?
                };

                // We have the barcode match region but to align against it we add some
                // padding on both sides of the barcode to anchor the alignment
                let barcode_region_start = barcode_region_start.saturating_sub(crate::PADDING);
                let barcode_region_end = (barcode_region_end + crate::PADDING).min(read.len());

                let barcode_region = &read[barcode_region_start..barcode_region_end];
                // Use v2 based searching now so we can search faster
                let barcode_len = barcode_group.barcodes[0].seq.len(); // includes padding like the encoded patterns
                // Lets for now take 0.4 * barcode_len this is already "random" treshold everything more is meaningless
                let k_cutoff = (barcode_len as f32 * 0.4) as usize;

                Self::collect_candidates_for_region(
                    &mut self.regular_searcher,
                    &mut self.best_per_pattern_buf,
                    &mut self.candidates_buf,
                    barcode_group,
                    flank_match,
                    barcode_region,
                    k_cutoff,
                );

                Self::score_and_push_result(
                    &mut self.lodhi,
                    &mut self.scored_buf,
                    &mut self.candidates_buf,
                    &mut self.results_buf,
                    read_id,
                    read.len(),
                    barcode_group,
                    flank_match,
                    barcode_region_start,
                    perfect_bar_score,
                    self.min_score_frac,
                    self.min_score_diff_frac,
                );
            }
        }

        collapse_overlapping_matches(&self.results_buf, 0.8)
    }
}
