use crate::annotate::barcodes::{Barcode, BarcodeGroup, BarcodeType};
use crate::annotate::cigar_parse::*;
use crate::annotate::interval::collapse_overlapping_matches;
use crate::annotate::model::Model;
use crate::filter::pattern::Cut;
use colored::*;
use itertools::Itertools;
use needletail::{FastxReader, Sequence, parse_fastx_file};
use pa_types::{CigarOp, Cost, CostModel, Pos};
use sassy::profiles::Iupac;
use sassy::profiles::Profile;
use sassy::{Match, Searcher, Strand};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::thread_local;
// (CigarElem, CigarOp) were only required by the previous, now removed implementation
// of `extract_bar_cost`.  They are no longer used, so the import has been removed.
use crate::WIGGLE_ROOM;
use cigar_lodhi_rs::*;
use pa_types::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{Arc, Mutex};

thread_local! {
    static OVERHANG_SEARCHER: std::cell::RefCell<Option<Searcher<Iupac>>> = std::cell::RefCell::new(None);
    static REGULAR_SEARCHER: std::cell::RefCell<Option<Searcher<Iupac>>> = std::cell::RefCell::new(None);

}

// tod what if barcodes get the same lowest edits?

#[derive(Clone)]
pub struct Demuxer {
    alpha: f32,
    queries: Vec<BarcodeGroup>,
    verbose: bool,
    top2_writer: Option<Arc<Mutex<BufWriter<File>>>>,
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
            verbose: false,
            top2_writer: None,
        }
    }

    pub fn new_with_verbose(alpha: f32, verbose: bool) -> Self {
        Self {
            alpha,
            queries: Vec::new(),
            verbose,
            top2_writer: None,
        }
    }

    pub fn new_with_verbose_and_top2(
        alpha: f32,
        verbose: bool,
        top2_writer: Option<Arc<Mutex<BufWriter<File>>>>,
    ) -> Self {
        Self {
            alpha,
            queries: Vec::new(),
            verbose,
            top2_writer,
        }
    }

    pub fn add_query_group(&mut self, query_group: BarcodeGroup) -> &mut Self {
        self.queries.push(query_group);
        self
    }

    pub fn demux(&mut self, read_id: &str, read: &[u8]) -> Vec<BarbellMatch> {
        let mut results: Vec<BarbellMatch> = Vec::new();

        // Initialize thread-local searcher if not already done
        OVERHANG_SEARCHER.with(|cell| {
            if cell.borrow().is_none() {
                *cell.borrow_mut() = Some(Searcher::<Iupac>::new_rc_with_overhang(self.alpha));
            }
        });

        REGULAR_SEARCHER.with(|cell| {
            if cell.borrow().is_none() {
                *cell.borrow_mut() = Some(Searcher::<Iupac>::new_rc());
            }
        });

        for (i, barcode_group) in self.queries.iter().enumerate() {
            let flank = &barcode_group.flank;
            let flank_k = barcode_group.k_cutoff.unwrap_or(0);
            let flank_matches = OVERHANG_SEARCHER.with(|cell| {
                if let Some(ref mut searcher) = *cell.borrow_mut() {
                    searcher.search(flank, &read, flank_k)
                } else {
                    unreachable!();
                }
            });

            let verbose = self.verbose;

            for (i, flank_match) in flank_matches.iter().enumerate() {
                if verbose {
                    println!(
                        "Flank {i}: {}-{}",
                        flank_match.text_start, flank_match.text_end
                    );
                }

                // Get barcode region
                let (mask_start, mask_end) = barcode_group.bar_region;

                // Extract read positions for barcode matchign region
                let Some((barcode_region_start, barcode_region_end)) =
                    get_matching_region(&flank_match, mask_start, mask_end)
                else {
                    continue; // no room for barcode?
                };
                // let barcode_region_start =
                //     (flank_match.text_start + barcode_region_start).saturating_sub(10);
                // let barcode_region_end = (flank_match.text_start + barcode_region_end)
                //     .saturating_add(10 + 6)
                //     .min(read.len());

                let barcode_region_start =
                    (flank_match.text_start + barcode_region_start).saturating_sub(10 + 5);
                let barcode_region_end =
                    (flank_match.text_start + barcode_region_end + 10 + 5).min(read.len());

                let barcode_region = &read[barcode_region_start..barcode_region_end];

                if verbose {
                    for _ in 0..10 {
                        println!(
                            "Barcode region: {}",
                            String::from_utf8_lossy(barcode_region)
                        );
                    }
                }

                let mut candidates: Vec<(Match, &Barcode)> = Vec::new();

                for barcode_and_flank in barcode_group.barcodes.iter() {
                    // println!(
                    //     "Looking for: {}",
                    //     String::from_utf8_lossy(&barcode_and_flank.seq)
                    // );
                    // Search the forward version
                    let fwd_best_hit = REGULAR_SEARCHER.with(|cell| {
                        if let Some(ref mut searcher) = *cell.borrow_mut() {
                            searcher
                                .search(&barcode_and_flank.seq, &barcode_region, 20)
                                .into_iter()
                                .filter(|m| m.strand == flank_match.strand)
                                // iter mutate and we do region_len - (m.text_start + m.text_start)
                                .min_by_key(|m| m.cost)
                        } else {
                            None
                        }
                    });

                    // if a fwd hit add it
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

                // Look at top candidates
                let top_candidates: Vec<(Match, &Barcode)> = candidates.into_iter().collect();

                let mut lodhi = Lodhi::new(3, 0.5);
                // Now only calculate likelihood scores for the top candidates
                let mut scored: Vec<(f64, Match, &Barcode)> = top_candidates
                    .into_iter()
                    .map(|(m, b)| {
                        // // Convert to per base with a convex reward for contiguous match runs
                        // // Tuneable: prefer longer runs without overwhelming the base model
                        // let coeff = 0.05_f64; // small coefficient
                        // let exponent = 1.6_f64; // convex (>1)
                        // let avg_cost_per_pattern_base =
                        //     simple_model.avg_cost_per_pattern_base(&m.cigar.ops, verbose);

                        let s = lodhi.compute(&m.cigar);

                        // let s = m.cost as f64;
                        (s, m, b)
                    })
                    .collect();

                // sort high to low
                scored.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

                // sort low to high
                //scored.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

                if verbose {
                    // print top 5 candidates
                    for (i, (score, m, b)) in scored.iter().take(5).enumerate() {
                        println!("Cost: {:?}", m.cost);
                        println!("Strand: {:?}", m.strand);
                        println!("Cigar: {}", m.cigar.to_string());
                        println!("Score: {}", score);
                        println!("--------------------------------");
                    }
                }

                let (mask_start, mask_end) = barcode_group.bar_region;
                // println!("Mask start: {mask_start}, mask end: {mask_end}");
                // println!(
                //     "bar slice: {}",
                //     String::from_utf8_lossy(&barcode_group.barcodes[0].seq[mask_start..=mask_end])
                // );
                let bar_read_region =
                    map_pat_to_text(&scored[0].1, mask_start as i32, mask_end as i32);
                let ((bar_start, bar_end), (read_bar_start, read_bar_end)) =
                    bar_read_region.expect("No barcode match region found; unusual");

                // 15 and 5 worked well

                let mut is_valid_barcode_match = scored[0].0 > 8.0;
                if scored.len() > 1 {
                    is_valid_barcode_match =
                        is_valid_barcode_match && scored[0].0 - scored[1].0 >= 2.0;
                }
                // } else {
                //     scored[0].0 >= 15.0
                // };

                let first_barcode_cost = scored[0].1.cost;
                let second_barcode_cost = if scored.len() > 1 {
                    scored[1].1.cost
                } else {
                    1000
                };

                // // Single acceptance rule: uniqueness LR
                // let (first_barcode_cost, second_barcode_cost, is_valid_barcode_match) = if scored
                //     .len()
                //     == 1
                // {
                //     let first_barcode_cost = scored[0].1.cost;
                //     if first_barcode_cost <= 6 || scored[0].0 <= 0.6 {
                //         (first_barcode_cost, 1000, true)
                //     } else {
                //         (first_barcode_cost, 1000, false)
                //     }
                // } else {
                //     let odds = (scored[0].0 - scored[1].0);
                //     //println!("odds: {}", odds);

                //     // Then we also check the barcode region to see if we have at most 8 edits
                //     let first_barcode_cost = scored[0].1.cost;
                //     let second_barcode_cost = scored[1].1.cost;

                //     let is_valid_barcode_match =
                //             // Rule 1: penalty <= 9 and gap from 2nd >= 3
                //             (first_barcode_cost <= 9 && (second_barcode_cost - first_barcode_cost) >= 3)
                //             ||
                //             // Rule 2: penalty <= 6 and gap from both 2nd & 3rd >= 6
                //             (second_barcode_cost - first_barcode_cost) >= 6
                //             ||
                //             odds > 1.2;

                //     // let is_valid_barcode_match = first_barcode_cost <= 6
                //     //     || (second_barcode_cost - first_barcode_cost >= 6);

                //     (
                //         first_barcode_cost,
                //         second_barcode_cost,
                //         is_valid_barcode_match,
                //     )
                // };

                if verbose {
                    for _ in 0..10 {
                        println!("Best label: {}", scored[0].2.label);
                        println!("Best match: {}", scored[0].1.cigar.to_string());
                        println!("Best score: {}", scored[0].0);
                        println!("Best cost: {}", first_barcode_cost);

                        if scored.len() > 1 {
                            println!("Second score: {}", scored[1].0);
                            println!("Second label: {}", scored[1].2.label);
                            println!("Second cost: {}", second_barcode_cost);
                            println!("Second cigar: {}", scored[1].1.cigar.to_string());
                            let odds = (scored[1].0 - scored[0].0).exp();
                            println!("Odds: {}", odds);
                        } else {
                            println!("No second match");
                        }
                        println!("Is valid: {}", is_valid_barcode_match);
                        println!("--------------------------------");
                    }
                }

                if is_valid_barcode_match {
                    results.push(BarbellMatch::new(
                        barcode_region_start + read_bar_start,
                        barcode_region_start + read_bar_end,
                        flank_match.text_start,
                        flank_match.text_end,
                        bar_start - barcode_group.flank_prefix.len(),
                        bar_end - barcode_group.flank_prefix.len(),
                        scored[0].2.match_type.clone(),
                        flank_match.cost,
                        scored[0].0.round() as i32,
                        scored[0].2.label.clone(),
                        scored[0].1.strand,
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
                        scored[0].0.round() as i32,
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
        collapse_overlapping_matches(&results, 0.8)
    }
}

// #[cfg(test)]
// mod tests {
//     use pa_types::Cigar;
//
//     use super::*;
//     use crate::annotate::interval::collapse_overlapping_matches;
//
//     /// Helper function to check if a value is within wiggle room of expected
//     fn assert_within_wiggle_room(actual: usize, expected: usize, wiggle_room: usize) {
//         let min = expected.saturating_sub(wiggle_room);
//         let max = expected + wiggle_room;
//         assert!(
//             actual >= min && actual <= max,
//             "Expected {actual} to be within [{min}, {max}] of expected {expected}"
//         );
//     }
// }

/*
5D6=4D4=2D
X3=5D5=3D2=3D
*/
