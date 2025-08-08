use pa_types::CostModel;
use pa_types::Pos;
use pa_types::*;
use sassy::Searcher;
use sassy::*;

pub struct ProbModel {
    /// Cost for one correctly matched base  (−ln (1−ε))
    pub match_cost: f64,
    /// Cost for the *first* base of an error run  (−ln ε  −ln (1−γ))
    pub err_open: f64,
    /// Cost for each additional base in the same error run  (−ln γ)
    pub err_ext: f64,
}

impl ProbModel {
    /// Create a new probability model given
    ///   * `epsilon` = probability that any given base is wrong
    ///   * `gamma`   = probability that an error run continues
    pub fn new(epsilon: f64, gamma: f64) -> Self {
        // Guard against invalid input that would panic on ln(0)
        let eps = epsilon.clamp(1e-6, 1.0 - 1e-6);
        let gam = gamma.clamp(1e-6, 1.0 - 1e-6);

        let match_cost = -((1.0 - eps).ln());
        let err_open = -eps.ln() - (1.0 - gam).ln();
        let err_ext = -gam.ln();

        Self {
            match_cost,
            err_open,
            err_ext,
        }
    }

    pub fn log_odds_ratio(&self, a: f64, b: f64) -> f64 {
        //assert!(a < b);
        (b - a).exp()
    }

    pub fn score_cigar(&self, cigar: &[CigarElem]) -> f64 {
        let mut pending_run: usize = 0; // current error-run length
        let mut total = 0.0_f64;

        let flush = |run_len: &mut usize, tot: &mut f64| {
            if *run_len > 0 {
                *tot += self.err_open + (*run_len as f64 - 1.0) * self.err_ext;
                *run_len = 0;
            }
        };

        for elem in cigar {
            match elem.op {
                CigarOp::Match => {
                    // For matches we may have a run of them; flush errors first
                    flush(&mut pending_run, &mut total);
                    total += elem.cnt as f64 * self.match_cost;
                }
                _ => {
                    // Any non-match contributes to the current error run
                    pending_run += elem.cnt as usize;
                }
            }
        }

        // trailing error run
        flush(&mut pending_run, &mut total);

        total
    }

    /// Convert an error budget m (number of error bases allowed in the barcode region)
    /// into a score cutoff for a region of length L.
    ///
    /// If `conservative=true`, assumes errors are isolated (each opens a run):
    ///   cutoff = L*match_cost + m*err_open
    /// If `conservative=false`, assumes one contiguous run if m>0:
    ///   cutoff = L*match_cost + err_open + (m-1)*err_ext
    pub fn score_cutoff_for_error_budget(&self, l: usize, m: usize, conservative: bool) -> f64 {
        let l_term = l as f64 * self.match_cost;
        if m == 0 {
            return l_term;
        }
        if conservative {
            l_term + (m as f64) * self.err_open
        } else {
            l_term + self.err_open + (m as f64 - 1.0) * self.err_ext
        }
    }

    /// Convert a target "fit" (e.g., 0.8 means allow up to 20% errors) for a region of length L
    /// into a score cutoff under the given run assumption.
    pub fn score_cutoff_for_fit(&self, l: usize, fit: f64, conservative: bool) -> f64 {
        let fit = fit.clamp(0.0, 1.0);
        let allowed_errors = ((1.0 - fit) * l as f64).ceil() as usize;
        self.score_cutoff_for_error_budget(l, allowed_errors, conservative)
    }
}

/// A probability model with separate affine costs for substitutions and indels.
///
/// - `match_cost` applies per matched base.
/// - Substitutions are scored with `sub_open + (len-1) * sub_ext`.
/// - Insertions and deletions are scored with `indel_open + (len-1) * indel_ext`.
pub struct AffineProbModel {
    pub match_cost: f64,
    pub sub_open: f64,
    pub sub_ext: f64,
    pub indel_open: f64,
    pub indel_ext: f64,
}

impl AffineProbModel {
    /// Create a model from separate substitution and indel parameters:
    ///   * `eps_sub`: probability that a substitution run opens
    ///   * `gam_sub`: probability that a substitution run continues
    ///   * `eps_indel`: probability that an indel run opens
    ///   * `gam_indel`: probability that an indel run continues
    ///
    /// The match probability is taken as 1 - eps_sub - eps_indel.
    pub fn new(eps_sub: f64, gam_sub: f64, eps_indel: f64, gam_indel: f64) -> Self {
        // Clamp to avoid ln(0) and invalid probabilities
        let es = eps_sub.clamp(1e-9, 1.0 - 1e-6);
        let ei = eps_indel.clamp(1e-9, 1.0 - 1e-6);
        let gs = gam_sub.clamp(1e-9, 1.0 - 1e-9);
        let gi = gam_indel.clamp(1e-9, 1.0 - 1e-9);
        let p_match = (1.0 - es - ei).max(1e-9);

        let match_cost = -p_match.ln();
        let sub_open = -es.ln() - (1.0 - gs).ln();
        let sub_ext = -gs.ln();
        let indel_open = -ei.ln() - (1.0 - gi).ln();
        let indel_ext = -gi.ln();

        Self {
            match_cost,
            sub_open,
            sub_ext,
            indel_open,
            indel_ext,
        }
    }

    pub fn log_odds_ratio(&self, a: f64, b: f64) -> f64 {
        (b - a).exp()
    }

    /// Score a CIGAR using affine costs per operation class.
    pub fn score_cigar(&self, cigar: &[CigarElem]) -> f64 {
        let mut total = 0.0_f64;
        for elem in cigar {
            match elem.op {
                CigarOp::Match => {
                    total += elem.cnt as f64 * self.match_cost;
                }
                CigarOp::Sub => {
                    let len = elem.cnt.max(1) as f64;
                    total += self.sub_open + (len - 1.0) * self.sub_ext;
                }
                CigarOp::Ins | CigarOp::Del => {
                    let len = elem.cnt.max(1) as f64;
                    total += self.indel_open + (len - 1.0) * self.indel_ext;
                }
            }
        }
        total
    }
}

/// Compute percent identity (as a fraction 0.0..1.0) from a CIGAR.
///
/// Definition used: identity = Matches / (Matches + Substitutions + Insertions + Deletions)
/// where counts are the total base lengths contributed by each op class.
pub fn percent_identity_from_cigar(cigar: &[CigarElem]) -> f64 {
    let mut matches: usize = 0;
    let mut subs: usize = 0;
    let mut ins: usize = 0;
    let mut dels: usize = 0;

    for elem in cigar.iter() {
        match elem.op {
            CigarOp::Match => matches += elem.cnt as usize,
            CigarOp::Sub => subs += elem.cnt as usize,
            CigarOp::Ins => ins += elem.cnt as usize,
            CigarOp::Del => dels += elem.cnt as usize,
        }
    }

    let denom = matches + subs + ins + dels;
    if denom == 0 {
        0.0
    } else {
        matches as f64 / denom as f64
    }
}

/// Map a pattern slice to the corresponding text slice.
/// Returns Some(((pat_start, pat_end), (text_start, text_end)))
/// or None if the pattern slice does not appear in the alignment.
pub fn map_pat_to_text(
    m: &Match,
    p_start: i32,
    p_end: i32,
) -> Option<((usize, usize), (usize, usize))> {
    let mut start_pair: Option<Pos> = None;
    let mut end_pair: Option<Pos> = None;

    // println!("P start: {}, end: {}", p_start, p_end);

    for Pos(i, j) in m.to_path().iter() {
        //  println!("i: {}, j: {}", *i, *j);
        if *i >= p_start && *i < p_end {
            if start_pair.is_none() {
                // println!("Setting start_pair");
                start_pair = Some(Pos(*i, *j));
            }
            end_pair = Some(Pos(*i, *j));
        }
    }

    match (start_pair, end_pair) {
        (Some(Pos(pi, pj)), Some(Pos(ei, ej))) => {
            // convert to half-open ranges
            Some((
                (pi as usize, (ei + 1) as usize),
                (pj as usize, (ej + 1) as usize),
            ))
        }
        _ => None,
    }
}

pub fn extract_cost_at_range_verbose(
    sassy_match: &Match,
    start: usize,
    end: usize,
    p: &[u8],
    t: &[u8],
    alpha: Option<f32>,
) -> Option<i32> {
    // Since overhang might be used we have an incomplete pattern match
    // for which we should use the alpha parameter
    let cost = sassy_match
        .cigar
        .to_path_with_costs(CostModel::unit())
        .iter()
        .skip(1) // skip
        .map(|x| x.1)
        .collect::<Vec<_>>();

    // println!("Cost: {:?}", cost);
    let path = sassy_match.to_path();
    //let cigar = sassy_match.cigar.clone();

    // REMOVE
    let start = 0;
    let end = p.len();

    let mut cost_in_region = vec![];
    let mut last_cost: Option<i32> = None;
    for (Pos(q_pos, r_pos), cost) in path.iter().zip(cost.iter()) {
        let mut marker = "";
        let q_pos = *q_pos as usize;
        if q_pos >= start && q_pos <= end {
            // Count when we transition to a different cost (indicating a new operation)
            if last_cost.is_some() && *cost != last_cost.unwrap() {
                cost_in_region.push(1);
            } else {
                cost_in_region.push(0)
            }
            marker = "*";
        }

        // println!(
        //     "q_pos: {}:{}, r_pos: {}:{} - cost: {} {} (last_cost: {:?})",
        //     q_pos, p[q_pos] as char, r_pos, t[*r_pos as usize] as char, cost, marker, last_cost
        // );

        last_cost = Some(*cost);
    }
    let left_overhang = sassy_match.pattern_start;
    // println!("Cost vector: {:?}", cost_in_region);
    // If any in range we have to add it using overhang cost
    let mut summed_cost = cost_in_region.iter().sum::<i32>();

    /*
    New cost
     */

    let n = cost_in_region.len();

    // Count transitions 0→1
    let mut transitions = cost_in_region
        .windows(2)
        .filter(|w| w[0] == 0 && w[1] == 1)
        .count();

    // Count total matches X
    let matches = cost_in_region.iter().filter(|&&b| b == 0).count();

    // Compute integer cost: ceil(transitions * n / matches)
    let mut trans_cost = if matches == 0 {
        // no matches → worst cost = n
        n
    } else {
        (transitions * n + matches - 1) / matches
    };

    if start < left_overhang {
        let overhang = left_overhang - start;
        let l_overhang_cost = (overhang as f32 * alpha.unwrap_or(0.0)).ceil() as i32;
        //  println!("Adding left overhang: {}", l_overhang_cost);
        summed_cost += l_overhang_cost;
        trans_cost += l_overhang_cost as usize;
    }

    if end > sassy_match.pattern_end {
        let right_overhang = (end + 1) - sassy_match.pattern_end;
        assert!(end < p.len());
        let r_overhang_cost = (right_overhang as f32 * alpha.unwrap_or(0.0)).ceil() as i32;
        //println!("Adding right overhang: {}", r_overhang_cost);
        summed_cost += r_overhang_cost;
        trans_cost += r_overhang_cost as usize;
    }

    // println!(
    //     "Transitions: {} and cost: {} with cost vector: {:?}",
    //     summed_cost, summed_cost, cost_in_region
    // );
    // println!("Summed cost: {}", );
    Some(trans_cost as i32)
}

pub fn extract_cost_at_range(
    sassy_match: &Match,
    start: usize,
    end: usize,
    alpha: Option<f32>,
) -> Option<i32> {
    let alpha = alpha.unwrap_or(0.0);

    // Count edits in the aligned region
    let mut region_cost: i32 = 0;
    let mut last_cost: Option<i32> = None;

    for (Pos(q_pos, _), cost) in sassy_match.to_path().iter().zip(
        sassy_match
            .cigar
            .to_path_with_costs(CostModel::unit())
            .iter()
            .skip(1) // skip the first position
            .map(|x| x.1),
    ) {
        let q_pos = *q_pos as usize;
        // Track cost changes from the beginning, not just in the window
        if q_pos >= start && q_pos <= end {
            // Count when we transition to a different cost (indicating a new operation)
            if last_cost.is_some() && cost != last_cost.unwrap() {
                region_cost += 1;
            }
        }
        last_cost = Some(cost);
    }

    // Add overhang penalties
    let left_overhang = sassy_match.pattern_start.saturating_sub(start);
    let right_overhang = (end + 1).saturating_sub(sassy_match.pattern_end);

    let overhang_cost = ((left_overhang + right_overhang) as f32 * alpha).ceil() as i32;

    Some(region_cost + overhang_cost)
}

/// Produce an owned CIGAR subregion covering the query half-open window [q_start, q_end).
/// This will split CIGAR elements if the boundaries fall inside an element.
pub fn subcigar_owned(cigar: &Cigar, q_start: usize, q_end: usize) -> Vec<CigarElem> {
    let q_start = q_start as i32;
    let q_end = q_end as i32;
    let mut i = 0;
    let mut j = 0;
    let mut expanded_cigar = vec![];
    for op in cigar.ops.iter() {
        for _ in 0..op.cnt {
            if i >= q_start && i < q_end {
                expanded_cigar.push(op.op);
            }
            match op.op {
                CigarOp::Match => {
                    i += 1;
                    j += 1;
                }
                CigarOp::Sub => {
                    i += 1;
                    j += 1;
                }
                CigarOp::Ins => {
                    j += 1;
                }
                CigarOp::Del => {
                    i += 1;
                }
            }
        }
    }
    collapse_ops_to_elems(&expanded_cigar)
}

/// Collapse a sequence of single-base CigarOps into run-length encoded CIGAR elements.
pub fn collapse_ops_to_elems(ops: &[CigarOp]) -> Vec<CigarElem> {
    let mut out: Vec<CigarElem> = Vec::new();
    for op in ops.iter() {
        if let Some(last) = out.last_mut() {
            if last.op == *op {
                last.cnt += 1;
                continue;
            }
        }
        out.push(CigarElem { op: *op, cnt: 1 });
    }
    out
}

/// Merge adjacent identical CIGAR elements by summing their counts.
pub fn collapse_cigar_elems(elems: &[CigarElem]) -> Vec<CigarElem> {
    let mut out: Vec<CigarElem> = Vec::new();
    for e in elems.iter() {
        if let Some(last) = out.last_mut() {
            if last.op == e.op {
                last.cnt += e.cnt;
                continue;
            }
        }
        out.push(*e);
    }
    out
}

#[cfg(test)]
mod subcigar_tests {
    use super::*;

    // #[test]
    // fn test_subcigar_trim_match_and_sub() {
    //     use pa_types::CigarOp::*;
    //     let ops = [
    //         CigarElem { op: Match, cnt: 3 }, // q:0..3
    //         CigarElem { op: Sub, cnt: 2 },   // q:3..5
    //         CigarElem { op: Ins, cnt: 1 },   // q:5..6
    //         CigarElem { op: Match, cnt: 4 }, // q:6..10
    //     ];
    //     let sub = subcigar(&ops, 2, 4); // overlaps last base of M3 and first of Sub2
    //     assert_eq!(
    //         sub,
    //         vec![
    //             CigarElem { op: Match, cnt: 1 },
    //             CigarElem { op: Sub, cnt: 1 },
    //         ]
    //     );
    // }

    // #[test]
    // fn test_subcigar_includes_deletion_if_anchored_in_window() {
    //     use pa_types::CigarOp::*;
    //     let ops = [
    //         CigarElem { op: Match, cnt: 3 }, // q:0..3
    //         CigarElem { op: Del, cnt: 2 },   // anchored at q:3
    //         CigarElem { op: Match, cnt: 2 }, // q:3..5
    //     ];
    //     let sub = subcigar(&ops, 3, 5); // window starts at q_pos 3 (deletion anchor)
    //     assert_eq!(
    //         sub,
    //         vec![
    //             CigarElem { op: Del, cnt: 2 },
    //             CigarElem { op: Match, cnt: 2 },
    //         ]
    //     );
    // }
}

#[cfg(test)]
mod collapse_tests {
    use super::*;
    use pa_types::CigarOp::*;

    #[test]
    fn test_collapse_ops_to_elems() {
        let ops = vec![Match, Match, Sub, Sub, Sub, Ins, Del, Del, Match];
        let got = collapse_ops_to_elems(&ops);
        let expect = vec![
            CigarElem { op: Match, cnt: 2 },
            CigarElem { op: Sub, cnt: 3 },
            CigarElem { op: Ins, cnt: 1 },
            CigarElem { op: Del, cnt: 2 },
            CigarElem { op: Match, cnt: 1 },
        ];
        assert_eq!(got, expect);
    }

    #[test]
    fn test_collapse_cigar_elems() {
        let elems = vec![
            CigarElem { op: Match, cnt: 2 },
            CigarElem { op: Match, cnt: 3 },
            CigarElem { op: Sub, cnt: 1 },
            CigarElem { op: Sub, cnt: 2 },
            CigarElem { op: Del, cnt: 1 },
            CigarElem { op: Del, cnt: 1 },
        ];
        let got = collapse_cigar_elems(&elems);
        let expect = vec![
            CigarElem { op: Match, cnt: 5 },
            CigarElem { op: Sub, cnt: 3 },
            CigarElem { op: Del, cnt: 2 },
        ];
        assert_eq!(got, expect);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use sassy::Searcher;
    use sassy::profiles::Iupac;

    fn rc(seq: &[u8]) -> Vec<u8> {
        let mut rc = seq.to_vec();
        rc.reverse();
        for b in rc.iter_mut() {
            *b = match *b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                _ => *b,
            };
        }
        rc
    }

    #[test]
    fn test_extract_cost_at_range() {
        let overhang = 0.5;
        let mut searcher = Searcher::<Iupac>::new_rc_with_overhang(overhang);

        let p = b"AAACACAAA".to_vec();
        let t = b"GGGGAAATATA".to_vec();
        let t = rc(&t);
        let matches = searcher.search(&p, &t, 3);

        println!("Pattern len: {}", p.len());

        let m = matches.first().unwrap();
        println!("match: {:?}", m);

        let cost_5_6 = extract_cost_at_range(m, 5, 6, Some(0.5));
        assert_eq!(cost_5_6, Some(2));

        let cost_5_6 = extract_cost_at_range(m, 5, 8, Some(0.5));
        // 2 more positions than before but with 0.5 cost, 2 * 0.5 = +1
        assert_eq!(cost_5_6, Some(3));

        let cost_7_8 = extract_cost_at_range(m, 7, 8, Some(0.5));
        assert_eq!(cost_7_8, Some(1));
    }

    #[test]
    fn test_prob_model() {
        let p = ProbModel::new(0.01, 0.15);
        let read = b"TGTATGTATACCTGAAATTCGGCTATGCCGACCGATCTGGCGCGCTTCGTAAGTCCTTGTTTTCGCATTTATCGTGAAACGCTTTCGCATTTTCGTGCGCCGCTTCAATGGAAATCCCCGGCATCGGCACGATCACGGCGCTGTCCTTCTATAGCGCGATGAGTGGACCCGACGCGGTTTACAGGCATGTCGACGATGTCGCGGCCTATCTGGGCCTGACCCCGCGCGTCTATCAGTCGGGCGAGAGCCTGACACATGGCGGGATCAGCAAGATGGGCAACCAGTTGACCCGCACGCATC";
        let bar = b"GCTTGGGTGTTTAACCTCTGCCACACACTCGTAAGTCCTTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";

        let mut searcher = Searcher::<Iupac>::new_rc_with_overhang(0.5);
        let matches = searcher.search(bar, &read, 27);
        // Get match with lowest cost
        let best_match = matches.iter().min_by_key(|m| m.cost).unwrap();

        let score = p.score_cigar(&best_match.cigar.ops);
        println!("Score: {}", score);
    }

    #[test]
    fn test_log_odds() {
        let model = ProbModel::new(0.01, 0.15);
        let a = 60.0;
        let b = 80.0;
        let odds = model.log_odds_ratio(a, b);
        println!("Odds: {}", odds);
        assert!(odds > 10000.0);
    }

    #[test]
    fn test_bug() {
        use pa_types::CigarOp::*;
        let ops = [
            CigarElem { op: Sub, cnt: 2 },
            CigarElem { op: Match, cnt: 1 },
            CigarElem { op: Sub, cnt: 2 },
            CigarElem { op: Match, cnt: 1 },
            CigarElem { op: Del, cnt: 1 },
            CigarElem { op: Sub, cnt: 1 },
            CigarElem { op: Match, cnt: 1 },
            CigarElem { op: Sub, cnt: 1 },
            CigarElem { op: Match, cnt: 1 },
            CigarElem { op: Sub, cnt: 1 },
            CigarElem { op: Match, cnt: 1 },
            CigarElem { op: Del, cnt: 2 },
            CigarElem { op: Match, cnt: 1 },
            CigarElem { op: Del, cnt: 1 },
            CigarElem { op: Match, cnt: 2 },
            CigarElem { op: Sub, cnt: 1 },
            CigarElem { op: Match, cnt: 1 },
            CigarElem { op: Sub, cnt: 1 },
            CigarElem { op: Match, cnt: 56 },
            CigarElem { op: Del, cnt: 4 },
            CigarElem { op: Match, cnt: 7 },
        ];
        let cigar = Cigar { ops: ops.to_vec() };
        println!("Cigar: {:?}", cigar.to_string());
        let score = ProbModel::new(0.01, 0.15).score_cigar(&cigar.ops);
        println!("Score: {}", score);
    }

    #[test]
    fn test_mapping_substring_in_cigar() {
        let p = b"AAAAACCCAAAA";
        let t = b"GGGGAAAAACCCAAAAGGGGG";
        let mut searcher = Searcher::<Iupac>::new_rc();
        let matches = searcher.search(p, &t, 0);
        let m = matches.first().unwrap();
        let ((ps, pe), (t_start, t_end)) = map_pat_to_text(m, 5, 7 + 1).unwrap();
        let match_slice = &t[t_start..t_end];
        println!(
            "Text slice: {:?}",
            String::from_utf8_lossy(&t[t_start..t_end])
        );
        assert_eq!(ps, 5);
        assert_eq!(pe, 8);
        assert!(match_slice == b"CCC");
    }

    #[test]
    fn test_mapping_substring_in_cigar_left_overhang() {
        let p = b"TTTCCCAAAA";
        let t = b"CAAAAGGGGG";
        let mut searcher = Searcher::<Iupac>::new_rc_with_overhang(0.5);
        let matches = searcher.search(p, &t, 3);
        let m = matches.first().unwrap();
        let ((ps, pe), (t_start, t_end)) = map_pat_to_text(m, 3, 5 + 1).unwrap();
        let match_slice = &t[t_start..t_end];
        println!(
            "Text slice: ({}, {}) - {:?}",
            t_start,
            t_end,
            String::from_utf8_lossy(&t[t_start..t_end])
        );
        assert_eq!(ps, 5);
        assert_eq!(pe, 6);
        assert!(match_slice == b"C");
    }

    #[test]
    fn test_get_fit() {
        let l = 24;
        let fit = 0.7;
        let cutoff = ProbModel::new(0.01, 0.15).score_cutoff_for_fit(l, fit, false);
        println!("Cutoff: {}", cutoff);

        let cutoff = ProbModel::new(0.01, 0.15).score_cutoff_for_fit(l, fit, true);
        println!("Cutoff: {}", cutoff);
    }
}
