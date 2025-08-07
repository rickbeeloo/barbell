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
}
