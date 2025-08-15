use cigar_lodhi_rs::Lodhi;
use pa_types::CostModel;
use pa_types::Pos;
use pa_types::*;
use sassy::Searcher;
use sassy::*;

#[derive(Copy, Clone, Debug)]
pub struct EmaParams {
    pub match_step: f64,       // raw +1 for a match (default 1.0)
    pub error_step: f64,       // raw -1 for an error (default -1.0)
    pub err_open: f64,         // extra penalty on first error of a run (>=0)
    pub err_ext_coef: f64,     // extra penalty added per additional error in run
    pub ema_alpha: f64,        // alpha in EMA (0..1). larger => more memory
    pub init_with_first: bool, // if true initialize ema with first raw step instead of 0.0
}

pub fn ema_score(cigar: &[CigarElem], params: &EmaParams) -> f64 {
    // build raw steps per base while tracking error runs
    let mut raw_steps: Vec<f64> = Vec::new();
    let mut match_run = 0usize;
    let mut err_run = 0usize;

    for elem in cigar {
        match elem.op {
            CigarOp::Match => {
                for _ in 0..elem.cnt {
                    match_run += 1;
                    err_run = 0;
                    raw_steps.push(params.match_step);
                }
            }
            CigarOp::Sub | CigarOp::Del | CigarOp::Ins => {
                for _ in 0..elem.cnt {
                    match_run = 0;
                    err_run += 1;
                    let step = if err_run == 1 {
                        params.error_step - params.err_open
                    } else {
                        // progressively harsher with run length: error_step - err_open - ext_coef*(run_len-1)
                        params.error_step
                            - params.err_open
                            - params.err_ext_coef * ((err_run - 1) as f64)
                    };
                    raw_steps.push(step);
                }
            }
        }
    }

    if raw_steps.is_empty() {
        return 0.0;
    }

    // apply EMA over raw_steps and accumulate the EMA as increments into score
    let mut score = 0.0_f64;
    let mut ema = if params.init_with_first {
        raw_steps[0]
    } else {
        0.0
    };

    for &s in raw_steps.iter() {
        ema = params.ema_alpha * ema + (1.0 - params.ema_alpha) * s;
        score += ema;
    }
    score
}

fn pattern_bases_in_cigar(cigar: &[CigarElem]) -> usize {
    let mut len: usize = 0;
    for e in cigar.iter() {
        match e.op {
            CigarOp::Match | CigarOp::Sub | CigarOp::Del => len += e.cnt as usize,
            CigarOp::Ins => {}
        }
    }
    len
}

pub fn triangle_score_both(cigar: &[CigarElem], match_step: f64, error_step: f64) -> f64 {
    if cigar.is_empty() {
        return 0.0;
    }

    let mut match_total: f64 = 0.0;
    let mut error_total: f64 = 0.0;
    let mut cur_match_run: usize = 0;
    let mut cur_error_run: usize = 0;

    for elem in cigar {
        match elem.op {
            CigarOp::Match => {
                // If there was an error run, flush it
                if cur_error_run > 0 {
                    let l = cur_error_run as f64;
                    error_total += error_step * (l * (l + 1.0) / 2.0);
                    cur_error_run = 0;
                }
                // extend current match run by elem.cnt
                cur_match_run += elem.cnt as usize;
            }

            // treat Sub, Ins, Del as errors; flush any pending match run first
            CigarOp::Sub | CigarOp::Ins | CigarOp::Del => {
                if cur_match_run > 0 {
                    let l = cur_match_run as f64;
                    match_total += match_step * (l * (l + 1.0) / 2.0);
                    cur_match_run = 0;
                }
                // extend current error run by elem.cnt
                cur_error_run += elem.cnt as usize;
            }
        }
    }

    // flush at end
    if cur_match_run > 0 {
        let l = cur_match_run as f64;
        match_total += match_step * (l * (l + 1.0) / 2.0);
    }
    if cur_error_run > 0 {
        let l = cur_error_run as f64;
        error_total += error_step * (l * (l + 1.0) / 2.0);
    }

    let raw_score = match_total - error_total;

    let qlen = pattern_bases_in_cigar(cigar);

    if qlen == 0 {
        unreachable!("CIGAR should at least have some match");
    } else {
        raw_score / (qlen as f64)
    }
}

pub fn triangle_score(cigar: &[CigarElem], match_step: f64, error_step: f64) -> f64 {
    if cigar.is_empty() {
        return 0.0;
    }

    let mut match_total: f64 = 0.0;
    let mut error_total: f64 = 0.0;
    let mut cur_match_run: usize = 0;

    for elem in cigar {
        match elem.op {
            CigarOp::Match => {
                // extend current match run by elem.cnt
                cur_match_run += elem.cnt as usize;
            }

            // treat Sub, Ins, Del as errors; flush any pending match run first
            CigarOp::Sub | CigarOp::Ins | CigarOp::Del => {
                if cur_match_run > 0 {
                    let l = cur_match_run as f64;
                    match_total += match_step * (l * (l + 1.0) / 2.0);
                    cur_match_run = 0;
                }
                // errors contribute per-base, we could do "streaks" here as well
                // as streaks of errors are likely more "okay" than single errors
                error_total += (elem.cnt as f64) * error_step;
            }
        }
    }

    // flush at end
    if cur_match_run > 0 {
        let l = cur_match_run as f64;
        match_total += match_step * (l * (l + 1.0) / 2.0);
    }

    let raw_score = match_total + error_total;

    let qlen = pattern_bases_in_cigar(cigar);

    if qlen == 0 {
        unreachable!("CIGAR should at least have some match");
    } else {
        raw_score / (qlen as f64)
    }
}

pub fn simple_kmer_score(cigar: &[CigarElem], match_step: f64, error_step: f64, k: usize) -> f64 {
    if cigar.is_empty() {
        return 0.0;
    }

    //const K: usize = 6;

    // total number of pattern positions spanned by the CIGAR (M+S+D)
    let pattern_len = pattern_bases_in_cigar(cigar);
    if pattern_len == 0 {
        return 0.0;
    }

    // Use absolute scales to enforce: matches add, errors subtract
    let m_scale = match_step.abs();
    let e_scale = error_step.abs();

    let kernel_sum_at = |anchor: usize| -> f64 {
        let left_span = anchor.min(k - 1);
        let right_span = (pattern_len - 1 - anchor).min(k - 1);
        let mut sum = k as f64; // d=0
        for d in 1..=left_span {
            sum += (k - d) as f64;
        }
        for d in 1..=right_span {
            sum += (k - d) as f64;
        }
        sum
    };

    let mut total_match: f64 = 0.0;
    let mut total_error: f64 = 0.0;

    let mut i_pat: usize = 0;

    for elem in cigar {
        match elem.op {
            CigarOp::Match => {
                for _ in 0..elem.cnt {
                    let anchor = i_pat;
                    let ksum = kernel_sum_at(anchor);
                    total_match += m_scale * ksum;
                    i_pat += 1;
                }
            }
            CigarOp::Sub => {
                for _ in 0..elem.cnt {
                    let anchor = i_pat;
                    let ksum = kernel_sum_at(anchor);
                    total_error += e_scale * ksum;
                    i_pat += 1;
                }
            }
            CigarOp::Del => {
                for _ in 0..elem.cnt {
                    let anchor = i_pat;
                    let ksum = kernel_sum_at(anchor);
                    total_error += e_scale * ksum;
                    i_pat += 1;
                }
            }
            CigarOp::Ins => {
                for _ in 0..elem.cnt {
                    let anchor = i_pat.min(pattern_len - 1);
                    let ksum = kernel_sum_at(anchor);
                    total_error += e_scale * ksum;
                }
            }
        }
    }

    let total_ops = cigar.iter().map(|e| e.cnt as usize).sum::<usize>();

    (total_match - total_error) / (total_ops as f64)
}

// pub fn ema_score_simple(cigar: &[CigarElem], match_step: f64, error_step: f64, alpha: f64) -> f64 {
//     if cigar.is_empty() {
//         return 0.0;
//     }

//     let mut ema = 0.0_f64;
//     let mut score = 0.0_f64;

//     for elem in cigar {
//         let step_val = match elem.op {
//             CigarOp::Match => match_step,
//             CigarOp::Sub | CigarOp::Del | CigarOp::Ins => error_step,
//         };
//         for _ in 0..elem.cnt {
//             ema = alpha * ema + (1.0 - alpha) * step_val;
//             score += ema;
//         }
//     }

//     score
// }

pub fn ema_area_between(a: &[CigarElem], b: &[CigarElem], params: &EmaParams) -> f64 {
    let curve_a = ema_curve(a, params);
    let curve_b = ema_curve(b, params);
    curve_a
        .iter()
        .zip(curve_b.iter())
        .map(|(x, y)| (x - y).abs())
        .sum()
}

/// Same streaming EMA but returns the whole curve for each base.
fn ema_curve(cigar: &[CigarElem], params: &EmaParams) -> Vec<f64> {
    let mut ema_vals = Vec::new();
    let mut ema = 0.0_f64;
    let mut first = true;
    let mut match_run = 0;
    let mut err_run = 0;

    for elem in cigar {
        match elem.op {
            CigarOp::Match => {
                for _ in 0..elem.cnt {
                    match_run += 1;
                    err_run = 0;
                    let raw = params.match_step;
                    if first {
                        ema = if params.init_with_first { raw } else { 0.0 };
                        first = false;
                    }
                    ema = params.ema_alpha * ema + (1.0 - params.ema_alpha) * raw;
                    ema_vals.push(ema);
                }
            }
            CigarOp::Sub | CigarOp::Del | CigarOp::Ins => {
                for _ in 0..elem.cnt {
                    match_run = 0;
                    err_run += 1;
                    let raw = if err_run == 1 {
                        params.error_step - params.err_open
                    } else {
                        params.error_step
                            - params.err_open
                            - params.err_ext_coef * ((err_run - 1) as f64)
                    };
                    if first {
                        ema = if params.init_with_first { raw } else { 0.0 };
                        first = false;
                    }
                    ema = params.ema_alpha * ema + (1.0 - params.ema_alpha) * raw;
                    ema_vals.push(ema);
                }
            }
        }
    }
    ema_vals
}

/// Compute difference score(a) - score(b)
pub fn ema_score_diff(a: &[CigarElem], b: &[CigarElem], params: &EmaParams) -> f64 {
    ema_score(a, params) - ema_score(b, params)
}

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

    pub fn score_cigar_errors_only(&self, cigar: &[CigarElem]) -> f64 {
        let mut pending_run: usize = 0;
        let mut total = 0.0;

        let flush = |run_len: &mut usize, tot: &mut f64| {
            if *run_len > 0 {
                *tot += self.err_open + (*run_len as f64 - 1.0) * self.err_ext;
                *run_len = 0;
            }
        };

        for elem in cigar {
            match elem.op {
                CigarOp::Match => flush(&mut pending_run, &mut total),
                _ => pending_run += elem.cnt as usize,
            }
        }

        flush(&mut pending_run, &mut total);
        total
    }

    pub fn score_cigar(&self, cigar: &[CigarElem], verbose: bool) -> f64 {
        let mut pending_run: usize = 0; // current error-run length
        let mut total = 0.0_f64;

        let flush = |run_len: &mut usize, tot: &mut f64| {
            if verbose {
                println!("Flushing error run of length {}", *run_len);
            }

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

    /// Score a CIGAR like `score_cigar` but subtract a small reward for contiguous match runs.
    /// The reward is affine per run: `match_open_bonus + (run_len-1) * match_ext_bonus`.
    /// Set both bonuses to zero to recover the original score.
    pub fn score_cigar_with_match_bonus(
        &self,
        cigar: &[CigarElem],
        verbose: bool,
        match_ext_bonus: f64,
    ) -> f64 {
        let base = self.score_cigar(cigar, verbose) - match_ext_bonus;
        if match_ext_bonus == 0.0 {
            return base;
        }
        // Sum rewards over contiguous Match runs (elements are already RLE-encoded)
        let mut bonus_total = 0.0_f64;
        for elem in cigar.iter() {
            if let CigarOp::Match = elem.op {
                let len = elem.cnt.max(1) as f64;
                bonus_total += (len - 1.0) * match_ext_bonus;
            }
        }
        base - bonus_total
    }

    /// Count how many pattern positions are spanned by this CIGAR: Match + Sub + Del.
    /// Insertions do not advance the pattern.
    fn pattern_bases_in_cigar(cigar: &[CigarElem]) -> usize {
        let mut len: usize = 0;
        for e in cigar.iter() {
            match e.op {
                CigarOp::Match | CigarOp::Sub | CigarOp::Del => len += e.cnt as usize,
                CigarOp::Ins => {}
            }
        }
        len
    }

    /// Average negative log-likelihood per pattern base (nats/base).
    /// Useful for using a single cutoff independent of query length.
    pub fn avg_cost_per_pattern_base(&self, cigar: &[CigarElem], verbose: bool) -> f64 {
        let total_len = Self::pattern_bases_in_cigar(cigar);
        if total_len == 0 {
            return 0.0;
        }
        if verbose {
            println!("Cigar: {:?}", cigar);
        }
        self.score_cigar(cigar, verbose) / total_len as f64
    }

    /// Per-base average using match-run bonuses.
    pub fn avg_cost_per_pattern_base_with_match_bonus(
        &self,
        cigar: &[CigarElem],
        verbose: bool,
        match_ext_bonus: f64,
    ) -> f64 {
        let total_len = Self::pattern_bases_in_cigar(cigar);
        if total_len == 0 {
            return 0.0;
        }
        self.score_cigar_with_match_bonus(cigar, verbose, match_ext_bonus) / total_len as f64
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

    /// Count how many pattern positions are spanned by this CIGAR: Match + Sub + Del.
    /// Insertions do not advance the pattern.
    fn pattern_bases_in_cigar(cigar: &[CigarElem]) -> usize {
        let mut len: usize = 0;
        for e in cigar.iter() {
            match e.op {
                CigarOp::Match | CigarOp::Sub | CigarOp::Del => len += e.cnt as usize,
                CigarOp::Ins => {}
            }
        }
        len
    }

    /// Average negative log-likelihood per pattern base (nats/base).
    /// Useful for using a single cutoff independent of query length.
    pub fn avg_cost_per_pattern_base(&self, cigar: &[CigarElem]) -> f64 {
        let total_len = Self::pattern_bases_in_cigar(cigar);
        if total_len == 0 {
            return 0.0;
        }

        self.score_cigar(cigar) / total_len as f64
    }
}

/// BLAST-style affine gap scoring model with bit-score computation.
pub struct BlastModel {
    pub match_reward: f64,
    pub mismatch_penalty: f64,
    pub gap_open_penalty: f64,
    pub gap_extend_penalty: f64,
    pub lambda: f64,
    pub k: f64,
    pub h: f64, // entropy parameter H
}

impl Default for BlastModel {
    fn default() -> Self {
        Self {
            match_reward: 2.0,
            mismatch_penalty: 3.0,
            gap_open_penalty: 5.0,
            gap_extend_penalty: 2.0,
            lambda: 1.37406,
            k: 0.710603,
            h: 1.30725,
        }
    }
}

impl BlastModel {
    pub fn new(
        match_reward: f64,
        mismatch_penalty: f64,
        gap_open_penalty: f64,
        gap_extend_penalty: f64,
        lambda: f64,
        k: f64,
        h: f64,
    ) -> Self {
        Self {
            match_reward,
            mismatch_penalty,
            gap_open_penalty,
            gap_extend_penalty,
            lambda,
            k,
            h,
        }
    }

    pub fn raw_score_cigar(&self, cigar: &[CigarElem]) -> f64 {
        let mut score: f64 = 0.0;
        for elem in cigar.iter() {
            match elem.op {
                CigarOp::Match => {
                    score += (elem.cnt as f64) * self.match_reward;
                }
                CigarOp::Sub => {
                    score -= (elem.cnt as f64) * self.mismatch_penalty;
                }
                CigarOp::Ins | CigarOp::Del => {
                    let len = elem.cnt.max(1) as f64;
                    score -= self.gap_open_penalty + (len - 1.0) * self.gap_extend_penalty;
                }
            }
        }
        score
    }

    pub fn bit_score(&self, raw_score: i32) -> f64 {
        let s = raw_score as f64;
        (self.lambda * s - self.k.ln()) / std::f64::consts::LN_2
    }

    /// Convenience: compute bit score directly from a CIGAR.
    pub fn bit_score_from_cigar(&self, cigar: &[CigarElem]) -> f64 {
        let raw = self.raw_score_cigar(cigar);
        self.bit_score(raw as i32)
    }

    /// Convenience: compute E-value directly from a CIGAR with given query (m) and subject (n) lengths.
    pub fn e_value_from_cigar(&self, cigar: &[CigarElem], m: usize, n: usize) -> f64 {
        let bits = self.bit_score_from_cigar(cigar);
        self.e_value(bits, m, n)
    }

    /// Solve length adjustment l by fixed-point iteration:
    /// l = (1/H) * ln(K * (m - l) * (n - l))
    pub fn calc_length_adjustment(&self, m: usize, n: usize) -> f64 {
        let m = m as f64;
        let n = n as f64;

        // BLAST initial guess
        let mut l = (self.k * m * n).ln() / self.h;

        let max_iter = 100usize;
        let eps = 1e-5_f64; // BLAST uses 1e-5

        for _ in 0..max_iter {
            let prev_l = l;

            // BLAST uses max(m - l, 1.0) and max(n - l, 1.0) INSIDE the iteration
            let m_eff = (m - prev_l).max(1.0);
            let n_eff = (n - prev_l).max(1.0);

            // compute new l
            l = (self.k * m_eff * n_eff).ln() / self.h;

            // ensure non-negative (BLAST behavior keeps it non-negative)
            if l < 0.0 {
                l = 0.0;
            }

            if (l - prev_l).abs() < eps {
                break;
            }
        }

        // final clamp: cannot exceed the smallest length and cannot be negative
        l.max(0.0).min(m.min(n))
    }
    /// Calculate effective length given raw length and adjustment
    pub fn effective_length(&self, raw_len: usize, length_adjustment: f64) -> f64 {
        (raw_len as f64 - length_adjustment).max(1.0) // minimum 1 to avoid zero or negative
    }

    /// Compute E-value using effective lengths
    pub fn e_value(&self, bit_score: f64, m: usize, n: usize) -> f64 {
        let l = self.calc_length_adjustment(m, n);
        let m_eff = self.effective_length(m, l);
        let n_eff = self.effective_length(n, l);
        //  println!("m eff: {}, n eff: {}", m_eff, n_eff);
        m_eff * n_eff * 2.0_f64.powf(-bit_score)
    }

    fn score_elem(&self, elem: &CigarElem) -> f64 {
        match elem.op {
            CigarOp::Match => (elem.cnt as f64) * self.match_reward,
            CigarOp::Sub => -((elem.cnt as f64) * self.mismatch_penalty),
            CigarOp::Ins | CigarOp::Del => {
                let len = elem.cnt.max(1) as f64;
                -(self.gap_open_penalty + (len - 1.0) * self.gap_extend_penalty)
            }
        }
    }

    pub fn top_vs_second_is_x_fold(&self, s1: i32, s2: i32, x: f64) -> (bool, bool, f64, f64) {
        assert!(x > 0.0, "x must be positive");

        // raw score difference
        let delta_raw = (s1 - s2) as f64;

        // required delta in raw scoring units: ln(x) / lambda
        let needed_delta_raw = x.ln() / self.lambda;

        // bits difference (for interpretability)
        let delta_bits = (self.lambda * delta_raw) / std::f64::consts::LN_2;
        let needed_bits = x.log2();

        // two equivalent tests (use whichever you prefer)
        let is_significant_raw = delta_raw >= needed_delta_raw;
        let is_significant_bits = delta_bits >= needed_bits;

        (
            is_significant_raw,
            is_significant_bits,
            delta_bits,
            delta_raw,
        )
    }

    // /// Element-granularity local sub-cigar (Kadane / Smith-Waterman style).
    // /// If `debug` is true, prints per-element scores and Kadane state transitions.
    // pub fn get_local_element_level(&self, cigar: &[CigarElem], debug: bool) -> Vec<CigarElem> {
    //     if cigar.is_empty() {
    //         return Vec::new();
    //     }

    //     if debug {
    //         println!("Index  Op   Cnt   Score");
    //         println!("------------------------");
    //         for (i, e) in cigar.iter().enumerate() {
    //             println!(
    //                 "{:5}  {:5?}  {:3}  {:4}",
    //                 i,
    //                 e.op,
    //                 e.cnt,
    //                 self.score_elem(e)
    //             );
    //         }
    //         println!("------------------------");
    //     }

    //     let mut current_sum: i32 = 0;
    //     let mut current_start: usize = 0;

    //     let mut best_sum: i32 = 0; // Smith-Waterman style: only positive maxima
    //     let mut best_start: usize = 0;
    //     let mut best_end: usize = 0; // exclusive

    //     for (i, elem) in cigar.iter().enumerate() {
    //         let score = self.score_elem(elem);

    //         if debug {
    //             println!("i={} score={} current_sum_before={}", i, score, current_sum);
    //         }

    //         if current_sum + score <= 0 {
    //             // reset after this element
    //             if debug {
    //                 println!("  RESET (current_sum + score <= 0). Next start = {}", i + 1);
    //             }
    //             current_sum = 0;
    //             current_start = i + 1;
    //         } else {
    //             current_sum += score;
    //             if debug {
    //                 println!("  EXTEND: current_sum_after={}", current_sum);
    //             }
    //             if current_sum > best_sum {
    //                 best_sum = current_sum;
    //                 best_start = current_start;
    //                 best_end = i + 1;
    //                 if debug {
    //                     println!(
    //                         "  NEW BEST: sum={} start={} end={}",
    //                         best_sum, best_start, best_end
    //                     );
    //                 }
    //             }
    //         }
    //     }

    //     if debug {
    //         println!(
    //             "final best_sum={} start={} end={}",
    //             best_sum, best_start, best_end
    //         );
    //     }

    //     if best_sum <= 0 {
    //         return Vec::new();
    //     }

    //     // Return compacted sub-cigar (clone)
    //     let sub = cigar[best_start..best_end].to_vec();
    //     self.compact_adjacent_same_op(sub)
    // }

    /// Merge adjacent runs of the same operation for tidier output.
    fn compact_adjacent_same_op(&self, mut v: Vec<CigarElem>) -> Vec<CigarElem> {
        if v.is_empty() {
            return v;
        }
        let mut out: Vec<CigarElem> = Vec::with_capacity(v.len());
        let mut cur = v.remove(0);
        for e in v.into_iter() {
            if e.op == cur.op {
                cur.cnt += e.cnt;
            } else {
                out.push(cur);
                cur = e;
            }
        }
        out.push(cur);
        out
    }
}

/// A minimal two-parameter run-length model.
///
/// - Each error element (Sub, Ins, Del) pays: open_cost + (len-1) * ext_cost
/// - Match elements contribute zero (kept simple by design)
pub struct RunCostModel {
    pub open_cost: f64,
    pub ext_cost: f64,
}

impl RunCostModel {
    pub fn new(open_cost: f64, ext_cost: f64) -> Self {
        Self {
            open_cost,
            ext_cost,
        }
    }

    pub fn log_odds_ratio(&self, a: f64, b: f64) -> f64 {
        (b - a).exp()
    }

    pub fn score_cigar(&self, cigar: &[CigarElem]) -> f64 {
        let mut total = 0.0_f64;
        let mut longest_match_run: usize = 0;

        for elem in cigar.iter() {
            match elem.op {
                CigarOp::Match => {
                    if (elem.cnt as usize) > longest_match_run {
                        longest_match_run = elem.cnt as usize;
                    }
                }
                CigarOp::Sub | CigarOp::Ins | CigarOp::Del => {
                    let len = elem.cnt.max(1) as f64;
                    total += self.open_cost + (len - 1.0) * self.ext_cost;
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

/// Quality keys extracted from a CIGAR for parameter-free, lexicographic ranking.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CigarQualityKeys {
    pub longest_match: usize,
    pub error_runs: usize,
    pub indel_runs: usize,
    pub total_errors: usize,
    pub tie_score: f64,
    pub total_len: usize,
}

impl CigarQualityKeys {
    /// Lexicographic comparison:
    /// - maximize longest_match (desc)
    /// - minimize error_runs (asc)
    /// - minimize indel_runs (asc)
    /// - minimize total_errors (asc)
    /// - minimize tie_score (asc)
    pub fn cmp_lex(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        match other.longest_match.cmp(&self.longest_match) {
            Ordering::Equal => {}
            ord => return ord,
        }
        match self.error_runs.cmp(&other.error_runs) {
            Ordering::Equal => {}
            ord => return ord,
        }
        match self.indel_runs.cmp(&other.indel_runs) {
            Ordering::Equal => {}
            ord => return ord,
        }
        match self.total_errors.cmp(&other.total_errors) {
            Ordering::Equal => {}
            ord => return ord,
        }
        self.tie_score
            .partial_cmp(&other.tie_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

/// Compute quality keys for a CIGAR using the provided run model for tie_score.
pub fn cigar_quality_keys(cigar: &[CigarElem], run_model: &RunCostModel) -> CigarQualityKeys {
    let mut longest_match: usize = 0;
    let mut error_runs: usize = 0;
    let mut indel_runs: usize = 0;
    let mut total_errors: usize = 0;
    let mut total_len: usize = 0;
    for e in cigar.iter() {
        total_len += e.cnt as usize;
        match e.op {
            CigarOp::Match => {
                if (e.cnt as usize) > longest_match {
                    longest_match = e.cnt as usize;
                }
            }
            CigarOp::Sub => {
                error_runs += 1;
                total_errors += e.cnt as usize;
            }
            CigarOp::Ins | CigarOp::Del => {
                error_runs += 1;
                indel_runs += 1;
                total_errors += e.cnt as usize;
            }
        }
    }
    let tie_score = run_model.score_cigar(cigar);
    CigarQualityKeys {
        longest_match,
        error_runs,
        indel_runs,
        total_errors,
        tie_score,
        total_len,
    }
}

/// Single-parameter, simple quality score: (Lmax - ErrBases) / L
/// Accept if >= tau (e.g., tau=0.10). Higher favors long clean runs.
pub fn simple_quality_score(keys: &CigarQualityKeys) -> f64 {
    if keys.total_len == 0 {
        return 0.0;
    }
    (keys.longest_match as f64 - keys.total_errors as f64) / keys.total_len as f64
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
        } else {
            marker = "";
        }

        // Fixme: bound check?
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
    Some(summed_cost as i32)
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
    println!("region cost: {}", region_cost);

    let overhang_cost = ((left_overhang + right_overhang) as f32 * alpha).ceil() as i32;

    Some(region_cost + overhang_cost)
}

pub fn get_matching_region(m: &Match, start: usize, end: usize) -> Option<(usize, usize)> {
    let path = m.to_path();
    let (start, end) = (start as i32, end as i32);

    let mut sub_path = path
        .iter()
        .filter(|Pos(q_pos, _)| *q_pos >= start && *q_pos <= end);
    let start = sub_path.next()?.0 as usize;
    let end = sub_path.last()?.0 as usize;

    Some((start, end))
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

/// Convert a raw Lodhi score to a fraction of the maximum possible score for a perfect
/// contiguous match of length `l` using the same `k` and `lambda_decay`.
pub fn score_to_frac(score: f64, l: usize, k: usize, lambda_decay: f64) -> f64 {
    if k == 0 || l < k {
        return 0.0;
    }
    let mut lodhi = Lodhi::new(k, lambda_decay);
    let perfect = Cigar {
        ops: vec![CigarElem {
            op: CigarOp::Match,
            cnt: l as i32,
        }],
    };
    let smax = lodhi.compute(&perfect);
    if smax <= 0.0 {
        0.0
    } else {
        (score / smax).clamp(0.0, 1.0)
    }
}

/// Convert a fraction (0..1) back to a raw Lodhi score scale for a target length `l`,
/// given `k` and `lambda_decay`.
pub fn frac_to_score(frac: f64, l: usize, k: usize, lambda_decay: f64) -> f64 {
    if k == 0 || l < k {
        return 0.0;
    }
    let mut lodhi = Lodhi::new(k, lambda_decay);
    let perfect = Cigar {
        ops: vec![CigarElem {
            op: CigarOp::Match,
            cnt: l as i32,
        }],
    };
    let smax = lodhi.compute(&perfect);
    (frac.clamp(0.0, 1.0)) * smax
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
    use pa_types::CigarOp::*;
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

        let score = p.score_cigar(&best_match.cigar.ops, true);
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
        let score = ProbModel::new(0.01, 0.15).score_cigar(&cigar.ops, true);
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
    fn test_simple_kmer_score_counts() {
        let l5 = vec![CigarElem { op: Match, cnt: 5 }];
        let l6 = vec![CigarElem { op: Match, cnt: 6 }];
        let l7 = vec![CigarElem { op: Match, cnt: 7 }];

        let s5 = simple_kmer_score(&l5, 1.0, 1.0, 6);
        let s6 = simple_kmer_score(&l6, 1.0, 1.0, 6);
        let s7 = simple_kmer_score(&l7, 1.0, 1.0, 6);

        // Longer contiguous matches should score higher
        assert!(s6 > s5);
        assert!(s7 > s6);

        let with_sub = vec![
            CigarElem { op: Match, cnt: 6 },
            CigarElem { op: Sub, cnt: 1 },
        ]; // pattern_len = 7
        let pure_match = vec![CigarElem { op: Match, cnt: 7 }]; // pattern_len = 7
        let s_err = simple_kmer_score(&with_sub, 1.0, 1.0, 6);
        let s_pure = simple_kmer_score(&pure_match, 1.0, 1.0, 6);
        assert!(s_err < s_pure);
    }

    #[test]
    fn test_two_cigars() {
        let c1 = &[
            CigarElem { op: Sub, cnt: 2 },
            CigarElem { op: Match, cnt: 3 },
        ];
        let c2 = &[
            CigarElem { op: Match, cnt: 1 },
            CigarElem { op: Sub, cnt: 1 },
            CigarElem { op: Match, cnt: 1 },
            CigarElem { op: Sub, cnt: 1 },
            CigarElem { op: Match, cnt: 1 },
        ];
        let s1 = simple_kmer_score(c1, 1.0, 1.0, 3);
        let s2 = simple_kmer_score(c2, 1.0, 1.0, 3);
        println!("s1: {}, s2: {}", s1, s2);
    }

    #[test]
    fn get_frac_for_score() {
        let score = 1.0;
        let frac = score_to_frac(score, 44, 3, 0.5);
        println!("Fraction: {}", frac);
    }

    #[test]
    fn get_score_for_frac() {
        let frac = 0.5;
        let score = frac_to_score(frac, 44, 3, 0.5);
        println!("Score: {}", score);
    }
}
