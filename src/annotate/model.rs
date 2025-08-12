use needletail::{FastxReader, Sequence, parse_fastx_file};
use pa_types::Pos;
use sassy::profiles::*;
use sassy::*;

#[derive(Clone, Debug)]
struct ErrorModel {
    // probs[row_ref][col_read] over A,C,G,T,- (index 4 is '-')
    probs: [[f64; 5]; 5],
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum StepType {
    MatchSubst = 0,
    Insertion = 1,
    Deletion = 2,
}

impl StepType {
    fn index(self) -> usize {
        self as usize
    }
}

#[derive(Clone, Debug)]
pub struct Model {
    alpha: f64,

    // Raw counts
    confusion_counts: [[u64; 5]; 5], // rows ref A,C,G,T,- ; cols read A,C,G,T,-
    step_counts: [u64; 3],           // counts of M, I, D
    init_counts: [u64; 3],           // first step type per path
    transition_counts: [[u64; 3]; 3], // prev_step -> next_step

    // Learned probabilities
    // Emissions
    subst_probs: [[f64; 4]; 4], // P(read_base | ref_base) for M steps
    ins_base_probs: [f64; 4],   // P(inserted_base) for I steps
    del_ref_probs: [f64; 4],    // P(deletion | ref_base) for D steps

    // Step process
    init_probs: [f64; 3],            // P(step_0)
    transition_probs: [[f64; 3]; 3], // P(step_t | step_{t-1})
}

impl Model {
    pub fn new(alpha: f64) -> Self {
        Self {
            alpha,
            confusion_counts: [[0; 5]; 5],
            step_counts: [0; 3],
            init_counts: [0; 3],
            transition_counts: [[0; 3]; 3],
            subst_probs: [[0.0; 4]; 4],
            ins_base_probs: [0.0; 4],
            del_ref_probs: [0.0; 4],
            init_probs: [0.0; 3],
            transition_probs: [[0.0; 3]; 3],
        }
    }

    /// Estimate empirical base frequencies from a reference slice.
    pub fn base_freqs(ref_seq: &[u8]) -> [f64; 4] {
        let mut counts = [0_u64; 4];
        for &b in ref_seq.iter() {
            if let Some(i) = base_to_index(b) {
                counts[i] += 1;
            }
        }
        let total: u64 = counts.iter().sum();
        if total == 0 {
            return [0.25, 0.25, 0.25, 0.25];
        }
        [
            counts[0] as f64 / total as f64,
            counts[1] as f64 / total as f64,
            counts[2] as f64 / total as f64,
            counts[3] as f64 / total as f64,
        ]
    }

    /// Sample a random DNA sequence of length `len` from given base frequencies.
    pub fn sample_random_seq(len: usize, freqs: [f64; 4]) -> Vec<u8> {
        use rand::distributions::{Distribution, WeightedIndex};
        use rand::thread_rng;
        let dist = WeightedIndex::new(&freqs)
            .unwrap_or_else(|_| WeightedIndex::new(&[1.0, 1.0, 1.0, 1.0]).unwrap());
        let bases = [b'A', b'C', b'G', b'T'];
        let mut rng = thread_rng();
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            let i = dist.sample(&mut rng);
            out.push(bases[i]);
        }
        out
    }

    /// Estimate null distribution statistics (mean, std) of log-likelihood by aligning `query`
    /// to random references sampled from `freqs` with given length.
    pub fn estimate_null_ll_stats(
        &self,
        query: &[u8],
        ref_len: usize,
        n_samples: usize,
        freqs: [f64; 4],
    ) -> Option<(f64, f64)> {
        if n_samples == 0 {
            return None;
        }
        let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(0.5);
        let mut lls: Vec<f64> = Vec::with_capacity(n_samples);
        let k = ((query.len() as f64) * 0.3) as usize;
        for _ in 0..n_samples {
            let ref_seq = Self::sample_random_seq(ref_len, freqs);
            let matches = searcher.search(query, &ref_seq, k);
            if matches.is_empty() {
                continue;
            }
            let best = matches.iter().min_by_key(|m| m.cost)?;
            let ll = self.score(&best.to_path(), query, &ref_seq);
            lls.push(ll);
        }
        if lls.is_empty() {
            return None;
        }
        let mean = lls.iter().copied().sum::<f64>() / (lls.len() as f64);
        let var = lls.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / (lls.len() as f64);
        let std = var.sqrt().max(1e-9);
        Some((mean, std))
    }

    /// Compute log-likelihood ratio against a null mean (Quick Bayes factor approx).
    pub fn llr_against_null_mean(
        &self,
        path: &[Pos],
        query: &[u8],
        ref_seq: &[u8],
        null_mean: f64,
    ) -> f64 {
        let ll = self.score(path, query, ref_seq);
        ll - null_mean
    }

    /// Convert ΔLL to likelihood ratio (exp(ΔLL)).
    pub fn lr_from_delta_ll(delta_ll: f64) -> f64 {
        delta_ll.exp()
    }

    pub fn train(&mut self, fastq_file: &str, query: &[u8]) -> Option<(sassy::Match, Vec<u8>)> {
        let mut reader = parse_fastx_file(&fastq_file).ok()?;
        let mut example: Option<(sassy::Match, Vec<u8>)> = None;

        let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(0.5);
        while let Some(record) = reader.next() {
            let seqrec = match record {
                Ok(r) => r,
                Err(_) => continue,
            };
            let norm_seq = seqrec.normalize(false);
            let ref_seq = norm_seq.as_ref();

            let k = (query.len() as f64 * 0.25) as usize;
            let matches = searcher.search(query, &norm_seq, k);
            if matches.is_empty() {
                continue;
            }

            let best = matches.iter().min_by_key(|m| m.cost).unwrap().clone();
            if example.is_none() {
                example = Some((best.clone(), ref_seq.to_vec()));
            }

            let path = best.to_path();
            if path.len() < 2 {
                continue;
            }

            self.tally_from_path(&path, query, ref_seq);
        }

        self.compute_probs();
        example
    }

    fn tally_from_path(&mut self, path: &[Pos], query: &[u8], ref_seq: &[u8]) {
        let mut prev_step: Option<StepType> = None;
        let mut first_marked = false;

        for w in path.windows(2) {
            let Pos(q_prev, r_prev) = w[0];
            let Pos(q_curr, r_curr) = w[1];
            let dq = q_curr - q_prev;
            let dr = r_curr - r_prev;

            let step = match (dq, dr) {
                (1, 1) => StepType::MatchSubst,
                (1, 0) => StepType::Insertion,
                (0, 1) => StepType::Deletion,
                _ => continue,
            };

            if !first_marked {
                self.init_counts[step.index()] += 1;
                first_marked = true;
            }

            self.step_counts[step.index()] += 1;

            if let Some(prev) = prev_step {
                self.transition_counts[prev.index()][step.index()] += 1;
            }
            prev_step = Some(step);

            match step {
                StepType::MatchSubst => {
                    let qb = query[q_curr as usize];
                    let rb = ref_seq[r_curr as usize];
                    if let (Some(ri), Some(ci)) = (base_to_index(rb), base_to_index(qb)) {
                        self.confusion_counts[ri][ci] += 1;
                    }
                }
                StepType::Insertion => {
                    let qb = query[q_curr as usize];
                    if let Some(ci) = base_to_index(qb) {
                        self.confusion_counts[4][ci] += 1;
                    }
                }
                StepType::Deletion => {
                    let rb = ref_seq[r_curr as usize];
                    if let Some(ri) = base_to_index(rb) {
                        self.confusion_counts[ri][4] += 1;
                    }
                }
            }
        }
    }

    fn compute_probs(&mut self) {
        let a = self.alpha;
        // Substitution emissions
        for r in 0..4 {
            let sum: f64 = (0..4).map(|c| self.confusion_counts[r][c] as f64).sum();
            let denom = sum + a * 4.0;
            for c in 0..4 {
                self.subst_probs[r][c] = ((self.confusion_counts[r][c] as f64) + a) / denom;
            }
        }
        // Insertion emissions (inserted base distribution)
        let sum_ins: f64 = (0..4).map(|c| self.confusion_counts[4][c] as f64).sum();
        let denom_ins = sum_ins + a * 4.0;
        for c in 0..4 {
            self.ins_base_probs[c] = ((self.confusion_counts[4][c] as f64) + a) / denom_ins;
        }
        // Deletion emissions (by ref base)
        let sum_del: f64 = (0..4).map(|r| self.confusion_counts[r][4] as f64).sum();
        let denom_del = sum_del + a * 4.0;
        for r in 0..4 {
            self.del_ref_probs[r] = ((self.confusion_counts[r][4] as f64) + a) / denom_del;
        }

        // Initial step distribution
        let sum_init: f64 = self.init_counts.iter().map(|&x| x as f64).sum();
        let denom_init = sum_init + a * 3.0;
        for s in 0..3 {
            self.init_probs[s] = ((self.init_counts[s] as f64) + a) / denom_init;
        }

        // Transition probabilities
        for s in 0..3 {
            let row_sum: f64 = (0..3).map(|t| self.transition_counts[s][t] as f64).sum();
            let denom = row_sum + a * 3.0;
            for t in 0..3 {
                self.transition_probs[s][t] = ((self.transition_counts[s][t] as f64) + a) / denom;
            }
        }
    }

    pub fn score(&self, path: &[Pos], query: &[u8], ref_seq: &[u8]) -> f64 {
        let mut ll = 0.0_f64;
        let mut prev_step: Option<StepType> = None;
        let mut first_applied = false;

        for w in path.windows(2) {
            let Pos(q_prev, r_prev) = w[0];
            let Pos(q_curr, r_curr) = w[1];
            let dq = q_curr - q_prev;
            let dr = r_curr - r_prev;

            let step = match (dq, dr) {
                (1, 1) => StepType::MatchSubst,
                (1, 0) => StepType::Insertion,
                (0, 1) => StepType::Deletion,
                _ => continue,
            };

            if !first_applied {
                // initial step probability
                let p0 = self.init_probs[step.index()].max(1e-300);
                ll += p0.ln();
                first_applied = true;
            }

            // transition
            if let Some(prev) = prev_step {
                let pt = self.transition_probs[prev.index()][step.index()].max(1e-300);
                ll += pt.ln();
            }
            prev_step = Some(step);

            // emission
            let pe = match step {
                StepType::MatchSubst => {
                    // Path nodes are between characters; the step consumes the previous bases
                    if q_curr <= 0 || r_curr <= 0 {
                        continue;
                    }
                    let qi = (q_curr - 1) as usize;
                    let ri = (r_curr - 1) as usize;
                    if qi >= query.len() || ri >= ref_seq.len() {
                        continue;
                    }
                    let qb = query[qi];
                    let rb = ref_seq[ri];
                    match (base_to_index(rb), base_to_index(qb)) {
                        (Some(ri), Some(ci)) => self.subst_probs[ri][ci],
                        _ => continue,
                    }
                }
                StepType::Insertion => {
                    // Consumes a query base only
                    if q_curr <= 0 {
                        continue;
                    }
                    let qi = (q_curr - 1) as usize;
                    if qi >= query.len() {
                        continue;
                    }
                    let qb = query[qi];
                    match base_to_index(qb) {
                        Some(ci) => self.ins_base_probs[ci],
                        None => continue,
                    }
                }
                StepType::Deletion => {
                    // Consumes a reference base only
                    if r_curr <= 0 {
                        continue;
                    }
                    let ri = (r_curr - 1) as usize;
                    if ri >= ref_seq.len() {
                        continue;
                    }
                    let rb = ref_seq[ri];
                    match base_to_index(rb) {
                        Some(ri) => self.del_ref_probs[ri],
                        None => continue,
                    }
                }
            };
            ll += pe.max(1e-300).ln();
        }

        ll
    }

    pub fn display(&self) {
        // Confusion matrix
        println!("\nConfusion matrix (rows=ref, cols=read):");
        println!("       A        C        G        T        -");
        for r in 0..5 {
            let row_label = idx_to_base(r);
            print!("{}", row_label);
            for c in 0..5 {
                print!("\t{:>8}", self.confusion_counts[r][c]);
            }
            println!();
        }

        // Per-base summary
        println!("\nPer-reference-base summary:");
        println!("ref\tcorrect\tsubst\tdel");
        for r in 0..4 {
            let correct = self.confusion_counts[r][r];
            let deletions = self.confusion_counts[r][4];
            let substitutions: u64 = (0..4)
                .filter(|&c| c != r)
                .map(|c| self.confusion_counts[r][c])
                .sum();
            println!(
                "{}\t{}\t{}\t{}",
                idx_to_base(r),
                correct,
                substitutions,
                deletions
            );
        }

        // Insertions summary
        let insertions_by_base: [u64; 4] = [
            self.confusion_counts[4][0],
            self.confusion_counts[4][1],
            self.confusion_counts[4][2],
            self.confusion_counts[4][3],
        ];
        let total_insertions: u64 = insertions_by_base.iter().sum();
        println!("\nInsertions by base (read vs gap in ref):");
        println!("A\tC\tG\tT\t(total)");
        println!(
            "{}\t{}\t{}\t{}\t{}",
            insertions_by_base[0],
            insertions_by_base[1],
            insertions_by_base[2],
            insertions_by_base[3],
            total_insertions
        );

        // Step process
        println!(
            "\nInitial step distribution (M,I,D):\n{:.4}\t{:.4}\t{:.4}",
            self.init_probs[0], self.init_probs[1], self.init_probs[2]
        );
        println!("\nTransition matrix P(next|prev) over (M,I,D):");
        println!("      M       I       D");
        for s in 0..3 {
            let label = match s {
                0 => 'M',
                1 => 'I',
                2 => 'D',
                _ => '?',
            };
            print!("{}", label);
            for t in 0..3 {
                print!("\t{:.4}", self.transition_probs[s][t]);
            }
            println!();
        }
    }
}

fn analyze(fastq_file: &str, query: &[u8]) -> (sassy::Match, Vec<u8>, ErrorModel) {
    let mut reader = parse_fastx_file(&fastq_file).expect("valid path/file");

    // Confusion matrix with gap '-': rows = reference {A,C,G,T,-}, cols = read/query {A,C,G,T,-}
    // counts[row][col]
    let mut counts: [[u64; 5]; 5] = [[0; 5]; 5];

    // Keep one example best match and its backing ref sequence for scoring demo
    let mut example_match: Option<sassy::Match> = None;
    let mut example_ref_seq: Option<Vec<u8>> = None;

    let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(0.5);
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let norm_seq = seqrec.normalize(false);
        let ref_seq = norm_seq.as_ref();

        // Locate, semi-global with high confidence
        let matches = searcher.search(query, &norm_seq, 8);

        if matches.is_empty() {
            continue;
        }

        // Get best match (owned)
        let best_match = matches.iter().min_by_key(|m| m.cost).unwrap().clone();

        // Save example for scoring
        if example_match.is_none() {
            example_match = Some(best_match.clone());
            example_ref_seq = Some(ref_seq.to_vec());
        }

        let path = best_match.to_path();
        if path.len() < 2 {
            continue;
        }

        // Tally counts from path
        tally_counts_from_path(&path, query, ref_seq, &mut counts);
    }

    print_confusion_matrix(&counts);
    print_per_base_summary(&counts);

    // Build model and show example scoring API (no-op if file had no matches)
    let model = build_error_model(&counts, 0.5);
    println!(
        "\nModel built (Laplace alpha=0.5). You can now score any path with score_alignment_path_log_likelihood(..).\n"
    );

    let m = example_match.expect("no matches found in input to construct example");
    let ref_seq_vec = example_ref_seq.expect("no reference sequence captured for example match");
    (m, ref_seq_vec, model)
}

fn tally_counts_from_path(path: &[Pos], query: &[u8], ref_seq: &[u8], counts: &mut [[u64; 5]; 5]) {
    for w in path.windows(2) {
        let Pos(q_prev, r_prev) = w[0];
        let Pos(q_curr, r_curr) = w[1];
        let dq = q_curr - q_prev;
        let dr = r_curr - r_prev;

        match (dq, dr) {
            // Diagonal: consume both → match/substitution
            (1, 1) => {
                let qb = query[q_curr as usize];
                let rb = ref_seq[r_curr as usize];
                if let (Some(ri), Some(ci)) = (base_to_index(rb), base_to_index(qb)) {
                    counts[ri][ci] += 1;
                }
            }
            // Horizontal: consume query only → insertion in read relative to reference
            (1, 0) => {
                let qb = query[q_curr as usize];
                if let Some(ci) = base_to_index(qb) {
                    counts[4][ci] += 1;
                }
            }
            // Vertical: consume reference only → deletion in read (gap opposite this ref base)
            (0, 1) => {
                let rb = ref_seq[r_curr as usize];
                if let Some(ri) = base_to_index(rb) {
                    counts[ri][4] += 1;
                }
            }
            _ => {}
        }
    }
}

fn build_error_model(counts: &[[u64; 5]; 5], alpha: f64) -> ErrorModel {
    let mut probs = [[0.0_f64; 5]; 5];

    for r in 0..5 {
        // For gap row (r=4), only columns A,C,G,T are meaningful
        let cols: &[usize] = if r == 4 {
            &[0, 1, 2, 3]
        } else {
            &[0, 1, 2, 3, 4]
        };
        let sum_counts: f64 = cols.iter().map(|&c| counts[r][c] as f64).sum::<f64>();
        let k = cols.len() as f64;
        let denom = sum_counts + alpha * k;

        for c in 0..5 {
            if r == 4 && c == 4 {
                probs[r][c] = 0.0; // unused
                continue;
            }
            if !cols.contains(&c) {
                probs[r][c] = 0.0;
                continue;
            }
            probs[r][c] = ((counts[r][c] as f64) + alpha) / denom;
        }
    }

    ErrorModel { probs }
}

fn score_alignment_path_log_likelihood(
    path: &[Pos],
    query: &[u8],
    ref_seq: &[u8],
    model: &ErrorModel,
) -> f64 {
    let mut ll = 0.0_f64;

    for w in path.windows(2) {
        let Pos(q_prev, r_prev) = w[0];
        let Pos(q_curr, r_curr) = w[1];
        let dq = q_curr - q_prev;
        let dr = r_curr - r_prev;

        let p = match (dq, dr) {
            (1, 1) => {
                let qb = query[q_curr as usize];
                let rb = ref_seq[r_curr as usize];
                match (base_to_index(rb), base_to_index(qb)) {
                    (Some(ri), Some(ci)) => model.probs[ri][ci],
                    _ => continue, // skip ambiguous bases
                }
            }
            (1, 0) => {
                let qb = query[q_curr as usize];
                match base_to_index(qb) {
                    Some(ci) => model.probs[4][ci], // insertion
                    None => continue,
                }
            }
            (0, 1) => {
                let rb = ref_seq[r_curr as usize];
                match base_to_index(rb) {
                    Some(ri) => model.probs[ri][4], // deletion
                    None => continue,
                }
            }
            _ => continue,
        };

        // Avoid -inf on numerical issues
        let p_clamped = if p <= 0.0 { 1e-300 } else { p };
        ll += p_clamped.ln();
    }

    ll
}

pub fn score_match_log_likelihood(
    m: &sassy::Match,
    query: &[u8],
    ref_seq: &[u8],
    model: &ErrorModel,
) -> f64 {
    let path = m.to_path();
    score_alignment_path_log_likelihood(&path, query, ref_seq, model)
}

fn base_to_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None, // skip ambiguous bases
    }
}

fn idx_to_base(i: usize) -> char {
    match i {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        4 => '-',
        _ => '?',
    }
}

fn print_confusion_matrix(counts: &[[u64; 5]; 5]) {
    println!("\nConfusion matrix (rows=ref, cols=read):");
    println!("       A        C        G        T        -");
    for r in 0..5 {
        let row_label = idx_to_base(r);
        print!("{}", row_label);
        for c in 0..5 {
            print!("\t{:>8}", counts[r][c]);
        }
        println!();
    }
}

fn print_per_base_summary(counts: &[[u64; 5]; 5]) {
    println!("\nPer-reference-base summary:");
    println!("ref\tcorrect\tsubst\tdel");
    for r in 0..4 {
        // only A,C,G, T rows
        let correct = counts[r][r];
        let deletions = counts[r][4];
        let substitutions: u64 = (0..4).filter(|&c| c != r).map(|c| counts[r][c]).sum();
        println!(
            "{}\t{}\t{}\t{}",
            idx_to_base(r),
            correct,
            substitutions,
            deletions
        );
    }

    // Insertions summary (independent of reference base)
    let insertions_by_base: [u64; 4] = [counts[4][0], counts[4][1], counts[4][2], counts[4][3]];
    let total_insertions: u64 = insertions_by_base.iter().sum();
    println!("\nInsertions by base (read vs gap in ref):");
    println!("A\tC\tG\tT\t(total)");
    println!(
        "{}\t{}\t{}\t{}\t{}",
        insertions_by_base[0],
        insertions_by_base[1],
        insertions_by_base[2],
        insertions_by_base[3],
        total_insertions
    );
}

fn main() {
    // Train the comprehensive model
    let query = b"GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";
    let mut model = Model::new(0.5);
    let (best_match, ref_seq) = model
        .train(
            "/home/solprof/PhD/barbell-sassy-rewrite/pass_sample.fastq",
            query,
        )
        .expect("training produced no matches");

    // Display learned stats
    model.display();

    // Example scoring on provided sequences
    let read = b"TTATGCATCAACTGACCGTTATGCAAAAAGCCGATATTGCCACACCAGACACGCTGACGAACACCACTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTGGTGCGCCACTTCAGCCGCCAACCCACCTGGGCCGGCCACATAGCCACCAAAGCCGGCTACCGCATCCACCTTCAGCTGTTTCATGGCTTTCATGGCGCTGAGCGTGGCGGTTAAAATTTTAAAAGGCGCAGCCAGTTTACGCACAATCCCATTGCCGCGTACACCCTGAATATTGATCTGATAAATTGGAATATTATGGTTTTTTAAGAGCCGGTTTTCCATACCCGCTGGAGTCGCCAGCCAGGAGACTTCGATGCCTTGTTGTTGTAAATCTTTTGCTACAGCTAATGCGGGAAATACATGTCCACCGGTACCTGCGGCCATCATCATGACATGTTTAGGCTGTTTCTTCTGAGCATCGGTCACGGTCTAATTCTTCAATTCTAAAATAAAAATGGGGGAATAGGATTTAAGTCATCAATTGTAATCACAAACGCCAGGGCTTGAAAGTTTATTTACTTTATTTCAACCAAGCATTTCATTGTTTTTTTTCACGTGTTTTTGAAGTTTAAACTGTATTTATCATCAGGTGTTCATTTTTAATGATTTTTTTATTGAAAAAGCTTTTTTTATGATCATTTCTGTCTACCTCTAAATTACCTCTGCTCTCTACTTACACCTCATCACTCCTCCCCATTGCCGCACTTTTCACCGAAAAATAAGCCAGAAAATTGCTATTCTATGTCAAATGAAATCCCGCAGAGATTCTGTACCTAAAAGATTAGAGATGCTTAAGATTCAGTGTTTATGGATTCTTGAGCTTTCTTTTTCTGAAACTCATGAAGTGTAGGAGTACTGTATTTATCCATCCAGGGTGATTCTCCATAAAAAGAACCATTAAGTTGTGAATAATTAATCAACCTTTGCCGGAACTTTTTAAAATAAATACGATCCACTTTACCCGGATTACGCCAGAACCTGGCCAACGGGCGCTGGTTGGTTAGTTTAGGCAAGTCGTATAATGCGGAACCCAAAACTTTGACTGGAATTTGATGATATAAGGCTTGAATCCCCGTGGTACTGTTTACCGCCACAAAGCCCAAAGCATGCTTTAATAAGGTAGGTAAATGAATATCACAAAAATAGTGAATCCGCCCCTCCACACCGTGCTTGCGCGCAAGCTTGCGGATCAAACTGGCATAGTTACGATAACCACGATCCATCGGGTGGTGCTTTAAAACCAAGTGATGAGAACCTTCAGCATGCTGGGCAAAATCCTCTATCACCATCTGGGTATGCTTTTTCATGCATTTCAGTTCAGAATGCGTCCTGATCTGAAAGTCATTATGTACCTGCAAGGCAAACACAAAATAGCGTTTTTGATGGTGCCGAATAAAATGGTCAAATTTACGGCTTTCAAAGCAGGAATTCCTTAAGCGGCGAATTCCCGCCCAGATCCAGTAATACAGCTCTGTCCAGACACTCATTTGCCGATGGTGTTTATAATGCCGATATTTCCAGAAAAAATCAGCATAAAAAAGTAATAGGCCATCACGCTGAATATTCTAAGGCTGGTACTGCTATGAATTATCTCAACCTGCTCATTGACCTGACAGCTTTTTGCAGTCGTCGTTTGGAAATGCGTTAAAAAATGAGAAAAGAAATTTACCCCATTTTGCTCAAAGGTAATATAGTTCGGACGAATATAGCCTTCTTCAAAAGCAAAAAACTTCACTTTATTTTTTACCGCTATTTGTTTCGCAGTCTGATGATACTTACGGCAGTCGCCAAAACAGACCACTGCATCAATGTCATGCGTTAACAGTAGATGCTCGACCCATTGGGCAAAATCCGTATAGGTTCCTGTATAGTCAAAGCTATATTCAAAATGCCGCGCATCGAACTGATCGCCGCCATTAAAATTCACCTTAAAGCATTGAACATCATGAGATTTTAGCCAGTTGGCAAAATTAGAAAAAAAACCACCCATAGGACCCTGCAATAGCAGGACCCGGTCTGAAAGTAATAAATTATTTACATGTACAGACATGAGAAATTAACCCAGTTTTTTTAATAAATTATATAAAGGCGCCAGTATACTCTTTTTGTTCTTGTTTGCCTTTACTGCTTGGTCTAATTGTTGGTTAATAGATGTAATTACATCTTCGGGACGCACCACCGGTATTCCCATTTTTTTTGTGGTGGCTAAATTATAAATCGGGTACTCCACCAAAGTCATATAAATCAGCTGTTCCAGGCTAATTTTTTTGTTACGTCGTTCTGATTTGTACAGGTCCTGGGTCAGCCCCCACCCCGCATAAAAAGGCAAACCATAGCAATATACCTTTAAGCCGCGAATAAGGGCCTCAAAACCACTGAGCGATGTGATGGTATGTAACTCTCGTACAATTTCGAAACATTCAAGAATAGAAGCATTTAACTCTATCTGATCAGCATAACGCAATATCTCGGCAGCAGGAATACTACCGACCCTTAAACCAGTTTGTACATCAGGATGAGGCTTATAAATAATATAAGCTTCACGATTCAGTTCTCTAACTTTTCGAAGTAGCTCTAAATTGGTCTTAATGCCGACACCACCCAATTGAATCGACATATCATCTTCAACCTGCCCAACTACCAAAATAACTTTATGGCAATTTGGGCGGGTCAAAGCTTTTGCTTCTCCTACATTATATTTCGAAATATTGAGATCTATCAGTTTCTGTTGTAAATCACGTGCTCGCTGACTTTGTTCTTCAGTCAACTTAACTTGATTTAGAAGATTCTCAATTCTTGATGGCTTGGTCGCATCGTAATAAATCCCGACATCATCAAATACAAGTGAACAAGTGAGAATTAAAGTGGCACCGAGGCCAACTGAGCGAATAAAGCCATCCTCTACAGTGGTGACATTTTTAAAATCTTGATCCTTGAGATATTTTGCTTTCTTACCCCATGCTAGGATTAATTGATTTTTCTTAGGCTTAAAATAACGTTGGAATTTTAAACTCACCTTTGGAAAGCTTAAAAAATCTGCAATAAATTTCTTCTTCCAAGGACTAAATCCATAAGCGATATAAGAATGACTTAAGCGTTTTTGAATTTTAATATTGGGAATAATCAGATCGATAATATCTTCACACTCACAGCGCTTGCCTGTTATTGGATTGACATAACGTGCATAGTGAAAATAAGCCGACGCAAATAAATGACTAAGTGAACGGTTAATGCCACGCCGCCCCTGTAAAATATGTAATGGCGCATATTGATCATCTGTGAGCCCCCATCCTGCATACCAAGGAACCCCAAAACAGTGAACCTTTTTCTGGCATAGTAAGGCTTCAAAACCTAATTGTGAACTCACAACGTAAACTTCATCTATATGTTGAAGCAATTCAATCGGATTATAGTTATCCGTAAGCAATTGAATGTTTGAATGACTTAAATCTTCAGTTGAAAAATGGGATTTTGCCTTACCTGCCAAGACATCCGGATGTGTTTTCACCCAAATGATTGCTTCGGGATGATCAGCAAATGCTTGGCTTAACATTTTTTTAAATGTATCTGAGGAAGCACCGGCATACTTAATCGACTGATCACCAAAAGTTTGATCAATGACTAAAATATGTTTTCTCTGGCTAGAATTTTTCAGGATCTAAAGCTTTAAATTTTTGGTTATATTTAGTAATTCCATACTTTAAAATTTTTGAAATAGCCTCGTCAGCCCTTAAATTTTGTTCAGGAAGCTCAGGTTTTAATATTAAGTTTTCTAAATCGGAAGGTTGTAGAGCATCGAAATAAATTCCTGTTTTGTCTACAATAAGTGAAAAAGGCGGGTAGCCATCTTTGCCTAACCCCAAAGAACGAATAAATCCCTCCTCTAGGCAAAGTACATTTATACCTTTTTTTAGAGCTAGGGCTTTAGCTTTGAAAAAGCTTTTTTTGCGACCCCAACCTGCTGTATAGGTGTGATATTCCTTCAATGGATCCACCCAATGATCAAGGGAATGTTGATCTAAATAATTAACTAACATCATATATTAATAAAAATCCACCTGAATAGCATGCAGTTTAGATTTTATTTTCCTCGAAAATTTACATAACCATTGAATGCGTTGTCTAAAGGATAAATTGCATAAATGATCCTGTTGATACCCAATCAAGTAACTATGATGTTTTTTAACTTCTTCCAATTGTTTTTGAAAATCGATATTTGCTGCATCTACATTTGCAAGAGAACAATTTATATTATAAAAATGCTGTATACGCTTAGGTATAAGCAAAGATTGTAAAACATTATCTTCAAAAAGCACTAAGGTTGGAATATGTTTATTAGTTCGTTACATTCCATGAAATATCAACAAATGCTGCCGAAATCATTATATTTTTTTGCTGCATAGGATGTCAGAAGTGTGCCAACCTTTATAATTAAAATGATCTGCAAAGCAGATAGACTGCTTTACAATCACTACTTCAATAGACAATTTATTAAATGCCTTAATCATCGTTTTAAAACGTTTTAAATAGCTTCGATCACGATGACCTAAATTTGGATGCGGTTGAATATAAAGCTGATGAATACCTTTCTTTTGAATTATATCTATCAAAAGTTTCGTAGCTTTCTCACCAATAGATGAATAGTTTTCAGAAGGTACGCCGCCTTCCCATGTTGGAGCATACAGTAAACTGTAGCTATCTGTGGAGTAGACATAAGGATTTTGCCCAATAAATGTATCTCCCATTTCATTTACTCGGAATTCCTGCTCAATATTATAAAGATTAAAAATTCCAGCTTTTAGAAAACGTTCAATTCCAACTTGACCAGAAGTAATCACAAAGTCATAAATGCGAAAAATTGGTTTAATAGATGCGAGTTTATGACTTTCACCATGACTAACAAAAATATGTGTGAGATCACGTCGTGCTGTAATTTTACAGTTCGATTGAGCATTAAATAAATAAAAGACTGTTTTGTGACTTTCAAATATCAAATCACTATTTTTTATGAATAAATTTAAAACGAATACCTGCGGTAACAGTTTTTCTTTATATAATTTTTAATATCTTTATAAAATTTAAAATATACACACAATACGATGCGTTGCATTAAACACATCGGTCTGAAAATAACTCATCAGTTGATCTATAGTCCCTTTCAATTCTTGAGTCTATATAGATGTAATACAAGTTATTTTCAGACATATAAGCTTCTTATCGACCGTCTTTTCGAACGTCGACAATTATACCCGTATTTTTAGTACCAAGCACCTTGGATTAAGGCTACTTCCTGTGCAGATAGAAGTAAATTGTATCACTGTATGCACCAAAGTTTGCTGTACGCATTGGTGTTGCAGTACGTTCAGGATTAATACAAGTTACGCAGATCCTCTAATAAGCCCACTCTTCAGCAATAGCTTGTGTAAAGTTTACAATTGCAGCATTGGTAGAAGAATACAGAGCATAGAACGCAGGCCACGCGTATATGAACTCGATGTAAAGTTCATTAGCATACCAAATGTTTCTAGCAAATATGGTTTAGCTGCGATAGCTACATTTACTGCGCCAGTATAATTCACACCAATTAAAGATGTTATTTCTTCAAGTGTCAGTGTATCAATTGGTTTTTTAATCAGTAGGCCTGCGGTATTCACGATATAGTCAATACGACCTAATTTTTGCTTCATTTCTGCAAGATAAGTTTTTACACTATCTAAATTACAGATGTCTACCCGATTATAAGAACGACTAGCAACCTCAACATGCGCACCGTGCAGTGTCGCAATGTTTTTGATTTCTAAACCAATACCACTAGTACCACCAAAAATGACGATATTTTTGCCATTTAAATAATCAAAACTATCACCCTCACGGAATGGAATTCCATGTGCAGTTGTTCAACTTTTTCAGCGATGAACAGATCTACAGGGTGCGTCACTTTCAAATTACGTTCACTACCTAGGACTGTTGCCACACGTGTTCTAGGCAACATGGCACCAATAACACCGCAGTCACATGTGAATACTTTGCGATTTTGTTCTAACGCTTTTTTATAAGCTAATTGAATCACACCCAATTTGAAAGCCTGAGGTGTTTGACCACGTCGCATATGCGAGCGTTGCGGAATATTTGCAATACATCCATCATCATATACTTCAACTAAAGTATCTGCTGAAGGAGTCACTACATCTACTGCATCAAAGTCTTTCAATGCTTCAATGCAATTACTAATAATTTCCTGATTCACCAATGGACGTACAGCATCATGAAATAGAATATTGCATGAATTATCATATTCTTTTAATGCCGTTAATGCCGAATAGGTAGAATCAAAACGCTCATCACCACCCGTCACAACTTTAGAAATTTTTGACCATTTATTTTTTTTGATCAAATCCCACGTACGATTAATATGATGTGTTTGTGATACGATGAAGATTTCATCAATATCATCTGAGTTTTGAAAAACATCAATGGTATACTCAATAATTGCTTTACCAGCGAGTTTTGTAAATTGTTTTGGCAGTGCACCGCCAAAGCGCGACCCTGAACCACCCGATAAAATTACTGCTATGTTCTTCATGTAATCCCTTGATAAGAAACCAAGATTCTCAATTAATATAACTATAAACTCACAAATTACTGTCTTACAGGCAATATATAGTTTGGAGTCAAACCGTACATACAAGATACATCATATAGTTTCTTTTCAAAAGTGTAATCACATGTTTCATTTGGTGAAATAATATACTCAATCGCTTTTTTACGATTTTCTGCCAAAAAATCATTGCCATTCAGGACCACCTCCTGAATTAACTCAGATAATTCAGCTATACTAGAAAAAGGGTAACAATAATCAGAATATTGCATATCACTTTGAGCAATCACAATTTCTTTATTTTTAGGAATATAGATAAAAATAGGGGCATCACTTGCGAGACACTCTGAAACAACTGCTGAAATATCACAAATATAAATATTTGAAGATTTAATTAATGTTTCAACCGAGACTTCCCTCGGATGATTAATCAAATTAATCTCTAATTCAGACAATGATAATCGAACCTTTTCATCATATCCCTTTAAATTAATATCTCGACTTCCTGTAACAGGATGATATTTGGACTGTAAGAGATAGTCTGGTAATAATGCTGAGAGCAACGTGAAAATGTCTTGGCTTAAAGAAACTGAACTATAGTTACTTTCTTCATATACCCCCTCCCAAGTAGGTAGGTATAACACACTTTTAGACTCACGTCGGTTCCATACAACTTCAGATGTTTTTAAAACTTCCCGAAGAGAAGGCCGGCCCACTTTAGCAAAAGACATCAGATCAACATCAAATCCTGCATTTTTAAAACGATCTATATGCGCTTGGCCAGCAACCCAAATTTCATCATATACTCTAAAGTACTTATGTGCACTGGCTGATTTATCACTATCACCATGACCTAAAAATATATGATTTACTTCATTAAATTTCAAAGTGTGTAATGTATTACCAGTATTTGAACTGTAATAAATATTTCTCACATAATAAAGTTTTGAAAATAATGAATCTATATCGCCGACGCCTTTTGCATAAGCAACATTAAACCTAGGATAATTTTTTCTGACCCACTGGTATAAATCCAAATTCCTTACTAAAATCAATGTATTTTCAGAAACAACACTAAAATATTTAAGCCAAAAATTTAAATGTGAAGGCCCATGAACCATCCCTTCTCCAGTATGGAATACATGGTGATACGGGTTAAAAAAATTTACATTTTCTATATTAATATTCACTTCAACACCATTTTAAAATTTTATTTAATTAAAAGAGATTTTAAATTTTCAGCAAAAACTTCTTTTGATTTCTCATAGAATCCAAATCTATATTTCAAAAATTCTTTATACTCAATTAAAAACTCTTCATTAGACTCTAATTGATGAATAGCATGATAAAACTCATCTACACTAGACAGTAAAAAAATATTTTCATTAAACTTAGAATTACTAATGTTAAATAAAGCATCTAAATTATTTGGTTTAAAATAAAATATTGGAATATTAATTGCTGCTATTTTCTGAAAGGAAAAATCATCTTCAAATCAGAAAATACAACATCAAAAATATTGATATTATTGGAAAAATTACCAAAGTTAACATCAATATTATTTACTTTCACAATATTTTTTTAAATCAGTCATTTGTCTAGATTTTATATTTTCAACATAAAAAATACATGATTCGATAATTTTTTATGAAGATTTAAAAATATAGAAAGAAGGCTATCTACAATTTGTAAGTTTGGATTAATCTTAATATAGATTGTTCGCACAATATTCGATTTTGAAATAAATGAATTACTGCGAGGATTACCAATATTTAGAACTTTCAAATGTCGCGTACTTAATTCATTTAAAAAGCAATCTTTTGTATATTCATTGTTTGTCCAAATTTCATCATAAGCTCTAAAAAATTTATGAACACCTGGCTGTCTATGCAGGTCAAAATTCCCACAAAAAATATGTTGATATTGATTAAATCGAATCAAATGAATATTATTTGAAGTATTTGAAACAAAGCAAATCTTTTTAAGATTGGGTAGGCGAAGCAGTACCTTTTCGATATCTGTAGGAGATTTTGCATATACAATATTAATATGACTATATTCTTTTCTTAAATTATTATAAAGATCTGAATTCCTAACCAAGATTACAGTATTTTTTCTGAAATCATAGGTATATTTAATCCAATCCCGACATATAGTCATTCCTACAGCAGCTTCTTCGCCACAATGAAAAATCGTGGTCACCTGATATTCATGGATAAGAGATGTTACATTAAAAGTAATTTGCTGGACAATATTAGGCTTCTTGGCTGCAGAAGGCTTACTTTTCTGCTCATTTGGCTTTGCCTTAACTTCAGCTTTTTTAGCCACCGCTAGTTTGGTCTGCTTATTAATAGTTACAGCAGTAGGTTTTGAGTTTAAAAAATTATTTTCTAAATACGTTACTCTCTTCTTAACAAAGTCTTTAAAAAACTTATCAGGAGAATTGATTAATTTCGTGATTTTACTCATTTTGAGGCCTTTAAAATAAAAATTATACAAATAGTGCCTTACCCTCAGATATATCATTAATTACTTCTTGATATTCTGAGTCATACATACTTAAATCAATTTTATTTTTCAACCCCATTCCTTGCAAAGCAAGTTTATATTTATATCTTAATTTAAACATTCCTAAAGTTTCTTCAAAAAGATTTAAAGAATCCTGATAGTTTTTAAAACCCAAATAAAGTTGATTCCTACGCTCTAATCTATTTTTTCTTAAAGGATCATCCACTATTCCGACCAAGAATTCATCCAATGATATTTTTTGATAATCCCATAAATAACAAAAATGATCTAATCTAGCTTTTAATATATGTTGATAAATCCATGGCTTATCTATTGAACTTACAGGTAATATAATTGGTTTATTAGTAAACAAATACTCTGACAATACCCCAGATACATCACAAATAATTGCATCAGACCAATTAAATTGATCAGATTTATCACCTCCTAATTTCATAGTTCTATTGC";

    let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(0.5);
    let matches = searcher.search(query, read, 30);
    let best = matches.iter().min_by_key(|m| m.cost).unwrap();
    let path = best.to_path();
    let ll = model.score(&path, query, ref_seq.as_slice());
    println!("log-likelihood: {}", ll);

    let read33 = b"ATGTCATCTACTCGTTCAGTTATGCAAAAACCGGGTATTACGGCCGAAATTGATTAGAGAAGTACACGGTTTTCGCATTTATCGTGAACAGCTTCGCGTTTTTCGTGCGCCGCTTCAGTTGCAGCAAGATTCAATGAATTAGAACTTAAAGACCATCAAAAAGCGGATTATAAAATCAAGGCAAAACAATTTGTGAAAATTTACGGACAATTAGCATCTATATTACCATATGAAATATTGGATTGGGAAAAAATGTTTTGGTTTTTGAAATTTTTAATTCCCAAACTTAAAATAAAGGATGAAAGTCAAGATACATTGGATAAACTAATGGAATCGATAGATTTATCCACTTATGCCATACAGAGAGTAAGGCTGAATGAAAAAATATCGCTTGATGAGGAAGAAAAAGAGGTAAATCCGCAAAACCCAAATCCGAGAGGTGCCCACGATGAAACGGAAAAAGATGAGTTGGATACCATCTTATCGTCTTTTAACGAACGATGGTTTCAAGGCTGGGAAAGTACCCCAGAAGAGCAGAGAATGATTTTGATTTCGTTAGGGAAGTCTGTGCAAAACCATCCAGATTTTATGACTAAGTATATG";
    let query33 = b"GCTTGGGTGTTTAACCGGAGTTCGTCCAGAGAAGTACACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";

    let matches2 = searcher.search(query33, read33, 30);
    let best2 = matches2.iter().min_by_key(|m| m.cost).unwrap();
    let path2 = best2.to_path();
    let ll2 = model.score(&path2, query33, ref_seq.as_slice());
    println!("log-likelihood: {}", ll2);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gap() {
        let query = b"AAAAAATTTTTT";
        let text = b"AAAAAAGGTTTTTT";
        let mut searcher = Searcher::<Iupac>::new_fwd_with_overhang(0.5);
        let matches = searcher.search(query, text, 2);

        let m = matches.iter().min_by_key(|m| m.cost).unwrap();
        let _path = m.to_path();
        println!("{:?}", m.to_path());
    }
}
