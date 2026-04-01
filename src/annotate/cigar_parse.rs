use pa_types::CigarOp;
use pa_types::Pos;
use sassy::*;

/// Given a match and the
pub fn map_pat_to_text_with_cost(
    m: &Match,
    p_start: i32,
    p_end: i32,
) -> Option<((usize, usize), (usize, usize), u32)> {
    // Track positions (Pos) as offsets in the pattern and text
    let mut start_pair: Option<Pos> = None;
    let mut end_pair: Option<Pos> = None;

    // Track index in path, e.g. we can have gaps in text/pattern which
    // which make true positions (above) and index unaligned
    let mut start_idx: Option<usize> = None;
    let mut end_idx: Option<usize> = None;

    let path = m.to_path();

    for (idx, &Pos(i, j)) in path.iter().enumerate() {
        if i >= p_start && i < p_end {
            if start_pair.is_none() {
                start_pair = Some(Pos(i, j));
                start_idx = Some(idx);
            }
            end_pair = Some(Pos(i, j));
            end_idx = Some(idx);
        }
    }

    match (start_pair, end_pair, start_idx, end_idx) {
        (Some(Pos(pi, pj)), Some(Pos(ei, ej)), Some(si), Some(ei_idx)) => {
            // Compute cost for the subpath
            let cost = compute_subpath_cost(m, si, ei_idx);
            Some((
                (pi as usize, (ei + 1) as usize),
                (pj as usize, (ej + 1) as usize),
                cost,
            ))
        }
        _ => None,
    }
}

fn compute_subpath_cost(m: &Match, start_idx: usize, end_idx: usize) -> u32 {
    let mut cost = 0u32;
    let mut current_idx = 0;

    for el in &m.cigar.ops {
        for _ in 0..el.cnt {
            if current_idx >= start_idx && current_idx <= end_idx {
                cost += match el.op {
                    CigarOp::Match => 0,
                    _ => 1, // ins, del, sub all are 1 edit
                };
            }
            current_idx += 1;

            if current_idx > end_idx {
                return cost;
            }
        }
    }

    cost
}

/// Get matching region for a given match and start/end positions
pub fn get_matching_region(m: &Match, start: usize, end: usize) -> Option<(usize, usize)> {
    let path = m.to_path();
    let (start, end) = (start as i32, end as i32);

    let mut sub_path = path
        .iter()
        .filter(|Pos(q_pos, _)| *q_pos >= start && *q_pos <= end);

    let start = sub_path.next()?.1 as usize;
    let end = sub_path.next_back()?.1 as usize;
    Some((start.min(end), end.max(start)))
}

#[cfg(test)]
mod test {
    use super::*;
    use sassy::Searcher;
    use sassy::profiles::Iupac;

    fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        let mut rev = Vec::new();
        for el in seq.iter().rev() {
            rev.push(match el {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'N',
            });
        }
        rev
    }

    #[test]
    fn test_cost_extraction_no_edits() {
        let p = b"AAAAACCCAAAA";
        let t = b"GGGGAAAAACCCAAAAGGGGG";
        let mut searcher = Searcher::<Iupac>::new_rc();
        let matches = searcher.search(p, &t, 0);
        let m = matches.first().unwrap();
        let ((_ps, _pe), (_t_start, _t_end), cost) =
            map_pat_to_text_with_cost(m, 5, 7 + 1).unwrap();
        assert_eq!(cost, 0);
        // Rc'ing it should not change the cost
        let p_rc = reverse_complement(p);
        let t_rc = reverse_complement(t);
        let mut searcher = Searcher::<Iupac>::new_rc();
        let matches = searcher.search(&p_rc, &t_rc, 0);
        let m = matches.first().unwrap();
        let ((_ps, _pe), (_t_start, _t_end), cost) =
            map_pat_to_text_with_cost(m, 5, 7 + 1).unwrap();
        assert_eq!(cost, 0);
    }

    #[test]
    fn test_cost_extraction_1_edits() {
        let p = b"AAAAACCCAAAA";
        let t = b"GGGGAAAAACGCAAAA";
        let mut searcher = Searcher::<Iupac>::new_rc();
        let matches = searcher.search(p, &t, 1);
        let m = matches.first().unwrap();
        let ((_ps, _pe), (_t_start, _t_end), cost) =
            map_pat_to_text_with_cost(m, 5, 7 + 1).unwrap();
        assert_eq!(cost, 1);
    }

    #[test]
    fn test_cost_extraction_1_edits_overhang_left_flank() {
        let p = b"AAAAACCCAAAA";
        let t = b"ACGCAAAAGGGGGGGGGGGG";
        let mut searcher = Searcher::<Iupac>::new_rc();
        let matches = searcher.search(p, &t, 5);
        let m = matches.first().unwrap();
        let ((_ps, _pe), (t_start, t_end), cost) = map_pat_to_text_with_cost(m, 5, 7 + 1).unwrap();
        assert_eq!(cost, 1);
        assert_eq!(t_start, 1);
        assert_eq!(t_end, 4);
    }

    #[test]
    fn test_cost_extraction_1_edits_overhang_right_flank() {
        let p = b"AAAAACCCAAAA";
        let t = b"GAAAAACGC";
        let mut searcher = Searcher::<Iupac>::new_rc();
        let matches = searcher.search(p, &t, 5);
        let m = matches.first().unwrap();
        let ((_ps, _pe), (t_start, t_end), cost) = map_pat_to_text_with_cost(m, 5, 7 + 1).unwrap();
        assert_eq!(cost, 1);
        assert_eq!(t_start, 6);
        assert_eq!(t_end, 9);
    }

    #[test]
    fn test_cost_overhang_including_bar() {
        // The tests before included overhang of the flank, but this will include overhang of
        // the barcode itself as well
        let p = b"AAAAACCCAAAA";
        let t = b"GCAAAAGGGGGGGGGGGG"; // Note first "C" is missing of the barcode
        let mut searcher = Searcher::<Iupac>::new_rc();
        let matches = searcher.search(p, &t, 8);
        let m = matches.first().unwrap();
        let ((_ps, _pe), (t_start, t_end), cost) = map_pat_to_text_with_cost(m, 5, 7 + 1).unwrap();
        assert_eq!(cost, 2);
        assert_eq!(t_start, 0);
        assert_eq!(t_end, 2);
    }
}
