use pa_types::Pos;
use sassy::*;

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

    for Pos(i, j) in m.to_path().iter() {
        if *i >= p_start && *i < p_end {
            if start_pair.is_none() {
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

/// Get matching region for a given match and start/end positions
pub fn get_matching_region(m: &Match, start: usize, end: usize) -> Option<(usize, usize)> {
    let path = m.to_path();
    let (start, end) = (start as i32, end as i32);

    let mut sub_path = path
        .iter()
        .filter(|Pos(q_pos, _)| *q_pos >= start && *q_pos <= end);
    let start = sub_path.next()?.0 as usize;
    let end = sub_path.next_back()?.0 as usize;

    Some((start, end))
}

#[cfg(test)]
mod test {
    use super::*;
    use sassy::Searcher;
    use sassy::profiles::Iupac;

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
}
