use pa_types::CostModel;
use pa_types::Pos;
use sassy::Searcher;
use sassy::*;

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
    }

    if end > sassy_match.pattern_end {
        let right_overhang = (end + 1) - sassy_match.pattern_end;
        assert!(end < p.len());
        let r_overhang_cost = (right_overhang as f32 * alpha.unwrap_or(0.0)).ceil() as i32;
        //println!("Adding right overhang: {}", r_overhang_cost);
        summed_cost += r_overhang_cost;
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
}
