use crate::annotate::barcodes::BarcodeType;
use crate::annotate::searcher::BarbellMatch;

pub fn collapse_overlapping_matches(
    matches: &[BarbellMatch],
    filter_overlap: f32,
) -> Vec<BarbellMatch> {
    if matches.is_empty() {
        return Vec::new();
    }
    let mut sorted = matches.to_vec();
    sorted.sort_by_key(|m| m.read_start_flank);

    let mut groups = Vec::new();
    let mut group = vec![sorted[0].clone()];

    for m in sorted.into_iter().skip(1) {
        if group.iter().any(|g| is_overlap(g, &m, filter_overlap)) {
            group.push(m);
        } else {
            groups.push(std::mem::take(&mut group));
            group.push(m);
        }
    }
    groups.push(group);

    groups.into_iter().map(select_best_match).collect()
}

fn is_overlap(a: &BarbellMatch, b: &BarbellMatch, threshold: f32) -> bool {
    let start = a.read_start_flank.max(b.read_start_flank);
    let end = a.read_end_flank.min(b.read_end_flank);

    if end <= start {
        return false;
    }

    let overlap = end - start;
    let min_len =
        (a.read_end_flank - a.read_start_flank).min(b.read_end_flank - b.read_start_flank);
    (overlap as f32 / min_len as f32) >= threshold
}

fn select_best_match(group: Vec<BarbellMatch>) -> BarbellMatch {
    let mut candidates: Vec<_> = group.into_iter().collect();

    // Priority order: 1) Ftag/Rtag (barcode matches), 2) Fflank/Rflank (flank matches)
    candidates.sort_by(|a, b| {
        // Define priority levels
        let priority_a = match a.match_type {
            BarcodeType::Ftag | BarcodeType::Rtag => 1, // Highest priority - actual barcode matches
            BarcodeType::Fflank | BarcodeType::Rflank => 2, // Medium priority - full flank matches
        };

        let priority_b = match b.match_type {
            BarcodeType::Ftag | BarcodeType::Rtag => 1,
            BarcodeType::Fflank | BarcodeType::Rflank => 2,
        };

        // First sort by priority level
        priority_a.cmp(&priority_b).then_with(|| {
            // If both are Ftag/Rtag, prioritize by barcode_cost, then flank_cost
            if priority_a == 1 && priority_b == 1 {
                if a.barcode_cost == b.barcode_cost && a.flank_cost == b.flank_cost {
                    println!(
                        "tie found, cost: {} at pos: {}",
                        a.barcode_cost, a.read_start_flank
                    );
                }
                a.barcode_cost
                    .cmp(&b.barcode_cost)
                    .then_with(|| a.flank_cost.cmp(&b.flank_cost))
            } else if priority_a == 2 && priority_b == 2 {
                // If both are Fflank/Rflank, prioritize by longest flank length
                let length_a = a.read_end_flank - a.read_start_flank;
                let length_b = b.read_end_flank - b.read_start_flank;
                length_b.cmp(&length_a) // Longer length first (reverse order)
            //a.flank_cost.cmp(&b.flank_cost)
            } else {
                std::cmp::Ordering::Equal
            }
        })
    });

    candidates[0].clone()
}

#[cfg(test)]
mod tests {

    use rand::seq::SliceRandom;

    use super::*;
    use crate::annotate::barcodes::BarcodeType;
    use crate::filter::pattern::Cut;
    use pa_types::Cost;
    use sassy::Strand;

    fn create_match(
        read_start_bar: usize,
        read_end_bar: usize,
        read_start_flank: usize,
        read_end_flank: usize,
        bar_start: usize,
        bar_end: usize,
        match_type: BarcodeType,
        flank_cost: Cost,
        barcode_cost: Cost,
        label: &str,
        strand: Strand,
        read_len: usize,
        read_id: String,
        rel_dist_to_end: isize,
        cuts: Option<Vec<(Cut, usize)>>,
    ) -> BarbellMatch {
        BarbellMatch::new(
            read_start_bar,
            read_end_bar,
            read_start_flank,
            read_end_flank,
            bar_start,
            bar_end,
            match_type,
            flank_cost,
            barcode_cost,
            label.to_string(),
            strand,
            read_len,
            read_id,
            rel_dist_to_end,
            cuts,
        )
    }

    fn create_match_template(
        start: usize,
        end: usize,
        match_type: BarcodeType,
        barcode_cost: usize,
        label: &str,
    ) -> BarbellMatch {
        create_match(
            start,
            end,
            start,
            end,
            0,
            10,
            match_type,
            0,
            barcode_cost as Cost,
            &label,
            Strand::Fwd,
            100,
            "test".to_string(),
            0,
            None,
        )
    }

    #[test]
    fn test_empty_input() {
        // Test empty input
        let matches: Vec<BarbellMatch> = Vec::new();
        let result = collapse_overlapping_matches(&matches, 0.5);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_single_match() {
        // Test single match
        let matches = vec![create_match_template(0, 10, BarcodeType::Ftag, 3, "test1")];
        let result = collapse_overlapping_matches(&matches, 0.5);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].label, "test1");
    }

    #[test]
    fn test_doulbe_no_overlap() {
        // Test non-overlapping matches
        let matches = vec![
            create_match_template(0, 10, BarcodeType::Ftag, 3, "test1"),
            create_match_template(10, 20, BarcodeType::Ftag, 3, "test2"),
        ];
        let result = collapse_overlapping_matches(&matches, 0.5);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].label, "test1");
        assert_eq!(result[1].label, "test2");
    }

    #[test]
    fn test_collapse_overlapping_matches() {
        // Test overlapping matches - should collapse to best one
        let matches = vec![
            create_match_template(0, 20, BarcodeType::Ftag, 0, "test1"),
            create_match_template(15, 20, BarcodeType::Ftag, 3, "test2"),
        ];
        let result = collapse_overlapping_matches(&matches, 0.5);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].label, "test1");
    }

    #[test]
    fn test_overlap_threshold() {
        let matches = vec![
            create_match_template(0, 20, BarcodeType::Ftag, 0, "test1"),
            create_match_template(10, 35, BarcodeType::Ftag, 3, "test2"),
        ];
        // They now overlap at 5, which is 50% of the smaller length
        let result = collapse_overlapping_matches(&matches, 0.5);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].label, "test1");
        // If we now change it to 0.6, it should not collapse
        let result = collapse_overlapping_matches(&matches, 0.6);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].label, "test1");
        assert_eq!(result[1].label, "test2");
    }

    #[test]
    fn test_correct_sorting() {
        let mut matches = vec![
            create_match_template(0, 10, BarcodeType::Ftag, 0, "test1"),
            create_match_template(10, 20, BarcodeType::Ftag, 3, "test2"),
            create_match_template(0, 15, BarcodeType::Ftag, 3, "test2"),
            create_match_template(100, 110, BarcodeType::Ftag, 3, "test3"),
        ];
        // Random shuffle the matches vector
        let mut rng = rand::thread_rng();
        for _ in 0..10 {
            matches.shuffle(&mut rng);
            let result = collapse_overlapping_matches(&matches, 0.5);
            assert_eq!(result.len(), 2);
            assert_eq!(result[0].label, "test1");
            assert_eq!(result[1].label, "test3");
        }
    }

    #[test]
    fn test_small_ovlerap() {
        let mut matches = vec![
            create_match_template(0, 10, BarcodeType::Ftag, 3, "test1"),
            create_match_template(10, 20, BarcodeType::Ftag, 1, "test2"),
        ];
        for _ in 0..4 {
            matches[1].read_start_flank -= 1;
            matches[1].read_end_flank -= 1;
            let result = collapse_overlapping_matches(&matches, 0.5);
            println!(
                "Match 1 from {} to {}",
                matches[1].read_start_flank, matches[1].read_end_flank
            );
            assert_eq!(result.len(), 2);
            assert_eq!(result[0].label, "test1");
            assert_eq!(result[1].label, "test2");
        }
        // Going one more would put it at 15, and thus trigger collapse
        matches[1].read_start_flank -= 1;
        matches[1].read_end_flank -= 1;
        let result = collapse_overlapping_matches(&matches, 0.5);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].label, "test2");
    }
}
