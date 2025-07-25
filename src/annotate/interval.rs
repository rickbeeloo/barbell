use crate::annotate::searcher::BarbellMatch;
use sassy::Strand;

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
    candidates.sort_by(|a, b| {
        a.match_type
            .cmp(&b.match_type)
            .then_with(|| (a.barcode_cost).cmp(&(b.barcode_cost)))
    });
    candidates[0].clone()
}

// #[cfg(test)]
// mod tests {

//     use super::*;
//     use crate::search::barcodes::BarcodeType;

//     fn create_match(
//         read_start_bar: usize,
//         read_end_bar: usize,
//         read_start_flank: usize,
//         read_end_flank: usize,
//         bar_start: usize,
//         bar_end: usize,
//         match_type: BarcodeType,
//         flank_cost: Cost,
//         barcode_cost: Cost,
//         label: &str,
//         strand: Strand,
//         read_len: usize,
//     ) -> BarbellMatch {
//         BarbellMatch::new(
//             read_start_bar,
//             read_end_bar,
//             read_start_flank,
//             read_end_flank,
//             bar_start,
//             bar_end,
//             match_type,
//             flank_cost,
//             barcode_cost,
//             label.to_string(),
//             strand,
//             read_len,
//         )
//     }

//     #[test]
//     fn test_collapse_overlapping_matches() {
//         // Test empty input
//         let matches: Vec<BarbellMatch> = Vec::new();
//         let result = collapse_overlapping_matches(&matches, 0.5);
//         assert_eq!(result.len(), 0);

//         // Test single match
//         let matches = vec![create_match(
//             0,
//             10,
//             10,
//             20,
//             0,
//             10,
//             BarcodeType::Fbar,
//             5,
//             3,
//             "test1",
//             Strand::Fwd,
//         )];
//         let result = collapse_overlapping_matches(&matches, 0.5);
//         assert_eq!(result.len(), 1);
//         assert_eq!(result[0].label, "test1");

//         // Test non-overlapping matches
//         let matches = vec![
//             create_match(
//                 0,
//                 10,
//                 10,
//                 20,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 5,
//                 3,
//                 "test1",
//                 Strand::Fwd,
//             ),
//             create_match(
//                 0,
//                 10,
//                 30,
//                 40,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 4,
//                 2,
//                 "test2",
//                 Strand::Fwd,
//             ),
//         ];
//         let result = collapse_overlapping_matches(&matches, 0.5);
//         assert_eq!(result.len(), 2);
//         assert_eq!(result[0].label, "test1");
//         assert_eq!(result[1].label, "test2");

//         // Test overlapping matches - should collapse to best one
//         let matches = vec![
//             create_match(
//                 0,
//                 10,
//                 10,
//                 20,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 5,
//                 3,
//                 "test1",
//                 Strand::Fwd,
//             ),
//             create_match(
//                 0,
//                 10,
//                 15,
//                 25,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 4,
//                 2,
//                 "test2",
//                 Strand::Fwd,
//             ), // Better cost
//         ];
//         let result = collapse_overlapping_matches(&matches, 0.5);
//         assert_eq!(result.len(), 1);
//         assert_eq!(result[0].label, "test2"); // Should select the one with lower cost

//         // Test multiple overlapping groups
//         let matches = vec![
//             create_match(
//                 0,
//                 10,
//                 10,
//                 20,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 5,
//                 3,
//                 "test1",
//                 Strand::Fwd,
//             ),
//             create_match(
//                 0,
//                 10,
//                 15,
//                 25,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 4,
//                 2,
//                 "test2",
//                 Strand::Fwd,
//             ),
//             create_match(
//                 0,
//                 10,
//                 40,
//                 50,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 3,
//                 1,
//                 "test3",
//                 Strand::Fwd,
//             ),
//             create_match(
//                 0,
//                 10,
//                 45,
//                 55,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 6,
//                 4,
//                 "test4",
//                 Strand::Fwd,
//             ),
//         ];
//         let result = collapse_overlapping_matches(&matches, 0.5);
//         assert_eq!(result.len(), 2);
//         assert_eq!(result[0].label, "test2"); // Best from first group
//         assert_eq!(result[1].label, "test3"); // Best from second group
//     }

//     #[test]
//     fn test_is_overlap() {
//         let match1 = create_match(
//             0,
//             10,
//             10,
//             20,
//             0,
//             10,
//             BarcodeType::Fbar,
//             5,
//             3,
//             "test1",
//             Strand::Fwd,
//         );
//         let match2 = create_match(
//             0,
//             10,
//             15,
//             25,
//             0,
//             10,
//             BarcodeType::Fbar,
//             4,
//             2,
//             "test2",
//             Strand::Fwd,
//         );

//         // Test overlapping matches
//         assert!(is_overlap(&match1, &match2, 0.5));
//         assert!(is_overlap(&match2, &match1, 0.5));

//         // Test non-overlapping matches
//         let match3 = create_match(
//             0,
//             10,
//             30,
//             40,
//             0,
//             10,
//             BarcodeType::Fbar,
//             3,
//             1,
//             "test3",
//             Strand::Fwd,
//         );
//         assert!(!is_overlap(&match1, &match3, 0.5));
//         assert!(!is_overlap(&match3, &match1, 0.5));

//         // Test adjacent matches (no overlap)
//         let match4 = create_match(
//             0,
//             10,
//             20,
//             30,
//             0,
//             10,
//             BarcodeType::Fbar,
//             2,
//             1,
//             "test4",
//             Strand::Fwd,
//         );
//         assert!(!is_overlap(&match1, &match4, 0.5));

//         // Test threshold behavior
//         let match5 = create_match(
//             0,
//             10,
//             18,
//             22,
//             0,
//             10,
//             BarcodeType::Fbar,
//             1,
//             1,
//             "test5",
//             Strand::Fwd,
//         );
//         // Overlap is 2, min_len is 4, so overlap ratio is 0.5
//         assert!(is_overlap(&match1, &match5, 0.5));
//         assert!(!is_overlap(&match1, &match5, 0.6)); // Higher threshold
//     }

//     #[test]
//     fn test_select_best_match() {
//         // Test single match
//         let matches = vec![create_match(
//             0,
//             10,
//             10,
//             20,
//             0,
//             10,
//             BarcodeType::Fbar,
//             5,
//             3,
//             "test1",
//             Strand::Fwd,
//         )];
//         let result = select_best_match(matches);
//         assert_eq!(result.label, "test1");

//         // Test multiple matches with different costs
//         let matches = vec![
//             create_match(
//                 0,
//                 10,
//                 10,
//                 20,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 5,
//                 3,
//                 "test1",
//                 Strand::Fwd,
//             ),
//             create_match(
//                 0,
//                 10,
//                 10,
//                 20,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 4,
//                 2,
//                 "test2",
//                 Strand::Fwd,
//             ), // Lower cost
//             create_match(
//                 0,
//                 10,
//                 10,
//                 20,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 6,
//                 4,
//                 "test3",
//                 Strand::Fwd,
//             ),
//         ];
//         let result = select_best_match(matches);
//         assert_eq!(result.label, "test2"); // Should select the one with lowest total cost

//         // Test matches with different types (Fbar < Rbar in ordering)
//         let matches = vec![
//             create_match(
//                 0,
//                 10,
//                 10,
//                 20,
//                 0,
//                 10,
//                 BarcodeType::Rbar,
//                 1,
//                 1,
//                 "test1",
//                 Strand::Fwd,
//             ),
//             create_match(
//                 0,
//                 10,
//                 10,
//                 20,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 5,
//                 5,
//                 "test2",
//                 Strand::Fwd,
//             ), // Lower type
//         ];
//         let result = select_best_match(matches);
//         assert_eq!(result.label, "test2"); // Should select Fbar over Rbar

//         // Test matches with same type and cost - should select first one
//         let matches = vec![
//             create_match(
//                 0,
//                 10,
//                 10,
//                 20,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 5,
//                 3,
//                 "test1",
//                 Strand::Fwd,
//             ),
//             create_match(
//                 0,
//                 10,
//                 10,
//                 20,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 5,
//                 3,
//                 "test2",
//                 Strand::Fwd,
//             ),
//         ];
//         let result = select_best_match(matches);
//         assert_eq!(result.label, "test1"); // Should select first one when tied
//     }

//     #[test]
//     fn test_edge_cases() {
//         // Test zero-length matches
//         let match1 = create_match(
//             0,
//             0,
//             10,
//             10,
//             0,
//             0,
//             BarcodeType::Fbar,
//             5,
//             3,
//             "test1",
//             Strand::Fwd,
//         );
//         let match2 = create_match(
//             0,
//             0,
//             10,
//             10,
//             0,
//             0,
//             BarcodeType::Fbar,
//             4,
//             2,
//             "test2",
//             Strand::Fwd,
//         );
//         assert!(!is_overlap(&match1, &match2, 0.5)); // No overlap for zero-length

//         // Test matches with very small overlap
//         let match1 = create_match(
//             0,
//             2,
//             10,
//             12,
//             0,
//             2,
//             BarcodeType::Fbar,
//             5,
//             3,
//             "test1",
//             Strand::Fwd,
//         );
//         let match2 = create_match(
//             0,
//             2,
//             11,
//             13,
//             0,
//             2,
//             BarcodeType::Fbar,
//             4,
//             2,
//             "test2",
//             Strand::Fwd,
//         );
//         // Overlap is 1, min_len is 2, so overlap ratio is 0.5
//         assert!(is_overlap(&match1, &match2, 0.5));
//         assert!(!is_overlap(&match1, &match2, 0.6)); // Higher threshold

//         // Test matches with different strands
//         let match1 = create_match(
//             0,
//             10,
//             10,
//             20,
//             0,
//             10,
//             BarcodeType::Fbar,
//             5,
//             3,
//             "test1",
//             Strand::Fwd,
//         );
//         let match2 = create_match(
//             0,
//             10,
//             15,
//             25,
//             0,
//             10,
//             BarcodeType::Fbar,
//             4,
//             2,
//             "test2",
//             Strand::Rc,
//         );
//         // Overlap detection should work regardless of strand
//         assert!(is_overlap(&match1, &match2, 0.5));
//     }

//     #[test]
//     fn test_sorting_behavior() {
//         // Test that matches are sorted by read_start before processing
//         let matches = vec![
//             create_match(
//                 0,
//                 10,
//                 30,
//                 40,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 5,
//                 3,
//                 "test3",
//                 Strand::Fwd,
//             ),
//             create_match(
//                 0,
//                 10,
//                 10,
//                 20,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 4,
//                 2,
//                 "test1",
//                 Strand::Fwd,
//             ),
//             create_match(
//                 0,
//                 10,
//                 20,
//                 30,
//                 0,
//                 10,
//                 BarcodeType::Fbar,
//                 3,
//                 1,
//                 "test2",
//                 Strand::Fwd,
//             ),
//         ];
//         let result = collapse_overlapping_matches(&matches, 0.5);
//         // Should still return 3 matches since they don't overlap
//         assert_eq!(result.len(), 3);
//         // Order should be preserved based on read_start
//         assert_eq!(result[0].label, "test1");
//         assert_eq!(result[1].label, "test2");
//         assert_eq!(result[2].label, "test3");
//     }
// }
