use crate::annotate::barcodes::BarcodeType;
use crate::annotate::searcher::BarbellMatch;
use sassy::Strand;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;

#[derive(Debug, Eq, PartialEq, Clone, Serialize, Deserialize)]
pub enum CutDirection {
    Before, // < cut at match.start
    After,  // > cut at match.end
}

// Records where..cut, in what direction, and also with id it belongs in case of paired cuts
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Cut {
    pub group_id: usize,
    pub direction: CutDirection,
}

#[derive(Debug, Eq, PartialEq)]
pub struct PatternElement {
    pub match_type: BarcodeType,
    pub orientation: Option<Strand>,
    pub label: Option<String>,      // Only used for Barcode matches
    pub placeholder: Option<usize>, // For referencing previous matches like $1, $2, etc.
    pub range: (isize, isize),      // Minimum and maximum positions (distance)
    pub relative_to: Option<RelativePosition>, // Relative..left or right ends, or previous match
    pub cuts: Option<Vec<Cut>>,     // Multiple cuts with identifiers
}

// Relative position of matches
#[derive(Debug, Eq, PartialEq)]
pub enum RelativePosition {
    Left,
    Right,
    PrevLeft,
    PrevRight, // unimplemented for now
}

#[derive(Debug, Eq, PartialEq)]
pub struct Pattern {
    pub elements: Vec<PatternElement>,
}

impl Cut {
    pub fn new(group_id: usize, direction: CutDirection) -> Self {
        Self {
            group_id,
            direction,
        }
    }

    pub fn from_string(s: &str) -> Option<Self> {
        let s = s.trim();
        if s.starts_with("Before(") && s.ends_with(")") {
            let content = &s[7..s.len() - 1];
            let group_id = content.parse().ok()?;
            Some(Cut::new(group_id, CutDirection::Before))
        } else if s.starts_with("After(") && s.ends_with(")") {
            let content = &s[6..s.len() - 1];
            let group_id = content.parse().ok()?;
            Some(Cut::new(group_id, CutDirection::After))
        } else {
            None
        }
    }

    pub fn from_pattern_string(pat_str: &str) -> Option<Self> {
        let two_char_prefix = &pat_str[..2];

        // Default..id = 0 if users does not specify it
        let id = if pat_str.len() == 2 {
            0
        } else {
            pat_str[2..].parse::<usize>().ok()?
        };

        match two_char_prefix {
            ">>" => Some(Cut::new(id, CutDirection::After)),
            "<<" => Some(Cut::new(id, CutDirection::Before)),
            _ => None,
        }
    }
}

impl fmt::Display for Cut {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.direction {
            CutDirection::Before => write!(f, "Before({})", self.group_id),
            CutDirection::After => write!(f, "After({})", self.group_id),
        }
    }
}

fn check_match_type_and_label(m: &BarbellMatch, pattern_element: &PatternElement) -> bool {
    // println!("Checking match type: {:?} against pattern type: {:?}",
    //     m.match_str.match_type, pattern_element.match_type);

    // First check if match types are exactly equal
    if m.match_type != pattern_element.match_type {
        // println!(
        //     "Match type mismatch: {:?} != {:?}",
        //     m.match_str.match_type, pattern_element.match_type
        // );
        return false;
    }

    match pattern_element.match_type {
        // If Ftag/Rtag, and we have a label in the pattern
        // the pattern label should match, otherwise we assume pattern
        // is label *, and can match any label
        BarcodeType::Ftag | BarcodeType::Rtag => {
            if let Some(ref expected_label) = pattern_element.label {
                // Labels are not the same
                //if let Some(ref m_label) = m.label {
                if expected_label.starts_with("~") {
                    // we do substring check
                    if let Some(substring) = expected_label.strip_prefix('~')
                        && !m.label.contains(substring)
                    {
                        // println!(
                        //     "Substring label mismatch: {:?} not in {:?}",
                        //     expected_label, &m.label
                        // );
                        return false;
                    }
                } else if expected_label != &m.label {
                    // println!("no substring prefix");
                    // println!(
                    //     "Explicit label mismatch: {:?} != {:?}",
                    //     expected_label, &m.label
                    // );
                    return false;
                }
            }
            true
        }
        // Flanks are always ok
        BarcodeType::Fflank | BarcodeType::Rflank => true,
    }
}

fn check_placeholder(
    m: &BarbellMatch,
    pattern_element: &PatternElement,
    matched_labels: &mut HashMap<usize, String>,
) -> bool {
    if let Some(placeholder_key) = pattern_element.placeholder {
        if let Some(stored_label) = matched_labels.get(&placeholder_key) {
            // Check against stored label
            if &m.label != stored_label {
                // println!(
                //     "Placeholder label mismatch: current {:?} != stored {:?}",
                //     &m.label, stored_label
                // );
                return false;
            }
        } else {
            // Store new label
            matched_labels.insert(placeholder_key, m.label.clone());
        }
    }
    true
}

fn check_orientation(m: &BarbellMatch, pattern_element: &PatternElement) -> bool {
    pattern_element.orientation.is_none() || pattern_element.orientation.as_ref() == Some(&m.strand)
}

fn check_relative_position(
    m: &BarbellMatch,
    pattern_element: &PatternElement,
    prev_end: Option<isize>,
    seq_len: isize,
) -> bool {
    if let Some(relative_to) = &pattern_element.relative_to {
        let m_start = m.read_start_bar as isize;
        let m_end = m.read_end_bar as isize;
        match relative_to {
            RelativePosition::Left => {
                let (left_bound, right_bound) = pattern_element.range;
                if m_start < left_bound || m_start > right_bound {
                    // println!(
                    //     "Relative Left mismatch: m.start {} not in [{}, {}]",
                    //     m_start, left_bound, right_bound
                    // );
                    return false;
                }
            }
            RelativePosition::Right => {
                let left_bound = seq_len.saturating_sub(pattern_element.range.1);
                let right_bound = seq_len.saturating_sub(pattern_element.range.0);
                if m_end < left_bound || m_end > right_bound {
                    return false;
                }
            }
            RelativePosition::PrevLeft => {
                if let Some(prev_end_pos) = prev_end {
                    let left_bound = prev_end_pos.saturating_add(pattern_element.range.0);
                    let right_bound = prev_end_pos.saturating_add(pattern_element.range.1);
                    if m_start < left_bound || m_start > right_bound {
                        for _ in 0..10 {
                            println!(
                                "PrevLeft gap mismatch: m.start {} not in [{}, {}]",
                                m_start, left_bound, right_bound
                            );
                        }

                        return false;
                    }
                }
            }
            RelativePosition::PrevRight => {
                unimplemented!(
                    "This does not make much sense when going left > right, unless we DP optimal"
                );
            }
        }
    }
    true
}

fn matches_pattern_element(
    m: &BarbellMatch,
    pattern_element: &PatternElement,
    prev_end: Option<isize>,
    matched_labels: &mut HashMap<usize, String>,
    seq_len: usize,
) -> bool {
    check_match_type_and_label(m, pattern_element)
        && check_placeholder(m, pattern_element, matched_labels)
        && check_orientation(m, pattern_element)
        && check_relative_position(m, pattern_element, prev_end, seq_len as isize)
}

pub fn match_pattern(matches: &[BarbellMatch], pattern: &Pattern) -> (bool, Vec<(usize, Cut)>) {
    let mut prev_end: Option<isize> = None;
    let mut matched_labels: HashMap<usize, String> = HashMap::new();
    let mut current_match_idx = 0;
    let mut cut_positions: Vec<(usize, Cut)> = Vec::new();

    // Early return if not enough matches
    if matches.len() < pattern.elements.len() {
        return (false, vec![]);
    }

    for pattern_element in pattern.elements.iter() {
        if current_match_idx >= matches.len() {
            return (false, vec![]);
        }

        let m = &matches[current_match_idx];
        let seq_len = m.read_len;

        if matches_pattern_element(m, pattern_element, prev_end, &mut matched_labels, seq_len) {
            // If this element has a cut marker, record the position and direction
            if let Some(cuts) = &pattern_element.cuts {
                for cut in cuts {
                    cut_positions.push((current_match_idx, cut.clone()));
                }
            }

            prev_end = Some(m.read_end_bar as isize);
            current_match_idx += 1;
        } else {
            return (false, vec![]);
        }
    }

    (true, cut_positions)
}

#[macro_export]
macro_rules! pattern_from_str {
    ($pattern:expr) => {{
        use sassy::Strand;
        use $crate::annotate::barcodes::BarcodeType;
        use $crate::filter::pattern::*;

        fn parse_range(range_str: &str) -> Option<(isize, isize)> {
            let parts: Vec<&str> = range_str
                .trim_matches(|p| p == '(' || p == ')')
                .split("..")
                .collect();
            if parts.len() == 2 {
                let start = parts[0].trim().parse::<isize>().ok()?;
                let end = parts[1].trim().parse::<isize>().ok()?;
                Some((start, end))
            } else {
                None
            }
        }

        fn parse_position(pos_str: &str) -> Option<(RelativePosition, (isize, isize))> {
            let pos_parts: Vec<&str> = pos_str.split('(').collect();

            if pos_parts.len() != 2 {
                return None;
            }

            let position = match pos_parts[0].trim_start_matches('@') {
                "left" => RelativePosition::Left,
                "right" => RelativePosition::Right,
                "prev_left" => RelativePosition::PrevLeft,
                _ => return None,
            };

            let range = parse_range(&pos_str[pos_parts[0].len()..].trim())?;
            Some((position, range))
        }

        fn parse_element(element_str: &str) -> Option<PatternElement> {
            // println!("element_str: {:?}", element_str);
            // First split the string into type and parameters
            let parts: Vec<&str> = element_str.splitn(2, '[').collect();
            if parts.len() != 2 {
                return None;
            }

            // Parse the match type
            let match_type = match parts[0].trim() {
                "Ftag" => BarcodeType::Ftag,
                "Rtag" => BarcodeType::Rtag,
                "Fflank" => BarcodeType::Fflank,
                "Rflank" => BarcodeType::Rflank,
                "Flank" | "flank" => panic!("Flank is not valid, use Fflank or Rflank"),
                _ => return None,
            };

            // Initialize defaults
            let mut label = None;
            let mut placeholder = None;
            let mut orientation = None;
            let mut relative_to = None;
            let mut range = (0, 0);
            let mut cuts: Vec<Cut> = vec![];

            // Parse parameters
            let params = parts[1].trim_end_matches(']');
            //     println!("params: {:?}", params);
            for param in params.split(',').map(|s| s.trim()) {
                //    println!("\tparam: {:?}", param);
                match param {
                    // Orientation
                    "fw" => orientation = Some(Strand::Fwd),
                    "rc" => orientation = Some(Strand::Rc),

                    // Position and range
                    p if p.starts_with('@') => {
                        if let Some((pos, r)) = parse_position(p) {
                            relative_to = Some(pos);
                            range = r;
                        }
                    }

                    // Placeholder
                    p if p.starts_with('?') => {
                        if let Ok(num) = p[1..].parse::<usize>() {
                            placeholder = Some(num);
                        }
                    }

                    p if p.starts_with('>') || p.starts_with('<') => {
                        if let Some(cut) = Cut::from_pattern_string(p) {
                            cuts.push(cut);
                        }
                    }

                    // Label
                    "*" => (), // Any label
                    p => label = Some(p.trim_matches('"').to_string()),
                }
            }

            // If cuts are empty, just None it
            let cuts: Option<Vec<Cut>> = if cuts.is_empty() { None } else { Some(cuts) };

            Some(PatternElement {
                match_type,
                orientation,
                label,
                placeholder,
                range,
                relative_to,
                cuts,
            })
        }

        // Split on __ and parse each element
        let elements: Vec<PatternElement> = $pattern
            .split("__")
            .filter_map(|s| parse_element(s.trim()))
            .collect();

        Pattern { elements }
    }};
}

mod tests {

    use super::*;

    #[test]
    fn test_pattern_macro() {
        let pattern = pattern_from_str!(
            "Ftag[fw, *, @left(0..250)]__Fflank[fw, @prev_left(5..100)]__Rtag[?1, fw, @right(0..20)]"
        );
        assert_eq!(pattern.elements.len(), 3);
        assert_eq!(
            pattern,
            Pattern {
                elements: vec![
                    PatternElement {
                        match_type: BarcodeType::Ftag,
                        orientation: Some(Strand::Fwd),
                        label: None,
                        placeholder: None,
                        range: (0, 250),
                        relative_to: Some(RelativePosition::Left),
                        cuts: None
                    },
                    PatternElement {
                        match_type: BarcodeType::Fflank,
                        orientation: Some(Strand::Fwd),
                        label: None,
                        placeholder: None,
                        range: (5, 100),
                        relative_to: Some(RelativePosition::PrevLeft),
                        cuts: None,
                    },
                    PatternElement {
                        match_type: BarcodeType::Rtag,
                        orientation: Some(Strand::Fwd),
                        label: None,
                        placeholder: Some(1),
                        range: (0, 20),
                        relative_to: Some(RelativePosition::Right),
                        cuts: None
                    }
                ]
            }
        );
    }

    #[test]
    fn test_distance_to_left_end() {
        // let pattern = pattern_from_str!("Ftag[fw, *, @left(0-250)]");
        let pattern = pattern_from_str!("Ftag[fw, *, @left(0..250)]");

        let mut matches = vec![BarbellMatch::new(
            0,   // read_start_bar
            100, // read_end_bar
            0,   // read_start_flank
            100, // read_end_flank
            0,   // bar_start
            24,  // bar_end
            BarcodeType::Ftag,
            0, // flank_cost
            0, // barcode_cost
            "XXX".to_string(),
            Strand::Fwd,
            500, // read_len
            "test".to_string(),
            0,
            None,
        )];

        matches[0].read_start_bar = 0;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(is_match);

        matches[0].read_start_bar = 100;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(is_match);

        matches[0].read_start_bar = 250;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(is_match);

        matches[0].read_start_bar = 251;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(!is_match);
    }

    #[test]
    fn test_distance_to_right_end() {
        let pattern = pattern_from_str!("Ftag[fw, *, @right(0..250)]");

        let mut matches = vec![BarbellMatch::new(
            0,   // read_start_bar
            100, // read_end_bar
            0,   // read_start_flank
            100, // read_end_flank
            0,   // bar_start
            100, // bar_end
            BarcodeType::Ftag,
            0, // flank_cost
            0, // barcode_cost
            "XXX".to_string(),
            Strand::Fwd,
            500, // read_len
            "test".to_string(),
            0,
            None,
        )];

        matches[0].read_end_bar = 500;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(is_match);

        matches[0].read_end_bar = 450;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(is_match);

        matches[0].read_end_bar = 250;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(is_match);

        matches[0].read_end_bar = 249;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(!is_match);
    }

    #[test]
    fn test_distance_to_prev_left() {
        let pattern =
            pattern_from_str!("Ftag[fw, *, @left(0..250)]__Fflank[fw, @prev_left(5..100)]");
        println!("Pattern: {:?}", pattern);

        let mut matches = vec![
            BarbellMatch::new(
                0,   // read_start_bar
                100, // read_end_bar
                0,   // read_start_flank
                100, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                500, // read_len
                "test".to_string(),
                0,
                None,
            ),
            BarbellMatch::new(
                100, // read_start_bar
                200, // read_end_bar
                100, // read_start_flank
                200, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Fflank,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                500, // read_len
                "test".to_string(),
                0,
                None,
            ),
        ];

        // Should not pass,
        matches[1].read_start_bar = 50;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(!is_match);

        // Same start as previous match, does not satisfy min distance of 5
        matches[1].read_start_bar = 100;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(!is_match);

        // Should pass, just satisfies min distance of 5
        matches[1].read_start_bar = 105;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(is_match);

        // Should still pass, but maxing out
        matches[1].read_start_bar = 200;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(is_match);

        // Should not pass anymore, too far away
        matches[1].read_start_bar = 201;
        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(!is_match);
    }

    #[test]
    fn test_placeholder() {
        // Means any label for front, but the last label should be the same as the first label
        let pattern =
            pattern_from_str!("Ftag[fw, ?1, @left(0..250)]__Rtag[fw, ?1, @right(0..250)]");
        println!("Pattern: {:?}", pattern);
        let mut matches = vec![
            BarbellMatch::new(
                0,   // read_start_bar
                100, // read_end_bar
                0,   // read_start_flank
                100, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                250, // read_len
                "test".to_string(),
                0,
                None,
            ),
            BarbellMatch::new(
                100, // read_start_bar
                200, // read_end_bar
                100, // read_start_flank
                200, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Rtag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                250, // read_len
                "test".to_string(),
                0,
                None,
            ),
        ];

        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(is_match);

        // Edit..different labels
        matches[1].label = "yyyy".to_string();

        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(!is_match);
    }

    #[test]
    fn test_placeholder_mixed_labels() {
        // Means any label for front, but the last label should be the same as the first label
        let pattern =
            pattern_from_str!("Ftag[fw, ?1, @left(0..250)]__Rtag[fw, ?2, @right(0..250)]");
        println!("Pattern: {:?}", pattern);
        let mut matches = vec![
            BarbellMatch::new(
                0,   // read_start_bar
                100, // read_end_bar
                0,   // read_start_flank
                100, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                250, // read_len
                "test".to_string(),
                0,
                None,
            ),
            BarbellMatch::new(
                100, // read_start_bar
                200, // read_end_bar
                100, // read_start_flank
                200, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Rtag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                250, // read_len
                "test".to_string(),
                0,
                None,
            ),
        ];

        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(is_match);
    }

    #[test]
    fn test_placeholder_not_ordered() {
        // Means any label for front, but the last label should be the same as the first label
        let pattern = pattern_from_str!(
            "Ftag[fw, ?1, @left(0..250)]__Ftag[fw, ?2, @prev_left(0..250)]__Ftag[fw, ?1, @left(0..250)]"
        );
        println!("Pattern: {:?}", pattern);
        let mut matches = vec![
            BarbellMatch::new(
                0,   // read_start_bar
                100, // read_end_bar
                0,   // read_start_flank
                100, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                600, // read_len
                "test".to_string(),
                0,
                None,
            ),
            BarbellMatch::new(
                100, // read_start_bar
                200, // read_end_bar
                100, // read_start_flank
                200, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "YYY".to_string(),
                Strand::Fwd,
                600, // read_len
                "test".to_string(),
                0,
                None,
            ),
            BarbellMatch::new(
                100, // read_start_bar
                200, // read_end_bar
                550, // read_start_flank
                600, // read_end_flank
                0,   // bar_start
                100, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                600, // read_len
                "test".to_string(),
                0,
                None,
            ),
        ];

        let (is_match, _) = match_pattern(&matches, &pattern);
        assert!(is_match);
    }

    #[test]
    fn test_pattern_with_cuts_default_fallback() {
        let pattern =
            pattern_from_str!("Ftag[fw, *, >>, @left(0..250)]__Fflank[fw, <<, @prev_left(5..100)]");

        let matches = vec![
            BarbellMatch::new(
                0,  // read_start_bar
                10, // read_end_bar
                0,  // read_start_flank
                10, // read_end_flank
                0,  // bar_start
                10, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                250, // read_len
                "test".to_string(),
                0,
                None,
            ),
            BarbellMatch::new(
                15, // read_start_bar
                20, // read_end_bar
                15, // read_start_flank
                20, // read_end_flank
                0,  // bar_start
                5,  // bar_end
                BarcodeType::Fflank,
                0, // flank_cost
                0, // barcode_cost
                "@Nothing".to_string(),
                Strand::Fwd,
                250, // read_len
                "test".to_string(),
                0,
                None,
            ),
        ];

        let (is_match, cut_positions) = match_pattern(&matches, &pattern);
        assert!(is_match);
        assert_eq!(
            cut_positions,
            vec![
                (0, Cut::new(0, CutDirection::After)),
                (1, Cut::new(0, CutDirection::Before))
            ]
        );
    }

    #[test]
    fn test_pattern_with_cuts_single_group() {
        let pattern = pattern_from_str!(
            "Ftag[fw, *, >>1, @left(0..250)]__Fflank[fw, <<1, @prev_left(5..100)]"
        );

        let matches = vec![
            BarbellMatch::new(
                0,  // read_start_bar
                10, // read_end_bar
                0,  // read_start_flank
                10, // read_end_flank
                0,  // bar_start
                10, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                250, // read_len
                "test".to_string(),
                0,
                None,
            ),
            BarbellMatch::new(
                15, // read_start_bar
                20, // read_end_bar
                15, // read_start_flank
                20, // read_end_flank
                0,  // bar_start
                5,  // bar_end
                BarcodeType::Fflank,
                0, // flank_cost
                0, // barcode_cost
                "@Nothing".to_string(),
                Strand::Fwd,
                250, // read_len
                "test".to_string(),
                0,
                None,
            ),
        ];

        let (is_match, cut_positions) = match_pattern(&matches, &pattern);
        assert!(is_match);
        assert_eq!(
            cut_positions,
            vec![
                (0, Cut::new(1, CutDirection::After)),
                (1, Cut::new(1, CutDirection::Before))
            ]
        );
    }

    #[test]
    fn test_pattern_with_multiple_cuts_fallback() {
        let pattern = pattern_from_str!(
            "Ftag[fw, *, >>1, @left(0..250)]__Fflank[fw, <<1, @prev_left(5..100)]__Rtag[fw, *, <<2, @right(0..20)]"
        );

        let matches = vec![
            BarbellMatch::new(
                0,  // read_start_bar
                10, // read_end_bar
                0,  // read_start_flank
                10, // read_end_flank
                0,  // bar_start
                10, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "XXX".to_string(),
                Strand::Fwd,
                50, // read_len
                "test".to_string(),
                0,
                None,
            ),
            BarbellMatch::new(
                15, // read_start_bar
                20, // read_end_bar
                15, // read_start_flank
                20, // read_end_flank
                0,  // bar_start
                5,  // bar_end
                BarcodeType::Fflank,
                0, // flank_cost
                0, // barcode_cost
                "@Nothing".to_string(),
                Strand::Fwd,
                50, // read_len
                "test".to_string(),
                0,
                None,
            ),
            BarbellMatch::new(
                30, // read_start_bar
                40, // read_end_bar
                30, // read_start_flank
                40, // read_end_flank
                0,  // bar_start
                10, // bar_end
                BarcodeType::Rtag,
                0, // flank_cost
                0, // barcode_cost
                "YYY".to_string(),
                Strand::Fwd,
                50, // read_len
                "test".to_string(),
                0,
                None,
            ),
        ];

        let (is_match, cut_positions) = match_pattern(&matches, &pattern);
        assert!(is_match);
        assert_eq!(
            cut_positions,
            vec![
                (0, Cut::new(1, CutDirection::After)),
                (1, Cut::new(1, CutDirection::Before)),
                (2, Cut::new(2, CutDirection::Before))
            ]
        );
    }

    #[test]
    fn test_cut_from_string() {
        let cut_str = "After(1)";
        let cut = Cut::from_string(cut_str).unwrap();
        assert_eq!(cut, Cut::new(1, CutDirection::After));

        let cut_str = "Before(2)";
        let cut = Cut::from_string(cut_str).unwrap();
        assert_eq!(cut, Cut::new(2, CutDirection::Before));

        // Test invalid input
        assert!(Cut::from_string("Invalid").is_none());
        assert!(Cut::from_string("After(abc)").is_none());
    }
}
