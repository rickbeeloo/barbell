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

impl CutDirection {
    pub fn rotate(&self) -> Self {
        match self {
            CutDirection::Before => CutDirection::After,
            CutDirection::After => CutDirection::Before,
        }
    }
}

// Records where..cut, in what direction, and also with id it belongs in case of paired cuts
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Cut {
    pub group_id: usize,
    pub direction: CutDirection,
}

#[derive(Debug, Clone, Eq, PartialEq)]
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
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum RelativePosition {
    Left,
    Right,
    PrevLeft,
    PrevRight, // unimplemented for now
    Any,
}

#[derive(Debug, Clone, Eq, PartialEq)]
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
    // First check if match types are exactly equal
    if m.match_type != pattern_element.match_type {
        return false;
    }

    match pattern_element.match_type {
        // If Ftag/Rtag, and we have a label in the pattern
        // the pattern label should match, otherwise we assume pattern
        // is label *, and can match any label
        BarcodeType::Ftag | BarcodeType::Rtag | BarcodeType::Fbar | BarcodeType::Rbar => {
            if let Some(ref expected_label) = pattern_element.label {
                // Labels are not the same
                if expected_label.starts_with("~") {
                    // we do substring check
                    if let Some(substring) = expected_label.strip_prefix('~')
                        && !m.label.contains(substring)
                    {
                        return false;
                    }
                } else if expected_label != &m.label {
                    return false;
                }
            }
            true
        }
        // Flanks are always ok
        BarcodeType::Fflank
        | BarcodeType::Rflank
        | BarcodeType::Fadapter
        | BarcodeType::Radapter => true,
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
                        return false;
                    }
                }
            }
            RelativePosition::PrevRight => {
                unimplemented!(
                    "This does not make much sense when going left > right, unless we DP optimal"
                );
            }
            RelativePosition::Any => {
                return true;
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

pub fn find_seeds(matches: &[BarbellMatch], pattern: &Pattern) -> Vec<usize> {
    // To match a pattern as subpattern, we can first search the first element
    // in the pattern and locate
    let mut seeds: Vec<usize> = Vec::new();
    let seed_element = pattern.elements.first().unwrap();

    for (idx, m) in matches.iter().enumerate() {
        if matches_pattern_element(m, seed_element, None, &mut HashMap::new(), m.read_len) {
            seeds.push(idx);
        }
    }
    seeds
}

pub fn check_from_seed(
    matches: &[BarbellMatch],
    pattern: &Pattern,
    seed_idx: usize,
) -> (bool, Vec<(usize, Cut)>) {
    // let matches = &matches[seed_idx..];

    // Early return if not enough matches
    if matches.len() < pattern.elements.len() {
        return (false, vec![]);
    }

    let mut prev_end: Option<isize> = None;
    let mut matched_labels: HashMap<usize, String> = HashMap::new();
    let mut current_match_idx = seed_idx;

    let mut cut_positions = Vec::with_capacity(2);

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

    // If the user specified only one cut positions (i.e. Ftag,>>) we cut until the
    // end of the read UNLESS there is (any) other match elements after it
    // If the user specified only one cut positions the other way (<<) we search
    // in front if there is anything left
    // if we already have two cuts we are fine and can return instantly
    if cut_positions.len() == 2 {
        // Bases covered here is between the cuts
        return (true, cut_positions);
    }

    // We have a single cut, time to check
    let cut = cut_positions.last().unwrap();
    let cut_idx = cut.0;

    let cut_direction = cut.1.direction.clone();
    match cut_direction {
        CutDirection::After => {
            // We cut until the end of the read UNLESS there is (any) other match elements after it
            if cut_idx < matches.len() - 1 {
                // We just get the next index
                let flip_cut = cut_direction.rotate(); // Have to flip cut sign from << to >>, or vice versa
                let new_cut = Cut::new(cut.1.group_id, flip_cut);
                cut_positions.push((cut_idx + 1, new_cut));
            }
        }
        CutDirection::Before => {
            if cut_idx > 0 {
                let flip_cut = cut_direction.rotate(); // Have to flip cut sign from << to >>, or vice versa
                let new_cut = Cut::new(cut.1.group_id, flip_cut);
                cut_positions.push((cut_idx - 1, new_cut));
            }
        }
    }

    (true, cut_positions)
}

pub fn match_sub_pattern(
    matches: &[BarbellMatch],
    pattern: &Pattern,
) -> (bool, Vec<Vec<(usize, Cut)>>) {
    // Early return if not enough matches
    if matches.len() < pattern.elements.len() {
        return (false, vec![]);
    }

    // First we look in the matches if any equals that of the beginnign of the pattern
    let seeds = find_seeds(matches, pattern);

    // Then from each seed we just check if everything "after" it matches the pattern
    let mut all_cuts = vec![];
    for seed_idx in seeds {
        let (is_match, cut_positions) = check_from_seed(matches, pattern, seed_idx);
        if is_match {
            all_cuts.push(cut_positions.clone());
        }
    }

    (!all_cuts.is_empty(), all_cuts)
}

pub fn match_full_pattern(
    matches: &[BarbellMatch],
    pattern: &Pattern,
) -> (bool, Vec<(usize, Cut)>) {
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
            // 1. If we don't match anymore we should abort
            return (false, vec![]);
        }
    }

    if current_match_idx < matches.len() {
        // 2. If we did not exhaust all matches by matching them
        // to the pattern, we should abort
        return (false, vec![]);
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
                "any" => RelativePosition::Any,
                _ => return None,
            };

            let range = parse_range(&pos_str[pos_parts[0].len()..].trim())?;
            Some((position, range))
        }

        fn basic_verify(pattern: &[PatternElement], patter_str: &str) -> bool {
            //Fixme: more complex verification would be nice
            let user_elems = patter_str.matches("__").count() + 1;
            user_elems == pattern.len()
        }

        fn parse_element(element_str: &str) -> Option<PatternElement> {
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
        let elements: Vec<PatternElement> = {
            let pattern_str: &str = $pattern;
            pattern_str
                .split("__")
                .filter_map(|s| parse_element(s.trim()))
                .collect()
        };

        // Do basic verification
        if !basic_verify(&elements, $pattern) {
            eprintln!("Seems we could not convert all your pattern elements to Barbell patterns, please compare your string to:");
            for (i, el) in elements.iter().enumerate() {
                eprintln!("  Element {}: {:#?}", i, el);
            }
            panic!("Pattern parse error for: {:?}", $pattern);
        }

        Pattern { elements }
    }};
}

#[cfg(test)]
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
    fn test_double_sub_pattern_matches() {
        let pattern = pattern_from_str!("Ftag[fw, *, @left(0..250), >>]");
        let matches = vec![
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
        ];
        let seeds = find_seeds(&matches, &pattern);
        assert_eq!(seeds, vec![0, 1]);
        // Using full pattern search this should not match
        let (is_match, _cut_positions) = match_full_pattern(&matches, &pattern);
        assert!(!is_match);

        // Using sub pattern search this should match twice
        // the firt slice from match1 --- match2, and match2-- end of read
        // First case:
        let first_slice = [
            (
                0,
                Cut {
                    group_id: 0,
                    direction: CutDirection::After,
                },
            ),
            (
                1,
                Cut {
                    group_id: 0,
                    direction: CutDirection::Before,
                },
            ),
        ]
        .to_vec();
        // Second case:
        let second_slice = [(
            1,
            Cut {
                group_id: 0,
                direction: CutDirection::After,
            },
        )]
        .to_vec();
        // Total exepcted cut positions:
        let expected_cut_positions = vec![first_slice, second_slice];
        let (is_match, cut_positions) = match_sub_pattern(&matches, &pattern);
        assert!(is_match);
        assert_eq!(cut_positions, expected_cut_positions);
    }

    #[test]
    fn test_simple_ftag_adapter_subpattern_matches() {
        let pattern = pattern_from_str!("Ftag[fw, *, @left(0..250), >>]");
        let matches = vec![
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
                BarcodeType::Fadapter,
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
        let seeds = find_seeds(&matches, &pattern);
        assert_eq!(seeds, vec![0]);
        // Using full pattern search this should not match
        let (is_match, _cut_positions) = match_full_pattern(&matches, &pattern);
        assert!(!is_match);

        // Using sub pattern search this should match twice
        // the firt slice from match1 --- match2, and match2-- end of read
        // First case:
        let first_slice = [
            (
                0,
                Cut {
                    group_id: 0,
                    direction: CutDirection::After,
                },
            ),
            (
                1,
                Cut {
                    group_id: 0,
                    direction: CutDirection::Before,
                },
            ),
        ]
        .to_vec();

        // Total exepcted cut positions:
        let expected_cut_positions = vec![first_slice];
        let (is_match, cut_positions) = match_sub_pattern(&matches, &pattern);
        assert!(is_match);
        assert_eq!(cut_positions, expected_cut_positions);
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
        let (is_match, _) = match_full_pattern(&matches, &pattern);
        assert!(is_match);

        matches[0].read_start_bar = 100;
        let (is_match, _) = match_full_pattern(&matches, &pattern);
        assert!(is_match);

        matches[0].read_start_bar = 250;
        let (is_match, _) = match_full_pattern(&matches, &pattern);
        assert!(is_match);

        matches[0].read_start_bar = 251;
        let (is_match, _) = match_full_pattern(&matches, &pattern);
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
        let (is_match, _) = match_full_pattern(&matches, &pattern);
        assert!(is_match);

        matches[0].read_end_bar = 450;
        let (is_match, _) = match_full_pattern(&matches, &pattern);
        assert!(is_match);

        matches[0].read_end_bar = 250;
        let (is_match, _) = match_full_pattern(&matches, &pattern);
        assert!(is_match);

        matches[0].read_end_bar = 249;
        let (is_match, _) = match_full_pattern(&matches, &pattern);
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
        let (is_match, _) = match_full_pattern(&matches, &pattern);
        assert!(!is_match);

        // Same start as previous match, does not satisfy min distance of 5
        matches[1].read_start_bar = 100;
        let (is_match, _) = match_full_pattern(&matches, &pattern);
        assert!(!is_match);

        // Should pass, just satisfies min distance of 5
        matches[1].read_start_bar = 105;
        let (is_match, _) = match_full_pattern(&matches, &pattern);
        assert!(is_match);

        // Should still pass, but maxing out
        matches[1].read_start_bar = 200;
        let (is_match, _) = match_full_pattern(&matches, &pattern);
        assert!(is_match);

        // Should not pass anymore, too far away
        matches[1].read_start_bar = 201;
        let (is_match, _) = match_full_pattern(&matches, &pattern);
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

        let (is_match, _) = match_full_pattern(&matches, &pattern);
        assert!(is_match);

        // Edit..different labels
        matches[1].label = "yyyy".to_string();

        let (is_match, _) = match_full_pattern(&matches, &pattern);
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

        let (is_match, _) = match_full_pattern(&matches, &pattern);
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

        let (is_match, _) = match_full_pattern(&matches, &pattern);
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

        let (is_match, cut_positions) = match_full_pattern(&matches, &pattern);
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

        let (is_match, cut_positions) = match_full_pattern(&matches, &pattern);
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

        let (is_match, cut_positions) = match_full_pattern(&matches, &pattern);
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
