use std::collections::HashMap;

#[derive(Debug, Eq, PartialEq, Clone)]
pub enum MatchType {
    Fbarcode, // Just specify the match type
    Rbarcode,
    Flank,
}

#[derive(Debug, PartialEq, Eq,  Clone)]
pub enum Orientation {
    Forward,
    ReverseComplement,
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Match {
    pub match_str: EncodedMatchStr,
    pub start: usize,
    pub end: usize,
    pub edits: i32,
    pub rel_dist_to_end: isize,
}



#[derive(Debug, Eq, PartialEq, Clone)]
pub struct EncodedMatchStr {
    pub match_type: MatchType,
    pub orientation: Orientation,
    pub label: Option<String>,
}

impl EncodedMatchStr {
    
    pub fn new(match_type: MatchType, orientation: Orientation, label: Option<String>) -> Self {
        Self { match_type, orientation, label }
    }

    pub fn stringify(&self) -> String {
        let label = if let Some(label) = &self.label {
            label.clone()
        } else {
            "Flank".to_string()
        };

        let ori  = match self.orientation {
            Orientation::Forward => "fw",
            Orientation::ReverseComplement => "rc",
        };

        let match_type = match self.match_type {
            MatchType::Fbarcode => "Fbar",
            MatchType::Rbarcode => "Rbar",
            MatchType::Flank => "Flank",
        };

        format!("{}#{}#{}", label, ori, match_type)
    }

    pub fn unstringify(s: &str) -> Self {
        let parts: Vec<&str> = s.split('#').collect();
        let label = parts[0].to_string();
        // Check label if label is flank we use label None, else the string
        let label = if label == "Flank" {
            None
        } else {
            Some(label)
        };

        let ori = match parts[1] {
            "fw" => Orientation::Forward,
            "rc" => Orientation::ReverseComplement,
            _ => panic!("Invalid orientation"),
        };
        let match_type = match parts[2] {
            "Fbar" => MatchType::Fbarcode,
            "Rbar" => MatchType::Rbarcode,
            "Flank" => MatchType::Flank,
            _ => panic!("Invalid match type"),
        };
        Self::new(match_type, ori, label)
    }
}



impl Match {
    
    pub fn new(
        match_str: EncodedMatchStr,
        start: usize,
        end: usize,
        edits: i32,
        rel_dist_to_end: isize, // negative for close to the right end
    ) -> Self {
        Self {
            match_str,
            start,
            end,
            edits,
            rel_dist_to_end,
        }
    }
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub enum CutDirection {
    Before, // < cut at match.start
    After,  // > cut at match.end
}

#[derive(Debug, Eq, PartialEq)]
pub struct PatternElement {
    pub match_type: MatchType,
    pub orientation: Option<Orientation>,
    pub label: Option<String>,                 // Only used for Barcode matches
    pub placeholder: Option<usize>,            // For referencing previous matches like $1, $2, etc.
    pub range: (isize, isize),                 // Minimum and maximum positions (distance)
    pub relative_to: Option<RelativePosition>, // Relative to left or right ends, or previous match
    pub cut: Option<CutDirection>,  // New field
}

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

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Cut {
    pub position: usize,
    pub direction: CutDirection,
}

impl Cut {
    pub fn new(position: usize, direction: CutDirection) -> Self {
        Self { position, direction }
    }

    pub fn to_string(&self) -> String {
        match self.direction {
            CutDirection::Before => format!("Before({})", self.position),
            CutDirection::After => format!("After({})", self.position),
        }
    }

    pub fn from_string(s: &str) -> Option<Self> {
        let s = s.trim();
        if s.starts_with("Before(") && s.ends_with(")") {
            let pos = s[7..s.len()-1].parse().ok()?;
            Some(Cut::new(pos, CutDirection::Before))
        } else if s.starts_with("After(") && s.ends_with(")") {
            let pos = s[6..s.len()-1].parse().ok()?;
            Some(Cut::new(pos, CutDirection::After))
        } else {
            None
        }
    }
}

fn check_match_type_and_label(m: &Match, pattern_element: &PatternElement) -> bool {
    // println!("Checking match type: {:?} against pattern type: {:?}", 
    //     m.match_str.match_type, pattern_element.match_type);
    
    // First check if match types are exactly equal
    if m.match_str.match_type != pattern_element.match_type {
        // println!(
        //     "Match type mismatch: {:?} != {:?}",
        //     m.match_str.match_type, pattern_element.match_type
        // );
        return false;
    }

    match pattern_element.match_type {
        MatchType::Fbarcode | MatchType::Rbarcode => {
            if let Some(ref expected_label) = pattern_element.label {
                // Labels are not the same
                if let Some(ref m_label) = m.match_str.label {
                    if expected_label != m_label {
                        // println!(
                        //     "Explicit label mismatch: {:?} != {:?}",
                        //     expected_label, m_label
                        // );
                        return false;
                    }
                // No label in match
                } else {
                    //  println!("Barcode without label in match");
                    return false;
                }
            }
            true
        }
        MatchType::Flank => {
           // println!("Checking Flank match");
            true
        }
    }
}

fn check_placeholder(
    m: &Match,
    pattern_element: &PatternElement,
    matched_labels: &mut HashMap<usize, String>,
) -> bool {
    if let Some(placeholder_key) = pattern_element.placeholder {
        if let Some(stored_label) = matched_labels.get(&placeholder_key) {
            // Check against stored label
            if let Some(ref m_label) = m.match_str.label {
                if m_label != stored_label {
                    // println!(
                    //     "Placeholder label mismatch: current {:?} != stored {:?}",
                    //     m_label, stored_label
                    // );
                    return false;
                }
            } else {
                println!("Current match has no label for placeholder comparison");
                return false;
            }
        } else {
            // Store new label
            if let Some(ref m_label) = m.match_str.label {
                matched_labels.insert(placeholder_key, m_label.clone());
            } else {
               // println!("Cannot save a placeholder because current match has no label");
                return false;
            }
        }
    }
    true
}

fn check_orientation(m: &Match, pattern_element: &PatternElement) -> bool {
    if let Some(ref orientation) = pattern_element.orientation {
        if &m.match_str.orientation != orientation {
            // println!(
            //     "Orientation mismatch: {:?} != {:?}",
            //     m.match_str.orientation, orientation
            // );
            return false;
        }
    }
    true
}

fn check_relative_position(
    m: &Match,
    pattern_element: &PatternElement,
    prev_end: Option<isize>,
    seq_len: isize,
) -> bool {
    if let Some(relative_to) = &pattern_element.relative_to {
        let m_start = m.start as isize;
        let m_end = m.end as isize;
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
                if m_start < left_bound || m_start > right_bound {
                    return false;
                }
            }
            RelativePosition::PrevLeft => {
                if let Some(prev_end_pos) = prev_end {
                    let left_bound = prev_end_pos.saturating_add(pattern_element.range.0);
                    let right_bound = prev_end_pos.saturating_add(pattern_element.range.1);
                    if m_start < left_bound || m_start > right_bound {
                        // println!(
                        //     "PrevLeft gap mismatch: m.start {} not in [{}, {}]",
                        //     m_start, left_bound, right_bound
                        // );
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
    m: &Match,
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

pub fn match_pattern(matches: &[Match], pattern: &Pattern, seq_len: usize) -> (bool, Vec<Cut>) {
    let mut prev_end: Option<isize> = None;
    let mut matched_labels: HashMap<usize, String> = HashMap::new();
    let mut current_match_idx = 0;
    let mut cut_positions: Vec<Cut> = Vec::new();

    // Early return if not enough matches
    if matches.len() < pattern.elements.len() {
        return (false, vec![]);
    }

    for (i, pattern_element) in pattern.elements.iter().enumerate() {
        if current_match_idx >= matches.len() {
            return (false, vec![]);
        }

        let m = &matches[current_match_idx];
        
        if matches_pattern_element(m, pattern_element, prev_end, &mut matched_labels, seq_len) {
            // If this element has a cut marker, record the position and direction
            if let Some(cut_dir) = &pattern_element.cut {
                let cut = match cut_dir {
                    CutDirection::Before => Cut::new(m.start, CutDirection::Before),
                    CutDirection::After => Cut::new(m.end, CutDirection::After),
                };
                cut_positions.push(cut);
            }
            
            prev_end = Some(m.end as isize);
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

        fn parse_range(range_str: &str) -> Option<(isize, isize)> {
            let parts: Vec<&str> = range_str
                .trim_matches(|p| p == '(' || p == ')')
                .split("to")
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
                "Fbarcode" => MatchType::Fbarcode,
                "Rbarcode" => MatchType::Rbarcode,
                "Flank" => MatchType::Flank,
                _ => return None,
            };

            // Initialize defaults
            let mut label = None;
            let mut placeholder = None;
            let mut orientation = None;
            let mut relative_to = None;
            let mut range = (0, 0);
            let mut cut = None;

            // Parse parameters
            let params = parts[1].trim_end_matches(']');
       //     println!("params: {:?}", params);
            for param in params.split(',').map(|s| s.trim()) {
            //    println!("\tparam: {:?}", param);
                match param {
                    // Orientation
                    "fw" => orientation = Some(Orientation::Forward),
                    "rc" => orientation = Some(Orientation::ReverseComplement),

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

                    // Cut direction
                    "<<" => cut = Some(CutDirection::Before),
                    ">>" => cut = Some(CutDirection::After),

                    // Label
                    "*" => (), // Any label
                    p => label = Some(p.trim_matches('"').to_string()),
                }
            }

            Some(PatternElement {
                match_type,
                orientation,
                label,
                placeholder,
                range,
                relative_to,
                cut,
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
            "Fbarcode[fw, *, @left(0 to 250)]__Flank[fw, @prev_left(5 to 100)]__Rbarcode[?1, fw, @right(0 to 20)]"
        );
        println!("Parsed pattern: {:?}", pattern);
        assert_eq!(pattern.elements.len(), 3);
        assert_eq!(
            pattern,
            Pattern {
                elements: vec![
                    PatternElement {
                        match_type: MatchType::Fbarcode,
                        orientation: Some(Orientation::Forward),
                        label: None,
                        placeholder: None,
                        range: (0, 250),
                        relative_to: Some(RelativePosition::Left),
                        cut: None
                    },
                    PatternElement {
                        match_type: MatchType::Flank,
                        orientation: Some(Orientation::Forward),
                        label: None,
                        placeholder: None,
                        range: (5, 100),
                        relative_to: Some(RelativePosition::PrevLeft),
                        cut: None
                    },
                    PatternElement {
                        match_type: MatchType::Rbarcode,
                        orientation: Some(Orientation::Forward),
                        label: None,
                        placeholder: Some(1),
                        range: (0, 20),
                        relative_to: Some(RelativePosition::Right),
                        cut: None
                    }
                ]
            }
        );
    }

    #[test]
    fn test_distance_to_left_end() {
       // let pattern = pattern_from_str!("Fbarcode[fw, *, @left(0-250)]");


        let pattern = pattern_from_str!("Fbarcode[fw, *, @left(0 to 250)]");

        let mut matches = vec![
            Match {
                match_str: EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("XXX".to_string())),
                start: 0,
                end: 100,
                edits: 0,
                rel_dist_to_end: 0,
            },
        ];

        matches[0].start = 0;
        let (is_match, cut_positions) = match_pattern(&matches, &pattern, 500);
        assert!(is_match);
      

        matches[0].start = 100;
        let (is_match, cut_positions) = match_pattern(&matches, &pattern, 500);
        assert!(is_match);


        matches[0].start = 250;
        let (is_match, cut_positions) = match_pattern(&matches, &pattern, 500);
        assert!(is_match);


        matches[0].start = 251;
        let (is_match, cut_positions) = match_pattern(&matches, &pattern, 500);
        assert!(!is_match);
       
    }

    #[test]
    fn test_distance_to_right_end() {
        let pattern = pattern_from_str!("Fbarcode[fw, *, @right(0 to 250)]");
       
        let mut matches = vec![ 
            Match {
                match_str: EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("XXX".to_string())),
                start: 0,
                end: 100,
                edits: 0,
                rel_dist_to_end: 0,
            },
        ];

        matches[0].start = 500;
        let (is_match, cut_positions) = match_pattern(&matches, &pattern, 500);
        assert!(is_match);

        matches[0].start = 450;
        let (is_match, cut_positions) = match_pattern(&matches, &pattern, 500);
        assert!(is_match);

        matches[0].start = 250;
        let (is_match, cut_positions) = match_pattern(&matches, &pattern, 500);
        assert!(is_match);

        matches[0].start = 249;
        let (is_match, cut_positions) = match_pattern(&matches, &pattern, 500);
        assert!(!is_match);

    }

    // #[test]
    // fn test_distance_to_prev_left() {
    //     let pattern = pattern_from_str!("Fbarcode[fw, *, @left(0 to 250)]__Flank[fw, @prev_left(5 to 100)]");
    //     println!("Pattern: {:?}", pattern);

    //     let mut matches = vec![
    //         Match {
    //             match_str: EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("XXX".to_string())),
    //             start: 0,
    //             end: 100,
    //             edits: 0,
    //             rel_dist_to_end: 0,
    //         },
    //         Match {
    //             match_str: EncodedMatchStr::new(MatchType::Flank, Orientation::Forward, None),
    //             start: 100,
    //             end: 200,
    //             edits: 0,
    //             rel_dist_to_end: 0,
    //         },
    //     ];
        
        

    //     // Should not pass, 
    //     matches[1].start = 50;
    //     let is_match = match_pattern(&matches, &pattern, 500);
    //     assert!(!is_match);

    //     // Same start as previous match, does not satisfy min distance of 5
    //     matches[1].start = 100;
    //     let is_match = match_pattern(&matches, &pattern, 500);
    //     assert!(!is_match);

    //     // Should pass, just satisfies min distance of 5
    //     matches[1].start = 105;
    //     let is_match = match_pattern(&matches, &pattern, 500);
    //     assert!(is_match);

    //     // Should still pass, but maxing out
    //     matches[1].start = 200;
    //     let is_match = match_pattern(&matches, &pattern, 500);
    //     assert!(is_match);

    //     // Should not pass anymore, too far away
    //     matches[1].start = 201;
    //     let is_match = match_pattern(&matches, &pattern, 500);
    //     assert!(!is_match);
        
    // }


    #[test]
    fn test_placeholder() {
        // Means any label for front, but the last label should be the same as the first label
        let pattern = pattern_from_str!("Fbarcode[fw, ?1, @left(0 to 250)]__Rbarcode[fw, ?1, @right(0 to 250)]");
        println!("Pattern: {:?}", pattern);

        let matches_samme_label = vec![
            Match {
                match_str: EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("XXX".to_string())),
                start: 0,
                end: 100,
                edits: 0,
                rel_dist_to_end: 0,
            },
            Match {
                match_str: EncodedMatchStr::new(MatchType::Rbarcode, Orientation::Forward, Some("XXX".to_string())),
                start: 450,
                end: 500,
                edits: 0,
                rel_dist_to_end: 0,
            },
        ];

        let matches_different_label = vec![
            Match {
                match_str: EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("XXX".to_string())),
                start: 0,
                end: 100,
                edits: 0,
                rel_dist_to_end: 0,
            },
            Match {
                match_str: EncodedMatchStr::new(MatchType::Rbarcode, Orientation::Forward, Some("YYY".to_string())),
                start: 450,
                end: 500,   
                edits: 0,
                rel_dist_to_end: 0,
            },
        ];

        let (is_match, cut_positions) = match_pattern(&matches_samme_label, &pattern, 500);
        assert!(is_match);


        let (is_match, cut_positions) = match_pattern(&matches_different_label, &pattern, 500);
        assert!(!is_match);
        
    }


    #[test]
    fn test_stringify_unstringify() {
        let m1 = EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("XXX".to_string()));
        let m1_str = m1.stringify();
        assert_eq!(m1_str, "XXX#fw#Fbar");
        let m1_unstr = EncodedMatchStr::unstringify(&m1_str);
        assert_eq!(m1, m1_unstr);
    }

    #[test]
    fn test_match_type_validation() {
        let pattern = pattern_from_str!(
            "Fbarcode[fw, ?1, @left(0 to 250)]__Fbarcode[fw, ?1, @prev_left(-100 to 200)]"
        );
        
        let matches = vec![
            Match {
                match_str: EncodedMatchStr::new(
                    MatchType::Fbarcode,
                    Orientation::Forward,
                    Some("barcode47".to_string())
                ),
                start: 0,
                end: 100,
                edits: 24,
                rel_dist_to_end: 1,
            },
            Match {
                match_str: EncodedMatchStr::new(
                    MatchType::Flank,  // This should cause the match to fail
                    Orientation::Forward,
                    None
                ),
                start: 47,
                end: 161,
                edits: 24,
                rel_dist_to_end: 47,
            },
        ];

        let (is_match, cut_positions) = match_pattern(&matches, &pattern, 500);
        assert!(!is_match, "Pattern should not match because second element is Flank when Fbarcode is required");
    }

    #[test]
    fn test_pattern_with_cuts() {
        let pattern = pattern_from_str!(
            "Fbarcode[fw, *, >>, @left(0 to 250)]__Flank[fw, <<, @prev_left(5 to 100)]"
        );
        
        let matches = vec![
            Match {
                match_str: EncodedMatchStr::new(MatchType::Fbarcode, Orientation::Forward, Some("XXX".to_string())),
                start: 10,
                end: 30,
                edits: 0,
                rel_dist_to_end: 0,
            },
            Match {
                match_str: EncodedMatchStr::new(MatchType::Flank, Orientation::Forward, None),
                start: 40,
                end: 60,
                edits: 0,
                rel_dist_to_end: 0,
            },
        ];

        let (is_match, cut_positions) = match_pattern(&matches, &pattern, 500);
        assert!(is_match);
        assert_eq!(cut_positions, vec![Cut::new(30, CutDirection::After), Cut::new(40, CutDirection::Before)]); 
    }

    #[test]
    fn test_cut_from_string() {
        let cut_str = "After(118)";
        let cut = Cut::from_string(cut_str).unwrap();
        assert_eq!(cut, Cut::new(118, CutDirection::After));

        let cut_str = "Before(50)";
        let cut = Cut::from_string(cut_str).unwrap();
        assert_eq!(cut, Cut::new(50, CutDirection::Before));

        // Test invalid input
        assert!(Cut::from_string("Invalid").is_none());
        assert!(Cut::from_string("After(abc)").is_none());
    }

}


