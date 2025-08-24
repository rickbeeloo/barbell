use crate::filter::pattern::*;
use crate::{
    pattern_from_str,
    trim::trim::{LabelConfig, LabelSide},
};
use std::collections::HashMap;
use std::sync::LazyLock;

//From https://github.com/nanoporetech/dorado/blob/e72f14925cd435fff823ebf244ce2195b135a863/dorado/utils/barcode_kits.cpp
const RAB_1st_FRONT: &str = "CCGTGAC";
const RAB_1st_REAR: &str = "AGAGTTTGATCATGGCTCAG";
const RAB_2nd_FRONT: &str = "CCGTGAC";
const RAB_2nd_REAR: &str = "CGGTTACCTTGTTACGACTT";

const RBK_FRONT: &str = "TATTGCT";
const RBK_REAR: &str = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";

const RBK4_FRONT: &str = "GCTTGGGTGTTTAACC";
const RBK4_REAR: &str = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";

const RBK4_kit14_FRONT: &str = "GCTTGGGTGTTTAACC";
const RBK4_kit14_REAR: &str = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";

const RLB_FRONT: &str = "CCGTGAC";
const RLB_REAR: &str = "CGTTTTTCGTGCGCCGCTTC";

const BC_1st_FRONT: &str = "GGTGCTG";
const BC_1st_REAR: &str = "TTAACCTTTCTGTTGGTGCTGATATTGC";
const BC_2nd_FRONT: &str = "GGTGCTG";
const BC_2nd_REAR: &str = "TTAACCTACTTGCCTGTCGCTCTATCTTC";

const NB_1st_FRONT: &str = "ATTGCTAAGGTTAA";
const NB_1st_REAR: &str = "CAGCACCT";
const NB_2nd_FRONT: &str = "ATTGCTAAGGTTAA";
const NB_2nd_REAR: &str = "CAGCACCT";

const LWB_1st_FRONT: &str = "CCGTGAC";
const LWB_1st_REAR: &str = "ACTTGCCTGTCGCTCTATCTTC";
const LWB_2nd_FRONT: &str = "CCGTGAC";
const LWB_2nd_REAR: &str = "TTTCTGTTGGTGCTGATATTGC";

const MAB_FRONT: &str = "GCTTGGGTGTTTAACC";
const MAB_REAR: &str = "CCATATCCGTGTCGCCCTT";

// Storing all options like dorado takes up too much of the binary so
// we just create views as subsets of the full list
const ALL_BARS: [&str; 96] = [
    "BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11", "BC12",
    "BC13", "BC14", "BC15", "BC16", "BC17", "BC18", "BC19", "BC20", "BC21", "BC22", "BC23", "BC24",
    "BC25", "BC26", "BC27", "BC28", "BC29", "BC30", "BC31", "BC32", "BC33", "BC34", "BC35", "BC36",
    "BC37", "BC38", "BC39", "BC40", "BC41", "BC42", "BC43", "BC44", "BC45", "BC46", "BC47", "BC48",
    "BC49", "BC50", "BC51", "BC52", "BC53", "BC54", "BC55", "BC56", "BC57", "BC58", "BC59", "BC60",
    "BC61", "BC62", "BC63", "BC64", "BC65", "BC66", "BC67", "BC68", "BC69", "BC70", "BC71", "BC72",
    "BC73", "BC74", "BC75", "BC76", "BC77", "BC78", "BC79", "BC80", "BC81", "BC82", "BC83", "BC84",
    "BC85", "BC86", "BC87", "BC88", "BC89", "BC90", "BC91", "BC92", "BC93", "BC94", "BC95", "BC96",
];

// RBK_1_96 is the same as BC_1_96 except for 26, 39, 40, 48, 54 and 60.

#[derive(Copy, Clone)]
pub struct LabelRange {
    pub from: &'static str,
    pub to: &'static str,
}

impl LabelRange {
    const fn new(from: &'static str, to: &'static str) -> Self {
        Self { from, to }
    }
}

pub struct KitConfig {
    pub name: &'static str,
    pub left_prefix: Option<&'static str>,
    pub left_suffix: Option<&'static str>,
    pub right_prefix: Option<&'static str>,
    pub right_suffix: Option<&'static str>,
    pub barcodes: LabelRange,
    pub barcodes2: Option<LabelRange>,
    pub label_config: LabelConfig,
    pub safe_patterns: fn() -> &'static [Pattern],
    pub maximize_patterns: fn() -> &'static [Pattern],
}

impl KitConfig {
    const fn new(
        name: &'static str,
        left_prefix: Option<&'static str>,
        left_suffix: Option<&'static str>,
        right_prefix: Option<&'static str>,
        right_suffix: Option<&'static str>,
        barcodes: LabelRange,
        barcodes2: Option<LabelRange>,
        label_config: LabelConfig,
        safe_patterns: fn() -> &'static [Pattern],
        maximize_patterns: fn() -> &'static [Pattern],
    ) -> Self {
        Self {
            name,
            left_prefix,
            left_suffix,
            right_prefix,
            right_suffix,
            barcodes,
            barcodes2,
            label_config,
            safe_patterns,
            maximize_patterns,
        }
    }

    fn barcodes(&self) -> Vec<String> {
        get_barcodes(self.barcodes.from, self.barcodes.to)
    }

    fn barcodes2(&self) -> Vec<String> {
        if let Some(r) = self.barcodes2 {
            get_barcodes(r.from, r.to)
        } else {
            Vec::new()
        }
    }
}

/*
 Labe configurations for the different kits

*/
const SINGLE_LABEL_CONFIG: LabelConfig = LabelConfig {
    include_label: true,
    include_orientation: false,
    include_flank: false,
    sort_labels: false,
    only_side: Some(LabelSide::Left),
};

const DOUBLE_LABEL_CONFIG_KEEP_SINGLE: LabelConfig = LabelConfig {
    include_label: true,
    include_orientation: false,
    include_flank: false,
    sort_labels: false,
    only_side: Some(LabelSide::Left),
};

// See when this makes sense, not sure if any protocol uses different
// barcodes on either side
#[allow(dead_code)]
const DOUBLE_LABEL_CONFIG_KEEP_DOUBLE: LabelConfig = LabelConfig {
    include_label: true,
    include_orientation: false,
    include_flank: false,
    sort_labels: false,
    only_side: None,
};

/*
    Patterns to use for each kit, the "safe" and "maximize" patterns
*/

// Lazy-initialized pattern sets built at runtime (avoid const-eval of parsing)
static SINGLE_LABEL_PATTERNS_SAFE: LazyLock<Vec<Pattern>> = LazyLock::new(|| {
    vec![
        // Single barcode on the left
        pattern_from_str!("Ftag[fw, *, @left(0..250), >>]"),
        // Double barcode on the left but with the same barcode label (within sample ligation)
        pattern_from_str!("Ftag[fw, ?1, @left(0..250)]__Ftag[fw, ?1, @prev_left(0..250), >>]"),
    ]
});

static SINGLE_LABEL_PATTERNS_MAXIMIZE: LazyLock<Vec<Pattern>> = LazyLock::new(|| {
    vec![
        // From safe
        pattern_from_str!("Ftag[fw, *, @left(0..250), >>]"),
        pattern_from_str!("Ftag[fw, ?1, @left(0..250)]__Ftag[fw, ?1, @prev_left(0..250), >>]"),
        // Slightly more risky (ignores that both labels should be identical and just uses left one)
        pattern_from_str!("Ftag[fw, *, @left(0..250)]__Ftag[fw, *, @prev_left(0..250), >>]"),
        // Would be odd to have barcode on the right, but can still assume left is fine and extract region within
        pattern_from_str!("Ftag[fw, *, @left(0..250), >>]__Ftag[<<, fw, *, @right(0..250)]"),
        // Same as above + double left
        pattern_from_str!(
            "Ftag[fw, *, @left(0..250)]__Ftag[fw, *, @prev_left(0..250), >>]__Ftag[<<, fw, *, @right(0..250)]"
        ),
    ]
});

static DOUBLE_LABEL_PATTERNS_SAFE: LazyLock<Vec<Pattern>> = LazyLock::new(|| {
    vec![
        // Single barcode on the left
        pattern_from_str!("Ftag[fw, *, @left(0..250), >>]"),
        // Single barcode on the right
        pattern_from_str!("Ftag[<<, rc, *, @right(0..250)]"),
        // Double barcode, left and right, as expected, with both the same barcode label
        pattern_from_str!("Ftag[fw, ?1, @left(0..250), >>]__Ftag[<<, rc, ?1, @right(0..250)]"),
    ]
});

static DOUBLE_LABEL_PATTERNS_MAXIMIZE: LazyLock<Vec<Pattern>> = LazyLock::new(|| {
    vec![
        // From safe
        pattern_from_str!("Ftag[fw, *, @left(0..250), >>]"),
        pattern_from_str!("Ftag[<<, rc, *, @right(0..250)]"),
        pattern_from_str!("Ftag[fw, ?1, @left(0..250), >>]__Ftag[<<, rc, ?1, @right(0..250)]"),
        // Same as ideal, but extra barcode on the left, we ensure two "inner" barcodes have the same label
        pattern_from_str!(
            "Ftag[fw, *, @left(0..250)]__Ftag[fw, ?1, @prev_left(0..250), >>]__Ftag[<<, rc, ?1, @right(0..250)]"
        ),
        // Barcode on the left, and just flank on the right, we can't with certainty prove that the barcodes were different, so we assume it's ok
        pattern_from_str!("Ftag[fw, *, @left(0..250), >>]__Fflank[<<, rc, *, @right(0..250)]"),
        // Same as above, but flipped
        pattern_from_str!("Fflank[fw, *, @left(0..250), >>]__Ftag[<<, rc, *, @right(0..250)]"),
        // Two barcodes on the left
        pattern_from_str!("Ftag[fw, *, @left(0..250)]__Ftag[fw, *, @prev_left(0..250), >>]"),
        // Weird chimeric pattern, where we have a double Ftag on the right...
        pattern_from_str!(
            "Ftag[fw, ?1, @left(0..250), >>]__Ftag[<<, fw, ?1, @right(0..250)]__Ftag[rc, *, @right(0..250)]"
        ),
        // Triple barcode on the left, again we assure inner barcodes are the same
        pattern_from_str!(
            "Ftag[fw, *, @left(0..250)]__Ftag[rc, *, @prev_left(0..250)]__Ftag[fw, ?1, @prev_left(0..250), >>]__Ftag[<<, rc, ?1, @right(0..250)]"
        ),
    ]
});

fn single_label_patterns_safe() -> &'static [Pattern] {
    &SINGLE_LABEL_PATTERNS_SAFE
}
fn single_label_patterns_maximize() -> &'static [Pattern] {
    &SINGLE_LABEL_PATTERNS_MAXIMIZE
}
fn double_label_patterns_safe() -> &'static [Pattern] {
    &DOUBLE_LABEL_PATTERNS_SAFE
}
fn double_label_patterns_maximize() -> &'static [Pattern] {
    &DOUBLE_LABEL_PATTERNS_MAXIMIZE
}

const KIT_16S: KitConfig = KitConfig::new(
    "16S",
    Some(RAB_1st_FRONT),
    Some(RAB_1st_REAR),
    Some(RAB_2nd_FRONT),
    Some(RAB_2nd_REAR),
    LabelRange::new("BC01", "BC24"),
    Some(LabelRange::new("BC01", "BC24")),
    DOUBLE_LABEL_CONFIG_KEEP_SINGLE,
    double_label_patterns_safe,
    double_label_patterns_maximize,
);

const KIT_LWB: KitConfig = KitConfig::new(
    "LWB",
    Some(LWB_1st_FRONT),
    Some(LWB_1st_REAR),
    Some(LWB_2nd_FRONT),
    Some(LWB_2nd_REAR),
    LabelRange::new("BC01", "BC12"),
    Some(LabelRange::new("BC01", "BC12")),
    DOUBLE_LABEL_CONFIG_KEEP_SINGLE,
    double_label_patterns_safe,
    double_label_patterns_maximize,
);

const KIT_LWB24: KitConfig = KitConfig::new(
    "LWB24",
    Some(LWB_1st_FRONT),
    Some(LWB_1st_REAR),
    Some(LWB_2nd_FRONT),
    Some(LWB_2nd_REAR),
    LabelRange::new("BC01", "BC24"),
    Some(LabelRange::new("BC01", "BC24")),
    DOUBLE_LABEL_CONFIG_KEEP_SINGLE,
    double_label_patterns_safe,
    double_label_patterns_maximize,
);

const KIT_NB12: KitConfig = KitConfig::new(
    "NB12",
    Some(NB_1st_FRONT),
    Some(NB_1st_REAR),
    None,
    None,
    LabelRange::new("NB01", "NB12"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

const KIT_NB24: KitConfig = KitConfig::new(
    "NB24",
    Some(NB_1st_FRONT),
    Some(NB_1st_REAR),
    None,
    None,
    LabelRange::new("NB01", "NB24"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

const KIT_NB96: KitConfig = KitConfig::new(
    "NB96",
    Some(NB_1st_FRONT),
    Some(NB_1st_REAR),
    None,
    None,
    LabelRange::new("NB01", "NB96"),
    None,
    SINGLE_LABEL_CONFIG,
    double_label_patterns_safe,
    double_label_patterns_maximize,
);

const KIT_RAB: KitConfig = KitConfig::new(
    "RAB",
    Some(RAB_1st_FRONT),
    Some(RAB_1st_REAR),
    Some(RAB_2nd_FRONT),
    Some(RAB_2nd_REAR),
    LabelRange::new("BC01", "BC12"),
    Some(LabelRange::new("BC01", "BC12")),
    DOUBLE_LABEL_CONFIG_KEEP_SINGLE,
    double_label_patterns_safe,
    double_label_patterns_maximize,
);

const KIT_RBK96: KitConfig = KitConfig::new(
    "RBK96",
    Some(RBK4_FRONT),
    Some(RBK4_REAR),
    None,
    None,
    LabelRange::new("RBK01", "RBK96"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

const KIT_RBK4: KitConfig = KitConfig::new(
    "RBK4",
    Some(RBK4_FRONT),
    Some(RBK4_REAR),
    None,
    None,
    LabelRange::new("BC01", "BC12"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

const KIT_RLB: KitConfig = KitConfig::new(
    "RLB",
    Some(RLB_FRONT),
    Some(RLB_REAR),
    None,
    None,
    LabelRange::new("BC01", "BC12A"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

const KIT_NB13_24: KitConfig = KitConfig::new(
    "NB13-24",
    Some(NB_1st_FRONT),
    Some(NB_1st_REAR),
    None,
    None,
    LabelRange::new("NB13", "NB24"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

const KIT_PCR12: KitConfig = KitConfig::new(
    "PCR12",
    Some(BC_1st_FRONT),
    Some(BC_1st_REAR),
    Some(BC_2nd_FRONT),
    Some(BC_2nd_REAR),
    LabelRange::new("BC01", "BC12"),
    Some(LabelRange::new("BC01", "BC12")),
    DOUBLE_LABEL_CONFIG_KEEP_SINGLE,
    double_label_patterns_safe,
    double_label_patterns_maximize,
);

// Unique kits
const KIT_PCR96: KitConfig = KitConfig::new(
    "PCR96",
    Some(BC_1st_FRONT),
    Some(BC_1st_REAR),
    Some(BC_2nd_FRONT),
    Some(BC_2nd_REAR),
    LabelRange::new("BC01", "BC96"),
    Some(LabelRange::new("BC01", "BC96")),
    DOUBLE_LABEL_CONFIG_KEEP_SINGLE,
    double_label_patterns_safe,
    double_label_patterns_maximize,
);

const KIT_RBK12: KitConfig = KitConfig::new(
    "RBK",
    Some(RBK_FRONT),
    Some(RBK_REAR),
    None,
    None,
    LabelRange::new("BC01", "BC12"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

const KIT_RBK24: KitConfig = KitConfig::new(
    "RBK24",
    Some(RBK4_FRONT),
    Some(RBK4_REAR),
    None,
    None,
    LabelRange::new("BC01", "BC24"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

const KIT_RBK96_kit14: KitConfig = KitConfig::new(
    "RBK096_kit14",
    Some(RBK4_kit14_FRONT),
    Some(RBK4_kit14_REAR),
    None,
    None,
    LabelRange::new("BC01", "BC96"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

const KIT_RBK24_kit14: KitConfig = KitConfig::new(
    "RBK24_kit14",
    Some(RBK4_kit14_FRONT),
    Some(RBK4_kit14_REAR),
    None,
    None,
    LabelRange::new("BC01", "BC24"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

const KIT_RPB24_kit14: KitConfig = KitConfig::new(
    "RPB24-Kit14",
    Some(RLB_FRONT),
    Some(RLB_REAR),
    None,
    None,
    LabelRange::new("BC01", "BC24"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

const KIT_VMK: KitConfig = KitConfig::new(
    "VMK",
    Some(RBK_FRONT),
    Some(RBK_REAR),
    None,
    None,
    LabelRange::new("BC01", "BC04"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

const KIT_VMK4: KitConfig = KitConfig::new(
    "VMK4",
    Some(RBK4_FRONT),
    Some(RBK4_REAR),
    None,
    None,
    LabelRange::new("BC01", "BC10"),
    None,
    SINGLE_LABEL_CONFIG,
    single_label_patterns_safe,
    single_label_patterns_maximize,
);

pub fn get_kit_info(kit: &str) -> KitConfig {
    match kit {
        // 16S
        "SQK-16S024" => KIT_16S,
        "SQK-16S114-24" => KIT_16S,
        // LWB
        "SQK-LWB001" => KIT_LWB,
        "SQK-PBK004" => KIT_LWB,
        "SQK-PCB109" => KIT_LWB,
        "SQK-PCB110" => KIT_LWB,
        // LWB24
        "SQK-PCB111-24" => KIT_LWB24,
        "SQK-PCB114-24" => KIT_LWB24,
        // NB12
        "EXP-NBD103" => KIT_NB12,
        "EXP-NBD104" => KIT_NB12,
        // NB13-24
        "EXP-NBD114" => KIT_NB13_24,
        "SQK-NBD111-24" => KIT_NB24,
        "SQK-NBD114-24" => KIT_NB24,
        "EXP-NBD114-24" => KIT_NB24,
        // NB96
        "EXP-NBD196" => KIT_NB96,
        "SQK-MLK111-96-XL" => KIT_NB96,
        "SQK-NBD111-96" => KIT_NB96,
        "SQK-NBD114-96" => KIT_NB96,
        "SQK-MLK114-96-XL" => KIT_NB96,
        // PCR12
        "EXP-PBC001" => KIT_PCR12,
        // PCR96
        "EXP-PBC096" => KIT_PCR96,
        // RAB
        "SQK-RAB204" => KIT_RAB,
        "SQK-RAB201" => KIT_RAB,
        // RBK
        "SQK-RBK001" => KIT_RBK12,
        // RBK096
        "SQK-RBK110-96" => KIT_RBK96,
        "SQK-RBK111-96" => KIT_RBK96,
        // RBK096_kit14
        "SQK-RBK114-96" => KIT_RBK96_kit14,
        "SQK-RBK111-24" => KIT_RBK24,
        // RBK24_kit14
        "SQK-RBK114-24" => KIT_RBK24_kit14,
        //  RBK4
        "SQK-RBK004" => KIT_RBK4,
        "VSK-PTC001" => KIT_RBK4,
        "VSK-VPS001" => KIT_RBK4,
        // RLB
        "SQK-RPB004" => KIT_RLB,
        "SQK-RLB001" => KIT_RLB,
        //
        "SQK-RPB114-24" => KIT_RPB24_kit14,
        // VMK
        "VSK-VMK001" => KIT_VMK,
        // VMK4
        "VSK-VMK004" => KIT_VMK4,
        // if name contains "." try to replace it and run again
        _ => {
            if kit.contains(".") {
                println!("Your kit name used '.' instead of '-' replaced it and trying again");
                let new_kit = kit.replace(".", "-");
                get_kit_info(&new_kit)
            } else {
                panic!("Unknown or unsupported kit: {kit}, please raise an issue")
            }
        }
    }
}

fn parse_label_simple(label: &str) -> (String, usize, bool) {
    // Formats like BC1, BC01, NB12A, RBK26, case-insensitive
    let upper = label.to_ascii_uppercase();
    let mut chars = upper.chars().peekable();
    let mut prefix = String::new();
    while let Some(&c) = chars.peek() {
        if c.is_ascii_alphabetic() {
            prefix.push(c);
            chars.next();
        } else {
            break;
        }
    }
    let mut num_str = String::new();
    while let Some(&c) = chars.peek() {
        if c.is_ascii_digit() {
            num_str.push(c);
            chars.next();
        } else {
            break;
        }
    }
    let a_flag = matches!(chars.peek(), Some('A'));
    if a_flag {
        chars.next();
    }

    let number: usize = num_str.parse().expect("Invalid numeric part in label");
    (prefix, number, a_flag)
}

pub fn get_barcodes(from_label: &str, to_label: &str) -> Vec<String> {
    let (pf_from, from_num, from_a) = parse_label_simple(from_label);
    let (pf_to, to_num, to_a) = parse_label_simple(to_label);

    assert!(
        pf_from == pf_to,
        "Mismatched label prefixes: {pf_from} vs {pf_to}"
    );

    let (start, end) = if from_num <= to_num {
        (from_num, to_num)
    } else {
        (to_num, from_num)
    };

    // Base slice uses BC labels 1..=96
    let slice_from = start - 1;
    let slice_to = end;
    let mut slice: Vec<String> = ALL_BARS[slice_from..slice_to]
        .iter()
        .map(|&s| s.to_string())
        .collect();

    // 12A handling: if either boundary has 'A' and the range includes 12
    let use_12a = (from_a || to_a) && (start <= 12 && 12 <= end);
    if use_12a {
        slice.iter_mut().for_each(|item| {
            if item == "BC12" {
                *item = "BC12A".to_string();
            }
        });
    }

    // NB kit: replace all BC with NB
    if pf_from == "NB" {
        slice.iter_mut().for_each(|item| {
            if item.starts_with("BC") {
                *item = item.replacen("BC", "NB", 1);
            }
        });
    }

    // RBK kit: relabel specific indices to RBK
    if pf_from == "RBK" {
        let special = [26usize, 39, 40, 48, 54, 60];
        slice.iter_mut().for_each(|item| {
            let label = item.as_str();
            if label.starts_with("BC") && label.len() >= 4 {
                let num_str = &label[2..4];
                if let Ok(n) = num_str.parse::<usize>()
                    && special.contains(&n)
                {
                    *item = item.replacen("BC", "RBK", 1);
                }
            }
        });
    }

    slice
}

// Compact, const-friendly storage of barcode sequences
const BC_SEQS: [&str; 96] = [
    "AAGAAAGTTGTCGGTGTCTTTGTG",
    "TCGATTCCGTTTGTAGTCGTCTGT",
    "GAGTCTTGTGTCCCAGTTACCAGG",
    "TTCGGATTCTATCGTGTTTCCCTA",
    "CTTGTCCAGGGTTTGTGTAACCTT",
    "TTCTCGCAAAGGCAGAAAGTAGTC",
    "GTGTTACCGTGGGAATGAATCCTT",
    "TTCAGGGAACAAACCAAGTTACGT",
    "AACTAGGCACAGCGAGTCTTGGTT",
    "AAGCGTTGAAACCTTTGTCCTCTC",
    "GTTTCATCTATCGGAGGGAATGGA",
    "CAGGTAGAAAGAAGCAGAATCGGA",
    "AGAACGACTTCCATACTCGTGTGA",
    "AACGAGTCTCTTGGGACCCATAGA",
    "AGGTCTACCTCGCTAACACCACTG",
    "CGTCAACTGACAGTGGTTCGTACT",
    "ACCCTCCAGGAAAGTACCTCTGAT",
    "CCAAACCCAACAACCTAGATAGGC",
    "GTTCCTCGTGCAGTGTCAAGAGAT",
    "TTGCGTCCTGTTACGAGAACTCAT",
    "GAGCCTCTCATTGTCCGTTCTCTA",
    "ACCACTGCCATGTATCAAAGTACG",
    "CTTACTACCCAGTGAACCTCCTCG",
    "GCATAGTTCTGCATGATGGGTTAG",
    "GTAAGTTGGGTATGCAACGCAATG",
    "CATACAGCGACTACGCATTCTCAT",
    "CGACGGTTAGATTCACCTCTTACA",
    "TGAAACCTAAGAAGGCACCGTATC",
    "CTAGACACCTTGGGTTGACAGACC",
    "TCAGTGAGGATCTACTTCGACCCA",
    "TGCGTACAGCAATCAGTTACATTG",
    "CCAGTAGAAGTCCGACAACGTCAT",
    "CAGACTTGGTACGGTTGGGTAACT",
    "GGACGAAGAACTCAAGTCAAAGGC",
    "CTACTTACGAAGCTGAGGGACTGC",
    "ATGTCCCAGTTAGAGGAGGAAACA",
    "GCTTGCGATTGATGCTTAGTATCA",
    "ACCACAGGAGGACGATACAGAGAA",
    "CCACAGTGTCAACTAGAGCCTCTC",
    "TAGTTTGGATGACCAAGGATAGCC",
    "GGAGTTCGTCCAGAGAAGTACACG",
    "CTACGTGTAAGGCATACCTGCCAG",
    "CTTTCGTTGTTGACTCGACGGTAG",
    "AGTAGAAAGGGTTCCTTCCCACTC",
    "GATCCAACAGAGATGCCTTCAGTG",
    "GCTGTGTTCCACTTCATTCTCCTG",
    "GTGCAACTTTCCCACAGGTAGTTC",
    "CATCTGGAACGTGGTACACCTGTA",
    "ACTGGTGCAGCTTTGAACATCTAG",
    "ATGGACTTTGGTAACTTCCTGCGT",
    "GTTGAATGAGCCTACTGGGTCCTC",
    "TGAGAGACAAGATTGTTCGTGGAC",
    "AGATTCAGACCGTCTCATGCAAAG",
    "CAAGAGCTTTGACTAAGGAGCATG",
    "TGGAAGATGAGACCCTGATCTACG",
    "TCACTACTCAACAGGTGGCATGAA",
    "GCTAGGTCAATCTCCTTCGGAAGT",
    "CAGGTTACTCCTCCGTGAGTCTGA",
    "TCAATCAAGAAGGGAAAGCAAGGT",
    "CATGTTCAACCAAGGCTTCTATGG",
    "AGAGGGTACTATGTGCCTCAGCAC",
    "CACCCACACTTACTTCAGGACGTA",
    "TTCTGAAGTTCCTGGGTCTTGAAC",
    "GACAGACACCGTTCATCGACTTTC",
    "TTCTCAGTCTTCCTCCAGACAAGG",
    "CCGATCCTTGTGGCTTCTAACTTC",
    "GTTTGTCATACTCGTGTGCTCACC",
    "GAATCTAAGCAAACACGAAGGTGG",
    "TACAGTCCGAGCCTCATGTGATCT",
    "ACCGAGATCCTACGAATGGAGTGT",
    "CCTGGGAGCATCAGGTAGTAACAG",
    "TAGCTGACTGTCTTCCATACCGAC",
    "AAGAAACAGGATGACAGAACCCTC",
    "TACAAGCATCCCAACACTTCCACT",
    "GACCATTGTGATGAACCCTGTTGT",
    "ATGCTTGTTACATCAACCCTGGAC",
    "CGACCTGTTTCTCAGGGATACAAC",
    "AACAACCGAACCTTTGAATCAGAA",
    "TCTCGGAGATAGTTCTCACTGCTG",
    "CGGATGAACATAGGATAGCGATTC",
    "CCTCATCTTGTGAAGTTGTTTCGG",
    "ACGGTATGTCGAGTTCCAGGACTA",
    "TGGCTTGATCTAGGTAAGGTCGAA",
    "GTAGTGGACCTAGAACCTGTGCCA",
    "AACGGAGGAGTTAGTTGGATGATC",
    "AGGTGATCCCAACAAGCGTAAGTA",
    "TACATGCTCCTGTTGTTAGGGAGG",
    "TCTTCTACTACCGATCCGAAGCAG",
    "ACAGCATCAATGTTTGGCTAGTTG",
    "GATGTAGAGGGTACGGTTTGAGGC",
    "GGCTCCATAGGAACTCACGCTACT",
    "TTGTGAGTGGAAAGATACAGGACC",
    "AGTTTCCATCACTTCAGACTTGGG",
    "GATTGTCCTCAAACTGCCACCTAC",
    "CCTGTCTGGAAGAAGAATGGACTT",
    "CTGAACGGTCATAGAGTCCACCAT",
];

const BP_SEQS: [&str; 24] = [
    "CAAGAAAGTTGTCGGTGTCTTTGTGAC",
    "CTCGATTCCGTTTGTAGTCGTCTGTAC",
    "CGAGTCTTGTGTCCCAGTTACCAGGAC",
    "CTTCGGATTCTATCGTGTTTCCCTAAC",
    "CCTTGTCCAGGGTTTGTGTAACCTTAC",
    "CTTCTCGCAAAGGCAGAAAGTAGTCAC",
    "CGTGTTACCGTGGGAATGAATCCTTAC",
    "CTTCAGGGAACAAACCAAGTTACGTAC",
    "CAACTAGGCACAGCGAGTCTTGGTTAC",
    "CAAGCGTTGAAACCTTTGTCCTCTCAC",
    "CGTTTCATCTATCGGAGGGAATGGAAC",
    "CCAGGTAGAAAGAAGCAGAATCGGAAC",
    "CAGAACGACTTCCATACTCGTGTGAAC",
    "CAACGAGTCTCTTGGGACCCATAGAAC",
    "CAGGTCTACCTCGCTAACACCACTGAC",
    "CCGTCAACTGACAGTGGTTCGTACTAC",
    "CACCCTCCAGGAAAGTACCTCTGATAC",
    "CCCAAACCCAACAACCTAGATAGGCAC",
    "CGTTCCTCGTGCAGTGTCAAGAGATAC",
    "CTTGCGTCCTGTTACGAGAACTCATAC",
    "CGAGCCTCTCATTGTCCGTTCTCTAAC",
    "CACCACTGCCATGTATCAAAGTACGAC",
    "CCTTACTACCCAGTGAACCTCCTCGAC",
    "CGCATAGTTCTGCATGATGGGTTAGAC",
];

const NB_SEQS: [&str; 96] = [
    "CACAAAGACACCGACAACTTTCTT",
    "ACAGACGACTACAAACGGAATCGA",
    "CCTGGTAACTGGGACACAAGACTC",
    "TAGGGAAACACGATAGAATCCGAA",
    "AAGGTTACACAAACCCTGGACAAG",
    "GACTACTTTCTGCCTTTGCGAGAA",
    "AAGGATTCATTCCCACGGTAACAC",
    "ACGTAACTTGGTTTGTTCCCTGAA",
    "AACCAAGACTCGCTGTGCCTAGTT",
    "GAGAGGACAAAGGTTTCAACGCTT",
    "TCCATTCCCTCCGATAGATGAAAC",
    "TCCGATTCTGCTTCTTTCTACCTG",
    "AGAACGACTTCCATACTCGTGTGA",
    "AACGAGTCTCTTGGGACCCATAGA",
    "AGGTCTACCTCGCTAACACCACTG",
    "CGTCAACTGACAGTGGTTCGTACT",
    "ACCCTCCAGGAAAGTACCTCTGAT",
    "CCAAACCCAACAACCTAGATAGGC",
    "GTTCCTCGTGCAGTGTCAAGAGAT",
    "TTGCGTCCTGTTACGAGAACTCAT",
    "GAGCCTCTCATTGTCCGTTCTCTA",
    "ACCACTGCCATGTATCAAAGTACG",
    "CTTACTACCCAGTGAACCTCCTCG",
    "GCATAGTTCTGCATGATGGGTTAG",
    "GTAAGTTGGGTATGCAACGCAATG",
    "CATACAGCGACTACGCATTCTCAT",
    "CGACGGTTAGATTCACCTCTTACA",
    "TGAAACCTAAGAAGGCACCGTATC",
    "CTAGACACCTTGGGTTGACAGACC",
    "TCAGTGAGGATCTACTTCGACCCA",
    "TGCGTACAGCAATCAGTTACATTG",
    "CCAGTAGAAGTCCGACAACGTCAT",
    "CAGACTTGGTACGGTTGGGTAACT",
    "GGACGAAGAACTCAAGTCAAAGGC",
    "CTACTTACGAAGCTGAGGGACTGC",
    "ATGTCCCAGTTAGAGGAGGAAACA",
    "GCTTGCGATTGATGCTTAGTATCA",
    "ACCACAGGAGGACGATACAGAGAA",
    "CCACAGTGTCAACTAGAGCCTCTC",
    "TAGTTTGGATGACCAAGGATAGCC",
    "GGAGTTCGTCCAGAGAAGTACACG",
    "CTACGTGTAAGGCATACCTGCCAG",
    "CTTTCGTTGTTGACTCGACGGTAG",
    "AGTAGAAAGGGTTCCTTCCCACTC",
    "GATCCAACAGAGATGCCTTCAGTG",
    "GCTGTGTTCCACTTCATTCTCCTG",
    "GTGCAACTTTCCCACAGGTAGTTC",
    "CATCTGGAACGTGGTACACCTGTA",
    "ACTGGTGCAGCTTTGAACATCTAG",
    "ATGGACTTTGGTAACTTCCTGCGT",
    "GTTGAATGAGCCTACTGGGTCCTC",
    "TGAGAGACAAGATTGTTCGTGGAC",
    "AGATTCAGACCGTCTCATGCAAAG",
    "CAAGAGCTTTGACTAAGGAGCATG",
    "TGGAAGATGAGACCCTGATCTACG",
    "TCACTACTCAACAGGTGGCATGAA",
    "GCTAGGTCAATCTCCTTCGGAAGT",
    "CAGGTTACTCCTCCGTGAGTCTGA",
    "TCAATCAAGAAGGGAAAGCAAGGT",
    "CATGTTCAACCAAGGCTTCTATGG",
    "AGAGGGTACTATGTGCCTCAGCAC",
    "CACCCACACTTACTTCAGGACGTA",
    "TTCTGAAGTTCCTGGGTCTTGAAC",
    "GACAGACACCGTTCATCGACTTTC",
    "TTCTCAGTCTTCCTCCAGACAAGG",
    "CCGATCCTTGTGGCTTCTAACTTC",
    "GTTTGTCATACTCGTGTGCTCACC",
    "GAATCTAAGCAAACACGAAGGTGG",
    "TACAGTCCGAGCCTCATGTGATCT",
    "ACCGAGATCCTACGAATGGAGTGT",
    "CCTGGGAGCATCAGGTAGTAACAG",
    "TAGCTGACTGTCTTCCATACCGAC",
    "AAGAAACAGGATGACAGAACCCTC",
    "TACAAGCATCCCAACACTTCCACT",
    "GACCATTGTGATGAACCCTGTTGT",
    "ATGCTTGTTACATCAACCCTGGAC",
    "CGACCTGTTTCTCAGGGATACAAC",
    "AACAACCGAACCTTTGAATCAGAA",
    "TCTCGGAGATAGTTCTCACTGCTG",
    "CGGATGAACATAGGATAGCGATTC",
    "CCTCATCTTGTGAAGTTGTTTCGG",
    "ACGGTATGTCGAGTTCCAGGACTA",
    "TGGCTTGATCTAGGTAAGGTCGAA",
    "GTAGTGGACCTAGAACCTGTGCCA",
    "AACGGAGGAGTTAGTTGGATGATC",
    "AGGTGATCCCAACAAGCGTAAGTA",
    "TACATGCTCCTGTTGTTAGGGAGG",
    "TCTTCTACTACCGATCCGAAGCAG",
    "ACAGCATCAATGTTTGGCTAGTTG",
    "GATGTAGAGGGTACGGTTTGAGGC",
    "GGCTCCATAGGAACTCACGCTACT",
    "TTGTGAGTGGAAAGATACAGGACC",
    "AGTTTCCATCACTTCAGACTTGGG",
    "GATTGTCCTCAAACTGCCACCTAC",
    "CCTGTCTGGAAGAAGAATGGACTT",
    "CTGAACGGTCATAGAGTCCACCAT",
];

// 12A special (used for RLB kit and as BC12A/NB12A when requested)
const BC12A_SEQ: &str = "GTTGAGTTACAAAGCACCGATCAG";

pub fn lookup_barcode_seq(label: &str) -> Option<&'static str> {
    let (prefix, number, is_a) = parse_label_simple(label);
    match prefix.as_str() {
        "BC" => {
            if is_a && number == 12 {
                return Some(BC12A_SEQ);
            }
            BC_SEQS.get(number.saturating_sub(1)).copied()
        }
        "NB" => {
            // NB12A handled as same as BC12A but with NB label upstream; sequence is same as BC12A orientation expected by pipeline
            if is_a && number == 12 {
                return Some(BC12A_SEQ);
            }
            NB_SEQS.get(number.saturating_sub(1)).copied()
        }
        "BP" => BP_SEQS.get(number.saturating_sub(1)).copied(),
        "RBK" => match number {
            26 => Some("ACTATGCCTTTCCGTGAAACAGTT"),
            39 => Some("TCTGCCACACACTCGTAAGTCCTT"),
            40 => Some("GTCGATACTGGACCTATCCCTTGG"),
            48 => Some("GAGTCCGTGACAACTTCTGAAAGC"),
            54 => Some("GGGTGCCAACTACATACCAAACCT"),
            60 => Some("GAACCCTACTTTGGACAGACACCT"),
            _ => BC_SEQS.get(number.saturating_sub(1)).copied(),
        },
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_barcodes_bc_1_to_12() {
        let barcodes = get_barcodes("BC01", "BC12");
        let expected = vec![
            "BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11",
            "BC12",
        ];
        assert_eq!(barcodes, expected);
    }

    #[test]
    fn test_get_barcodes_bc_1_to_12_with_12a() {
        let barcodes = get_barcodes("BC1A", "BC12A");
        let expected = vec![
            "BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11",
            "BC12A",
        ];
        assert_eq!(barcodes, expected);
    }

    #[test]
    fn test_get_barcodes_bc_1_to_13_with_12a() {
        let barcodes = get_barcodes("BC1A", "BC13A");
        let expected = vec![
            "BC01", "BC02", "BC03", "BC04", "BC05", "BC06", "BC07", "BC08", "BC09", "BC10", "BC11",
            "BC12A", "BC13",
        ];
        assert_eq!(barcodes, expected);
    }

    #[test]
    fn test_get_barcodes_nb_1_to_12() {
        let barcodes = get_barcodes("NB01", "NB12");
        let expected = vec![
            "NB01", "NB02", "NB03", "NB04", "NB05", "NB06", "NB07", "NB08", "NB09", "NB10", "NB11",
            "NB12",
        ];
        assert_eq!(barcodes, expected);
    }

    #[test]
    fn test_get_barcodes_rbk_special_relabel() {
        let barcodes = get_barcodes("RBK24", "RBK28");
        let expected = vec!["BC24", "BC25", "RBK26", "BC27", "BC28"];
        assert_eq!(barcodes, expected);
    }

    #[test]
    fn lookup_barcode_seq_bc_12a() {
        let seq = lookup_barcode_seq("BC12A");
        assert_eq!(seq, Some("GTTGAGTTACAAAGCACCGATCAG"));
    }
}

/*
struct KitInfo {
    std::string name;
   # bool double_ends;
   # bool ends_different;
   # bool rear_only_barcodes;
   # bool rna_barcodes;

    # Top pattern
    std::string top_front_flank;
    std::string top_rear_flank;

    # Bottom pattern
    std::string bottom_front_flank;
    std::string bottom_rear_flank;

    # Top barcode list
    std::vector<std::string> barcodes;

    # Bottom barcode list
    std::vector<std::string> barcodes2;

    #BarcodeKitScoringParams scoring_params;
};
*/
