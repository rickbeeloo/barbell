pub mod annotate;
pub mod filter;
pub mod inspect;
pub mod preset;
pub mod trim;
pub mod tune;

/// Extra padding around barcode to ensure correct "alignment"
/// we on purpose give more room then needed so we can discard matches
/// that "cheat" their score by using the overhang cost
pub(crate) const WIGGLE_ROOM: usize = 50;
// pub(crate) const PADDING: usize = 15; // Should be smaller than wiggle room ideally

// For presets we compile the strings, this costs quite some space thouhg in the binary
// better ideas?
const RAPID_BARS_CONTENT: &str = include_str!("../examples/rapid_bars.fasta");
// const RAPID_BARS_CONTENT: &str = include_str!("../examples/debug.fasta");
