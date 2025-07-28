pub mod annotate;
pub mod filter;
pub mod inspect;
pub mod preset;
pub mod trim;

/// Extra padding around barcode to ensure correct "alignment"
/// when comparing to masked region
pub(crate) const WIGGLE_ROOM: usize = 5;

// For presets we compile the strings, this costs quite some space thouhg in the binary
// better ideas?
const RAPID_BARS_CONTENT: &str = include_str!("../examples/rapid_bars.fasta");
