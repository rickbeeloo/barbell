pub mod annotate;
pub mod filter;
pub mod inspect;
pub mod preset;
pub mod trim;

// fixme: not compile it in, maybe auto download?
const RAPID_BARS_CONTENT: &str = include_str!("../examples/rapid_bars.fasta");
