pub mod filter;
pub mod inspect;
pub mod progress;
pub mod search;
pub mod trim;

/// Extra padding around barcode to ensure correct "alignment"
/// when comparing to masked region
pub(crate) const WIGGLE_ROOM: usize = 5;
