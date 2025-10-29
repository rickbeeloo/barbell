use std::fmt;

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchMode {
    Adapter = 0,
    JustBars = 1,
    BarsAndFlanks = 2,
}

impl SearchMode {
    pub fn parse_str(s: &str) -> Option<Self> {
        // Try parsing as number first
        if let Ok(num) = s.parse::<u8>() {
            return match num {
                0 => Some(SearchMode::Adapter),
                1 => Some(SearchMode::JustBars),
                2 => Some(SearchMode::BarsAndFlanks),
                _ => None,
            };
        }

        // Fall back to string matching
        match s.to_lowercase().as_str() {
            "adapter" => Some(SearchMode::Adapter),
            "just-bars" | "justbars" | "bars" => Some(SearchMode::JustBars),
            "bars-and-flanks" | "barsandflanks" | "flanks" => Some(SearchMode::BarsAndFlanks),
            _ => None,
        }
    }

    pub fn help_text() -> &'static str {
        "Search strategy modes (can use numbers or names, comma-separated):
    0. adapter         - Search for adapter sequences
    1. just-bars       - Search for barcode sequences only (4-5x slower!)
    2. bars-and-flanks - Search for barcodes with flanking regions
  
  Examples: --search-strategy adapter,just-bars
            --search-strategy 0,1
            --search-strategy bars-and-flanks\n"
    }

    pub fn description(&self) -> &'static str {
        match self {
            SearchMode::Adapter => "Adapter search",
            SearchMode::JustBars => "Barcode-only search",
            SearchMode::BarsAndFlanks => "Barcodes with flanks",
        }
    }
}

impl fmt::Display for SearchMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SearchMode::Adapter => write!(f, "adapter"),
            SearchMode::JustBars => write!(f, "just-bars"),
            SearchMode::BarsAndFlanks => write!(f, "bars-and-flanks"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SearchStrategy {
    search: u64,
}

impl SearchStrategy {
    pub fn new() -> Self {
        Self { search: 0 }
    }

    pub fn enable(&mut self, mode: SearchMode) {
        self.search |= 1u64 << (mode as u8);
    }

    pub fn disable(&mut self, mode: SearchMode) {
        self.search &= !(1u64 << (mode as u8));
    }

    pub fn is_enabled(&self, mode: SearchMode) -> bool {
        (self.search & (1u64 << (mode as u8))) != 0
    }

    pub fn enable_adapter(&mut self) {
        self.enable(SearchMode::Adapter);
    }

    pub fn enable_just_bars(&mut self) {
        self.enable(SearchMode::JustBars);
    }

    pub fn enable_bars_and_flanks(&mut self) {
        self.enable(SearchMode::BarsAndFlanks);
    }

    pub fn adapter_enabled(&self) -> bool {
        self.is_enabled(SearchMode::Adapter)
    }

    pub fn bars_and_flanks_enabled(&self) -> bool {
        self.is_enabled(SearchMode::BarsAndFlanks)
    }

    pub fn just_bars_enabled(&self) -> bool {
        self.is_enabled(SearchMode::JustBars)
    }

    /// Parse user input so we can set the bits
    pub fn from_modes(modes: &[SearchMode]) -> Self {
        let mut strategy = Self::new();
        for &mode in modes {
            strategy.enable(mode);
        }
        strategy
    }

    /// Display the search strategy in a human-readable format
    pub fn display(&self) -> String {
        let modes = [
            SearchMode::Adapter,
            SearchMode::JustBars,
            SearchMode::BarsAndFlanks,
        ];

        let mut lines = vec!["Search Strategy:".to_string()];
        for mode in &modes {
            let symbol = if self.is_enabled(*mode) { "✓" } else { "✗" };
            lines.push(format!("  {} {}", symbol, mode.description()));
        }
        lines.join("\n")
    }

    /// Get a compact one-line representation
    pub fn display_compact(&self) -> String {
        let modes = [
            SearchMode::Adapter,
            SearchMode::JustBars,
            SearchMode::BarsAndFlanks,
        ];

        let enabled: Vec<_> = modes
            .iter()
            .filter(|&&mode| self.is_enabled(mode))
            .map(|mode| mode.to_string())
            .collect();

        if enabled.is_empty() {
            "none".to_string()
        } else {
            enabled.join(", ")
        }
    }
}

impl Default for SearchStrategy {
    fn default() -> Self {
        Self::new()
    }
}

// CLI parse
impl std::str::FromStr for SearchStrategy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut strategy = SearchStrategy::new();
        for part in s.split(',') {
            let mode = SearchMode::parse_str(part.trim()).ok_or_else(|| {
                format!("Invalid search mode: '{part}'. {}", SearchMode::help_text())
            })?;
            strategy.enable(mode);
        }
        Ok(strategy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_operations() {
        let mut strat = SearchStrategy::new();
        assert!(!strat.adapter_enabled());

        strat.enable_adapter();
        assert!(strat.adapter_enabled());
        assert!(!strat.just_bars_enabled());

        strat.enable_just_bars();
        assert!(strat.adapter_enabled());
        assert!(strat.just_bars_enabled());
    }

    #[test]
    fn test_from_str() {
        let strat: SearchStrategy = "adapter,just-bars".parse().unwrap();
        assert!(strat.adapter_enabled());
        assert!(strat.just_bars_enabled());
        assert!(!strat.bars_and_flanks_enabled());
    }

    #[test]
    fn test_from_numbers() {
        let strat: SearchStrategy = "0,1".parse().unwrap();
        assert!(strat.adapter_enabled());
        assert!(strat.just_bars_enabled());
        assert!(!strat.bars_and_flanks_enabled());
    }

    #[test]
    fn test_mixed_format() {
        let strat: SearchStrategy = "0,bars-and-flanks".parse().unwrap();
        assert!(strat.adapter_enabled());
        assert!(strat.bars_and_flanks_enabled());
        assert!(!strat.just_bars_enabled());
    }

    #[test]
    fn test_disable() {
        let mut strat = SearchStrategy::from_modes(&[SearchMode::Adapter, SearchMode::JustBars]);
        assert!(strat.adapter_enabled());
        assert!(strat.just_bars_enabled());

        strat.disable(SearchMode::Adapter);
        assert!(!strat.adapter_enabled());
        assert!(strat.just_bars_enabled());
    }

    #[test]
    fn test_display() {
        let mut strat = SearchStrategy::new();
        strat.enable_adapter();
        strat.enable_just_bars();

        let output = strat.display();
        assert!(output.contains("✓ Adapter search"));
        assert!(output.contains("✓ Barcode-only search"));
        assert!(output.contains("✗ Barcodes with flanks"));
    }

    #[test]
    fn test_display_compact() {
        let mut strat = SearchStrategy::new();
        strat.enable_adapter();
        strat.enable_just_bars();

        assert_eq!(strat.display_compact(), "adapter, just-bars");

        let empty = SearchStrategy::new();
        assert_eq!(empty.display_compact(), "none");
    }
}
