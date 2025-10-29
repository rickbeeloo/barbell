use std::fmt;

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterMode {
    Exact = 0,
    Flexible = 1,
    Terminal = 2,
    UniqueLabels = 3,
    Clean = 4,
}

impl FilterMode {
    pub fn from_str(s: &str) -> Option<Self> {
        // Try parsing as number first
        if let Ok(num) = s.parse::<u8>() {
            return match num {
                0 => Some(FilterMode::Exact),
                1 => Some(FilterMode::Flexible),
                2 => Some(FilterMode::Terminal),
                3 => Some(FilterMode::UniqueLabels),
                4 => Some(FilterMode::Clean),
                _ => None,
            };
        }

        // Fall back to string matching
        match s.to_lowercase().as_str() {
            "exact" => Some(FilterMode::Exact),
            "flexible" | "flex" => Some(FilterMode::Flexible),
            "terminal" | "term" => Some(FilterMode::Terminal),
            "unique-labels" | "uniquelabels" | "unique" => Some(FilterMode::UniqueLabels),
            "clean" => Some(FilterMode::Clean),
            _ => None,
        }
    }

    pub fn help_text() -> &'static str {
        "Filter strategy modes (can use numbers or names, comma-separated):
    0. exact          - Exact matching only
    1. flexible       - Use sub-pattern matching without prefilters
    2. terminal       - Discard reads with internal matches
    3. unique-labels  - Only keep reads where all labels are identical (e.g. NB01---NB01, but not NB01---NB02)
    4. clean          - Discard reads that are 75% or more contaminated
  
  Examples: --filter-strategy exact,terminal
            --filter-strategy 0,2,4
            --filter-strategy flexible,clean\n"
    }

    pub fn description(&self) -> &'static str {
        match self {
            FilterMode::Exact => "Exact matching only",
            FilterMode::Flexible => "Flexible sub-pattern matching",
            FilterMode::Terminal => "Discard reads with internal matches",
            FilterMode::UniqueLabels => "Keep reads with identical labels only",
            FilterMode::Clean => "Discard contaminated reads (≥75%)",
        }
    }
}

impl fmt::Display for FilterMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FilterMode::Exact => write!(f, "exact"),
            FilterMode::Flexible => write!(f, "flexible"),
            FilterMode::Terminal => write!(f, "terminal"),
            FilterMode::UniqueLabels => write!(f, "unique-labels"),
            FilterMode::Clean => write!(f, "clean"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FilterStrategy {
    filters: u64,
}

impl FilterStrategy {
    pub fn new() -> Self {
        Self { filters: 0 }
    }

    pub fn enable(&mut self, mode: FilterMode) {
        self.filters |= 1u64 << (mode as u8);
    }

    pub fn disable(&mut self, mode: FilterMode) {
        self.filters &= !(1u64 << (mode as u8));
    }

    pub fn is_enabled(&self, mode: FilterMode) -> bool {
        (self.filters & (1u64 << (mode as u8))) != 0
    }

    pub fn enable_exact(&mut self) {
        self.enable(FilterMode::Exact);
    }

    pub fn enable_flexible(&mut self) {
        self.enable(FilterMode::Flexible);
    }

    pub fn enable_terminal(&mut self) {
        self.enable(FilterMode::Terminal);
    }

    pub fn enable_unique_labels(&mut self) {
        self.enable(FilterMode::UniqueLabels);
    }

    pub fn enable_clean(&mut self) {
        self.enable(FilterMode::Clean);
    }

    pub fn exact_enabled(&self) -> bool {
        self.is_enabled(FilterMode::Exact)
    }

    pub fn flexible_enabled(&self) -> bool {
        self.is_enabled(FilterMode::Flexible)
    }

    pub fn terminal_enabled(&self) -> bool {
        self.is_enabled(FilterMode::Terminal)
    }

    pub fn unique_labels_enabled(&self) -> bool {
        self.is_enabled(FilterMode::UniqueLabels)
    }

    pub fn clean_enabled(&self) -> bool {
        self.is_enabled(FilterMode::Clean)
    }

    /// Parse user input to set the filter bits
    pub fn from_modes(modes: &[FilterMode]) -> Self {
        let mut strategy = Self::new();
        for &mode in modes {
            strategy.enable(mode);
        }
        strategy
    }

    /// Display the filter strategy in a human-readable format
    pub fn display(&self) -> String {
        let modes = [
            FilterMode::Exact,
            FilterMode::Flexible,
            FilterMode::Terminal,
            FilterMode::UniqueLabels,
            FilterMode::Clean,
        ];

        let mut lines = vec!["Filter Strategy:".to_string()];
        for mode in &modes {
            let symbol = if self.is_enabled(*mode) { "✓" } else { "✗" };
            lines.push(format!("  {} {}", symbol, mode.description()));
        }
        lines.join("\n")
    }

    /// Get a compact one-line representation
    pub fn display_compact(&self) -> String {
        let modes = [
            FilterMode::Exact,
            FilterMode::Flexible,
            FilterMode::Terminal,
            FilterMode::UniqueLabels,
            FilterMode::Clean,
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

impl Default for FilterStrategy {
    fn default() -> Self {
        Self::new()
    }
}

// CLI parse
impl std::str::FromStr for FilterStrategy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut strategy = FilterStrategy::new();
        for part in s.split(',') {
            let mode = FilterMode::from_str(part.trim()).ok_or_else(|| {
                format!("Invalid filter mode: '{part}'. {}", FilterMode::help_text())
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
        let mut strat = FilterStrategy::new();
        assert!(!strat.exact_enabled());

        strat.enable_exact();
        assert!(strat.exact_enabled());
        assert!(!strat.flexible_enabled());

        strat.enable_flexible();
        assert!(strat.exact_enabled());
        assert!(strat.flexible_enabled());
    }

    #[test]
    fn test_from_str() {
        let strat: FilterStrategy = "exact,terminal,clean".parse().unwrap();
        assert!(strat.exact_enabled());
        assert!(strat.terminal_enabled());
        assert!(strat.clean_enabled());
        assert!(!strat.flexible_enabled());
        assert!(!strat.unique_labels_enabled());
    }

    #[test]
    fn test_from_numbers() {
        let strat: FilterStrategy = "0,2,4".parse().unwrap();
        assert!(strat.exact_enabled());
        assert!(strat.terminal_enabled());
        assert!(strat.clean_enabled());
        assert!(!strat.flexible_enabled());
        assert!(!strat.unique_labels_enabled());
    }

    #[test]
    fn test_mixed_format() {
        let strat: FilterStrategy = "0,flexible,4".parse().unwrap();
        assert!(strat.exact_enabled());
        assert!(strat.flexible_enabled());
        assert!(strat.clean_enabled());
    }

    #[test]
    fn test_disable() {
        let mut strat = FilterStrategy::from_modes(&[FilterMode::Exact, FilterMode::Clean]);
        assert!(strat.exact_enabled());
        assert!(strat.clean_enabled());

        strat.disable(FilterMode::Exact);
        assert!(!strat.exact_enabled());
        assert!(strat.clean_enabled());
    }

    #[test]
    fn test_display() {
        let mut strat = FilterStrategy::new();
        strat.enable_exact();
        strat.enable_terminal();
        strat.enable_clean();

        let output = strat.display();
        assert!(output.contains("✓ Exact matching only"));
        assert!(output.contains("✗ Flexible sub-pattern matching"));
        assert!(output.contains("✓ Discard reads with internal matches"));
        assert!(output.contains("✗ Keep reads with identical labels only"));
        assert!(output.contains("✓ Discard contaminated reads (≥75%)"));
    }

    #[test]
    fn test_display_compact() {
        let mut strat = FilterStrategy::new();
        strat.enable_exact();
        strat.enable_clean();

        assert_eq!(strat.display_compact(), "exact, clean");

        let empty = FilterStrategy::new();
        assert_eq!(empty.display_compact(), "none");
    }
}
