use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Duration;

#[derive(Clone, Copy)]
pub(crate) struct ProgressSpec {
    pub prefix: &'static str,
    pub color: &'static str,
    pub tick_ms: u64,
    pub finish_label: &'static str,
}

pub(crate) const ANNOTATION_PROGRESS_SPECS: [ProgressSpec; 3] = [
    ProgressSpec {
        prefix: "Total:",
        color: "cyan",
        tick_ms: 100,
        finish_label: "Done",
    },
    ProgressSpec {
        prefix: "Found:",
        color: "green",
        tick_ms: 120,
        finish_label: "Found",
    },
    ProgressSpec {
        prefix: "Missed:",
        color: "red",
        tick_ms: 140,
        finish_label: "Missed",
    },
];

pub(crate) const TRIM_PROGRESS_SPECS: [ProgressSpec; 4] = [
    ProgressSpec {
        prefix: "Total:",
        color: "cyan",
        tick_ms: 100,
        finish_label: "Total",
    },
    ProgressSpec {
        prefix: "Trimmed:",
        color: "green",
        tick_ms: 120,
        finish_label: "Trimmed",
    },
    ProgressSpec {
        prefix: "Trimmed split:",
        color: "green",
        tick_ms: 140,
        finish_label: "Trimmed split",
    },
    ProgressSpec {
        prefix: "Failed:",
        color: "red",
        tick_ms: 160,
        finish_label: "Failed",
    },
];

fn create_spinner_bar(multi_progress: &MultiProgress, spec: ProgressSpec) -> ProgressBar {
    let bar = multi_progress.add(ProgressBar::new_spinner());
    let template = format!(
        "{{spinner:.{color}}} {{prefix:.bold.white:<8}} {{msg:.bold.{color}:>6}} {{elapsed:.dim}}",
        color = spec.color
    );
    bar.set_style(
        ProgressStyle::with_template(&template)
            .unwrap()
            .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    bar.enable_steady_tick(Duration::from_millis(spec.tick_ms));
    bar.set_prefix(spec.prefix.to_string());
    bar
}

pub(crate) struct ProgressTracker {
    bars: Vec<ProgressBar>,
    counts: Vec<AtomicUsize>,
    specs: Vec<ProgressSpec>,
    error_msg: Mutex<Option<String>>,
}

impl ProgressTracker {
    pub(crate) fn new(specs: &[ProgressSpec]) -> Self {
        let multi_progress = MultiProgress::new();
        let bars = specs
            .iter()
            .map(|spec| create_spinner_bar(&multi_progress, *spec))
            .collect();
        let counts = specs.iter().map(|_| AtomicUsize::new(0)).collect();
        Self {
            bars,
            counts,
            specs: specs.to_vec(),
            error_msg: Mutex::new(None),
        }
    }

    #[inline(always)]
    pub(crate) fn add(&self, idx: usize, count: usize) {
        self.counts[idx].fetch_add(count, Ordering::Relaxed);
    }

    #[inline(always)]
    pub(crate) fn inc(&self, idx: usize) {
        self.add(idx, 1);
    }

    #[inline(always)]
    pub(crate) fn refresh(&self) {
        for (bar, count) in self.bars.iter().zip(self.counts.iter()) {
            bar.set_message(count.load(Ordering::Relaxed).to_string());
        }
    }

    pub(crate) fn store_error(&self, msg: impl Into<String>) {
        let msg = msg.into();
        self.bars[0].println(msg.clone());
        if let Ok(mut err) = self.error_msg.lock() {
            *err = Some(msg);
        }
    }

    pub(crate) fn print_error(&self, msg: impl Into<String>) {
        self.bars[0].println(msg.into());
    }

    pub(crate) fn take_error(&self) -> Option<String> {
        self.error_msg.lock().ok().and_then(|mut err| err.take())
    }

    pub(crate) fn clear(&self) {
        for bar in &self.bars {
            bar.finish_and_clear();
        }
    }

    pub(crate) fn finish(&self, unit: &str) {
        self.refresh();
        for ((bar, count), spec) in self
            .bars
            .iter()
            .zip(self.counts.iter())
            .zip(self.specs.iter())
        {
            let count = count.load(Ordering::Relaxed);
            bar.finish_with_message(format!("{}: {count} {unit}", spec.finish_label));
        }
    }
}
