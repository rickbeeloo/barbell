use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, SystemTime, UNIX_EPOCH};

#[derive(Clone, Copy)]
pub(crate) struct ProgressSpec {
    pub prefix: &'static str,
    pub color: &'static str,
    pub tick_ms: u64,
}

pub(crate) const ANNOTATION_PROGRESS_SPECS: [ProgressSpec; 3] = [
    ProgressSpec {
        prefix: "Total:",
        color: "cyan",
        tick_ms: 100,
    },
    ProgressSpec {
        prefix: "Kept:",
        color: "green",
        tick_ms: 120,
    },
    ProgressSpec {
        prefix: "Dropped:",
        color: "red",
        tick_ms: 140,
    },
];

pub(crate) const FILTER_PROGRESS_SPECS: [ProgressSpec; 3] = [
    ProgressSpec {
        prefix: "Total:",
        color: "cyan",
        tick_ms: 100,
    },
    ProgressSpec {
        prefix: "Kept:",
        color: "green",
        tick_ms: 120,
    },
    ProgressSpec {
        prefix: "Dropped:",
        color: "red",
        tick_ms: 140,
    },
];

pub(crate) const TRIM_PROGRESS_SPECS: [ProgressSpec; 4] = [
    ProgressSpec {
        prefix: "Total:",
        color: "cyan",
        tick_ms: 100,
    },
    ProgressSpec {
        prefix: "Kept:",
        color: "green",
        tick_ms: 120,
    },
    ProgressSpec {
        prefix: "Kept split:",
        color: "green",
        tick_ms: 140,
    },
    ProgressSpec {
        prefix: "Failed:",
        color: "red",
        tick_ms: 160,
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
    log: Option<ProgressLog>,
}

struct ProgressLog {
    path: PathBuf,
    step: String,
}

impl ProgressLog {
    fn write(&self, counts: &[AtomicUsize], specs: &[ProgressSpec]) -> Result<(), String> {
        let file = File::create(&self.path)
            .map_err(|e| format!("Failed to create log file '{}': {e}", self.path.display()))?;

        let mut w = BufWriter::new(file);
        writeln!(w, "step\tmetric\tcount")
            .map_err(|_| format!("Failed to write progress log '{}'", self.path.display()))?;

        for (count, spec) in counts.iter().zip(specs.iter()) {
            writeln!(
                w,
                "{}\t{}\t{}",
                self.step,
                spec.prefix,
                count.load(Ordering::Relaxed)
            )
            .map_err(|_| format!("Failed to write progress log '{}'", self.path.display()))?;
        }

        w.flush()
            .map_err(|_| format!("Failed to flush progress log '{}'", self.path.display()))
    }
}

impl ProgressTracker {
    pub(crate) fn new(specs: &[ProgressSpec]) -> Self {
        Self::new_inner(specs, None)
    }

    pub(crate) fn new_with_logging(
        specs: &[ProgressSpec],
        step: impl Into<String>,
        log_dir: impl AsRef<Path>,
    ) -> Self {
        let step = step.into();
        let ts = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_millis();
        let path = log_dir.as_ref().join(format!("{step}.{ts}.log"));
        Self::new_inner(specs, Some(ProgressLog { path, step }))
    }

    fn new_inner(specs: &[ProgressSpec], log: Option<ProgressLog>) -> Self {
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
            log,
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

        if let Some(log) = &self.log
            && let Err(e) = log.write(&self.counts, &self.specs)
        {
            self.print_error(e);
        }

        for (bar, count) in self.bars.iter().zip(self.counts.iter()) {
            let count = count.load(Ordering::Relaxed);
            // Prefix already shows the label (e.g. "Total:"), so don't repeat it in the message.
            bar.finish_with_message(format!("{count} {unit}"));
        }
    }
}
