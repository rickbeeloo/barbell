use crate::search::searcher::BarbellMatch;
use sassy::search::Strand;

pub struct LabelConfig {
    include_label: bool,
    include_orientation: bool,
    include_flank: bool,
    sort_labels: bool,
}

impl LabelConfig {
    pub fn new(
        include_label: bool,
        include_orientation: bool,
        include_flank: bool,
        sort_labels: bool,
    ) -> Self {
        Self {
            include_label,
            include_orientation,
            include_flank,
            sort_labels,
        }
    }

    fn create_label(&self, annotations: &[BarbellMatch]) -> String {
        if !self.include_label {
            return "none".to_string();
        }

        let mut label_parts: Vec<String> = annotations
            .iter()
            .filter_map(|m| {
                let label = m.label.clone();

                // Skip if it's a flank and we don't want flanks
                if !self.include_flank && label.contains("flank") {
                    return None;
                }

                let mut result = label;

                if self.include_orientation {
                    let ori = match m.strand {
                        Strand::Fwd => "fw",
                        Strand::Rc => "rc",
                    };
                    result = format!("{}_{}", result, ori);
                }

                Some(result)
            })
            .collect();

        if label_parts.is_empty() {
            "none".to_string()
        } else if self.sort_labels {
            label_parts.sort();
            label_parts.join("__")
        } else {
            label_parts.join("__")
        }
    }
}
