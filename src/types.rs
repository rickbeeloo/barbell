use serde::{Serialize, Deserialize};

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum Orientation {
    Forward,
    ReverseComplement,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum MatchType {
    Fbarcode, 
    Rbarcode,
    Flank,  
} 

#[derive(Debug, Eq, PartialEq, Clone, Serialize, Deserialize)]
pub enum CutDirection {
    Before, // < cut at match.start
    After,  // > cut at match.end
}


// Records where to cut, in what direction, and also with id it belongs in case of paired cuts
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Cut {
    pub group_id: usize,
    pub direction: CutDirection,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Match {
    #[serde(default)]
    pub read: Option<String>,
    #[serde(serialize_with = "serialize_label", deserialize_with = "deserialize_label")]
    pub label: EncodedMatchStr,
    pub start: usize,
    pub end: usize,
    pub edit_dist: Option<i32>,
    #[serde(default)]
    pub read_len: Option<usize>,
    pub rel_dist_to_end: isize,
    #[serde(default)]
    pub record_set_idx: Option<usize>,
    #[serde(default)]
    pub record_idx: Option<usize>,
    #[serde(default)]
    #[serde(serialize_with = "serialize_cuts", deserialize_with = "deserialize_cuts")] //, skip_serializing_if = "Option::is_none"
    pub cuts: Option<Vec<Cut>>,
}

// We can implement Eq because we're using approximate equality for floats
// impl Eq for Match {}

impl Match {
    
    pub fn new(
        label: EncodedMatchStr,
        start: usize,
        end: usize,
        edit_dist: Option<i32>,
        rel_dist_to_end: isize, // negative for close to the right end
    ) -> Self {
        Self {
            label,
            start,
            end,
            edit_dist,
            rel_dist_to_end,
            // Below only modified before writing to csv so we can serialize/deserialize easily
            cuts: None, 
            read: None, 
            read_len: None, 
            record_set_idx: None, 
            record_idx: None, 
        }
    }

    pub fn add_read_info(&mut self, read: String, read_len: usize, record_set_idx: usize, record_idx: usize) {
        self.read = Some(read);
        self.read_len = Some(read_len);
        self.record_set_idx = Some(record_set_idx);
        self.record_idx = Some(record_idx);
    }
}

fn deserialize_cuts<'de, D>(deserializer: D) -> Result<Option<Vec<Cut>>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let s: String = String::deserialize(deserializer)?;
    
    if s.is_empty() || s == "-" {
        return Ok(None);
    }
    
    let cuts = if s.contains(',') {
        s.split(',')
            .map(str::trim)
            .flat_map(Cut::from_string)
            .collect::<Vec<_>>()
    } else {
        Cut::from_string(&s).into_iter().collect()
    };
    
    Ok(Some(cuts))
}

fn serialize_cuts<S>(cuts: &Option<Vec<Cut>>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    match cuts {
        None => serializer.serialize_str(""),
        Some(cuts) if cuts.is_empty() => serializer.serialize_str(""),
        Some(cuts) => {
            let cuts_str = cuts.iter()
                .map(|cut| cut.to_string())
                .collect::<Vec<_>>()
                .join(",");
            serializer.serialize_str(&cuts_str)
        }
    }
}

// Function to encode match information
#[derive(Debug, Eq, PartialEq, Clone, Serialize, Deserialize)]
pub struct EncodedMatchStr {
    pub match_type: MatchType,
    pub orientation: Orientation,
    pub label: Option<String>,
}

// This is only for  the label of the match
impl EncodedMatchStr {
    
    pub fn new(match_type: MatchType, orientation: Orientation, label: Option<String>) -> Self {
        Self { match_type, orientation, label }
    }

    pub fn stringify(&self) -> String {
        let label = if let Some(label) = &self.label {
            label.clone()
        } else {
            "Flank".to_string()
        };

        let ori  = match self.orientation {
            Orientation::Forward => "fw",
            Orientation::ReverseComplement => "rc",
        };

        let match_type = match self.match_type {
            MatchType::Fbarcode => "Fbarcode",
            MatchType::Rbarcode => "Rbarcode",
            MatchType::Flank => "Flank",
        };

        format!("{}#{}#{}", label, ori, match_type)
    }

    pub fn unstringify(s: &str) -> Self {
        let parts: Vec<&str> = s.split('#').collect();
        let label = parts[0].to_string();
        // Check label if label is flank we use label None, else the string
        let label = if label == "Flank" {
            None
        } else {
            Some(label)
        };

        let ori = match parts[1] {
            "fw" => Orientation::Forward,
            "rc" => Orientation::ReverseComplement,
            _ => panic!("Invalid orientation"),
        };
        let match_type = match parts[2] {
            "Fbarcode" => MatchType::Fbarcode,
            "Rbarcode" => MatchType::Rbarcode,
            "Flank" => MatchType::Flank,
            _ => panic!("Invalid match type"),
        };
        Self::new(match_type, ori, label)
    }
}

// Add these serialization functions
fn serialize_label<S>(label: &EncodedMatchStr, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    serializer.serialize_str(&label.stringify())
}

fn deserialize_label<'de, D>(deserializer: D) -> Result<EncodedMatchStr, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    Ok(EncodedMatchStr::unstringify(&s))
}

fn serialize_log_prob<S>(val: &Option<u8>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    match val {
        Some(x) => {
            serializer.serialize_some(x)
        }
        None => serializer.serialize_none(),
    }
}