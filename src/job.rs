//! Shared reporting for the long-running jobs.
//!
//! Fetching and folding both take long enough to need saying so. They report
//! through one channel, and where the underlying tool tells us how far along it
//! is — curl's transfer table, a tqdm bar, an HTTP content length — that is
//! turned into a real fraction rather than another line of scrolling text.

use crate::protein::Structure;

#[derive(Clone, Debug, Default)]
pub struct Progress {
    pub text: String,
    /// 0..=1 when the job knows how far along it is.
    pub fraction: Option<f32>,
}

/// What a job produced. A UniProt accession with no deposited model resolves
/// to a sequence rather than a structure, which the app then folds.
pub enum Outcome {
    Structure(Structure),
    Sequence(String),
    /// An experimental structure and a prediction folded from its own
    /// sequence, already superimposed.
    #[cfg_attr(not(feature = "simplefold"), allow(dead_code))]
    Comparison {
        reference: Structure,
        prediction: Structure,
        rmsd: f32,
        matched: usize,
    },
}

pub enum Event {
    /// The phase the job has reached, in our own words.
    Stage(String),
    /// A line of output worth keeping.
    #[cfg_attr(not(feature = "simplefold"), allow(dead_code))]
    Log(String),
    /// Transient status, replacing whatever came before it.
    Progress(Progress),
    Finished(anyhow::Result<Outcome>),
}

impl Progress {
    #[cfg_attr(not(feature = "simplefold"), allow(dead_code))]
    pub fn indeterminate(text: impl Into<String>) -> Self {
        Progress {
            text: text.into(),
            fraction: None,
        }
    }

    pub fn bytes(read: u64, total: Option<u64>) -> Self {
        Progress {
            text: match total {
                Some(t) => format!("{} of {}", size(read), size(t)),
                None => size(read),
            },
            fraction: total.filter(|t| *t > 0).map(|t| read as f32 / t as f32),
        }
    }
}

pub fn size(bytes: u64) -> String {
    const UNITS: [&str; 4] = ["B", "KB", "MB", "GB"];
    let mut v = bytes as f64;
    let mut unit = 0;
    while v >= 1024.0 && unit < UNITS.len() - 1 {
        v /= 1024.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{} B", bytes)
    } else {
        format!("{v:.1} {}", UNITS[unit])
    }
}

/// Recover a percentage from a subprocess's own progress display.
#[cfg_attr(not(feature = "simplefold"), allow(dead_code))]
pub fn parse(line: &str) -> Progress {
    let trimmed = line.trim();

    // tqdm, and anything else that writes ` 42%|###...`.
    if let Some(bar) = trimmed.find("%|") {
        let digits: String = trimmed[..bar]
            .chars()
            .rev()
            .take_while(|c| c.is_ascii_digit())
            .collect();
        if let Ok(pct) = digits.chars().rev().collect::<String>().parse::<f32>() {
            let detail = trimmed
                .rsplit('|')
                .next()
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .unwrap_or(trimmed);
            return Progress {
                text: detail.to_string(),
                fraction: Some((pct / 100.0).clamp(0.0, 1.0)),
            };
        }
    }

    // curl's transfer table: the percentage leads, and the three trailing
    // clock fields are what make the row recognisable.
    let f: Vec<&str> = trimmed.split_whitespace().collect();
    if f.len() >= 10
        && f[f.len() - 2].matches(':').count() == 2
        && f[f.len() - 3].matches(':').count() == 2
    {
        if let Ok(pct) = f[0].parse::<f32>() {
            return Progress {
                text: format!("{} of {}", f[3], f[1]),
                fraction: Some((pct / 100.0).clamp(0.0, 1.0)),
            };
        }
    }

    Progress::indeterminate(trimmed)
}
