//! Input classification and structure retrieval.
//!
//! One text field takes everything, so the classifier here decides whether a
//! paste is a file on disk, a PDB entry, a UniProt accession, or a sequence
//! to fold.

use std::io::Read;
use std::path::PathBuf;
use std::sync::mpsc::Sender;
use std::time::{Duration, Instant};

use crate::formats::parse_any;
use crate::job::{Event, Outcome, Progress};
use crate::protein::{Source, Structure};
use anyhow::{anyhow, Context, Result};

#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Input {
    File(PathBuf),
    Sequence(String),
    Pdb(String),
    UniProt(String),
    Empty,
    Unknown,
}

impl Input {
    pub fn kind(&self) -> &'static str {
        match self {
            Input::File(_) => "FILE",
            Input::Sequence(_) => "SEQUENCE",
            Input::Pdb(_) => "PDB ID",
            Input::UniProt(_) => "UNIPROT",
            Input::Empty => "",
            Input::Unknown => "UNRECOGNISED",
        }
    }

    pub fn detail(&self) -> String {
        match self {
            Input::File(p) => p
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default(),
            Input::Sequence(s) => format!("{} aa", s.len()),
            Input::Pdb(id) => id.clone(),
            Input::UniProt(id) => id.clone(),
            Input::Empty => String::new(),
            Input::Unknown => "not a file, PDB ID, accession or sequence".into(),
        }
    }

    pub fn action(&self) -> &'static str {
        match self {
            Input::File(_) => "Open",
            Input::Sequence(_) => "Fold",
            _ => "Fetch",
        }
    }
}

const AA: &str = "ACDEFGHIKLMNPQRSTVWYBXZUO";

/// Strip FASTA headers, digits and whitespace, leaving bare residue letters.
pub fn clean_sequence(raw: &str) -> String {
    raw.lines()
        .filter(|l| !l.trim_start().starts_with('>'))
        .flat_map(|l| l.chars())
        .filter(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_ascii_uppercase())
        .collect()
}

fn is_pdb_id(s: &str) -> bool {
    let b = s.as_bytes();
    b.len() == 4 && b[0].is_ascii_digit() && b[1..].iter().all(|c| c.is_ascii_alphanumeric())
}

fn is_uniprot(s: &str) -> bool {
    let b = s.as_bytes();
    let digit = |i: usize| b[i].is_ascii_digit();
    let alnum = |i: usize| b[i].is_ascii_alphanumeric();
    let upper = |i: usize| b[i].is_ascii_uppercase();

    // [OPQ][0-9][A-Z0-9]{3}[0-9]
    if b.len() == 6 && matches!(b[0], b'O' | b'P' | b'Q') {
        if digit(1) && (2..5).all(alnum) && digit(5) {
            return true;
        }
    }
    // [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
    if (b.len() == 6 || b.len() == 10) && upper(0) && !matches!(b[0], b'O' | b'P' | b'Q') {
        if digit(1) {
            let ok = (0..(b.len() - 2) / 4).all(|k| {
                let o = 2 + k * 4;
                upper(o) && alnum(o + 1) && alnum(o + 2) && digit(o + 3)
            });
            if ok {
                return true;
            }
        }
    }
    false
}

/// A path typed or pasted into the box, with `~` for home and the quotes a
/// copied path sometimes carries.
fn as_file(text: &str) -> Option<PathBuf> {
    let text = text.trim_matches(|c: char| c == '"' || c == '\'');
    if !(text.contains('/') || text.contains('.') || text.starts_with('~')) {
        return None;
    }
    let path = match text.strip_prefix("~/") {
        Some(rest) => dirs::home_dir()?.join(rest),
        None => PathBuf::from(text),
    };
    path.is_file().then_some(path)
}

pub fn classify(raw: &str) -> Input {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Input::Empty;
    }
    if !trimmed.contains('\n') {
        if let Some(path) = as_file(trimmed) {
            return Input::File(path);
        }
    }

    // A single bare token can be an identifier.
    if !trimmed.contains(char::is_whitespace) && !trimmed.starts_with('>') {
        let token = trimmed.trim_matches(|c: char| c == '"' || c == '\'');
        let upper = token.to_ascii_uppercase();

        // AlphaFold DB filenames, e.g. AF-P69905-F1.
        if let Some(rest) = upper.strip_prefix("AF-") {
            let acc = rest.split('-').next().unwrap_or("");
            if is_uniprot(acc) {
                return Input::UniProt(acc.to_string());
            }
        }
        if is_pdb_id(&upper) {
            return Input::Pdb(upper);
        }
        if is_uniprot(&upper) {
            return Input::UniProt(upper);
        }
    }

    let seq = clean_sequence(trimmed);
    if seq.len() >= 12 && seq.chars().all(|c| AA.contains(c)) {
        return Input::Sequence(seq);
    }
    Input::Unknown
}

// ---------------------------------------------------------------------------
// Network
// ---------------------------------------------------------------------------

/// Structures run to tens of megabytes, well past ureq's 10 MB default.
const MAX_BODY: u64 = 256 * 1024 * 1024;

/// Fetch a document, reporting progress against the content length when the
/// server declares one.
fn get(url: &str, tx: Option<&Sender<Event>>) -> Result<String> {
    let mut resp = ureq::get(url)
        // UniProt serves gzip without a `content-encoding` header, which no
        // client can decode; ask for plain bytes instead.
        .header("accept-encoding", "identity")
        .header("user-agent", "profold/0.1")
        .call()
        .with_context(|| format!("requesting {url}"))?;

    let total: Option<u64> = resp
        .headers()
        .get("content-length")
        .and_then(|v| v.to_str().ok())
        .and_then(|v| v.parse().ok());

    let mut reader = resp.body_mut().with_config().limit(MAX_BODY).reader();
    let mut out: Vec<u8> = Vec::with_capacity(total.unwrap_or(64 * 1024) as usize);
    let mut buf = vec![0u8; 64 * 1024];
    let mut read = 0u64;
    let mut last = Instant::now();

    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        out.extend_from_slice(&buf[..n]);
        read += n as u64;
        // Throttled, so a fast transfer does not flood the channel.
        if let Some(tx) = tx {
            if last.elapsed() > Duration::from_millis(80) {
                last = Instant::now();
                let _ = tx.send(Event::Progress(Progress::bytes(read, total)));
            }
        }
    }
    if let Some(tx) = tx {
        let _ = tx.send(Event::Progress(Progress::bytes(read, total.or(Some(read)))));
    }
    Ok(String::from_utf8_lossy(&out).into_owned())
}

pub fn fetch_pdb(id: &str, tx: &Sender<Event>) -> Result<Structure> {
    let _ = tx.send(Event::Stage(format!("Fetching {id} from the PDB")));
    let url = format!("https://files.rcsb.org/download/{id}.cif");
    let text = get(&url, Some(tx)).map_err(|e| anyhow!("{id} not found in the PDB ({e})"))?;
    let _ = tx.send(Event::Stage("Reading structure".into()));
    parse_any(&text, id.to_string(), Source::Pdb)
}

pub fn fetch_alphafold(acc: &str, tx: &Sender<Event>) -> Result<Structure> {
    let _ = tx.send(Event::Stage(format!("Fetching {acc} from AlphaFold DB")));
    let meta = get(&format!("https://alphafold.ebi.ac.uk/api/prediction/{acc}"), None)
        .map_err(|e| anyhow!("no AlphaFold model for {acc} ({e})"))?;
    let json: serde_json::Value = serde_json::from_str(&meta)?;
    let entry = json
        .get(0)
        .ok_or_else(|| anyhow!("no AlphaFold model for {acc}"))?;

    let url = entry
        .get("cifUrl")
        .or_else(|| entry.get("pdbUrl"))
        .and_then(|v| v.as_str())
        .ok_or_else(|| anyhow!("AlphaFold response had no structure URL"))?;

    let text = get(url, Some(tx))?;
    let _ = tx.send(Event::Stage("Reading structure".into()));
    let mut s = parse_any(&text, acc.to_string(), Source::AlphaFold)?;
    if let Some(desc) = entry.get("uniprotDescription").and_then(|v| v.as_str()) {
        s.title = desc.to_string();
    }
    Ok(s)
}

/// Used when an accession has no AlphaFold model but can still be folded.
pub fn fetch_uniprot_sequence(acc: &str) -> Result<String> {
    let fasta = get(&format!("https://rest.uniprot.org/uniprotkb/{acc}.fasta"), None)?;
    let seq = clean_sequence(&fasta);
    if seq.is_empty() {
        anyhow::bail!("no sequence for {acc}");
    }
    Ok(seq)
}

/// Something to fetch: an entry in the PDB, or an AlphaFold model.
#[derive(Clone)]
pub enum Retrieval {
    Pdb(String),
    AlphaFold(String),
}

impl Retrieval {
    pub fn run(self, tx: &Sender<Event>) -> Result<Outcome> {
        match self {
            Retrieval::Pdb(id) => fetch_pdb(&id, tx).map(Outcome::Structure),
            Retrieval::AlphaFold(acc) => match fetch_alphafold(&acc, tx) {
                Ok(s) => Ok(Outcome::Structure(s)),
                // An accession without a model can still be folded, by a
                // plugin that folds.
                Err(e) => match fetch_uniprot_sequence(&acc) {
                    Ok(seq) => Ok(Outcome::Sequence(seq)),
                    Err(_) => Err(e),
                },
            },
        }
    }
}
