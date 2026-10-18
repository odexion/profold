//! Input classification and structure retrieval.
//!
//! One text field takes everything, so the classifier here decides whether a
//! paste is a file on disk, a PDB entry, a UniProt accession, a PubChem
//! compound, a sequence, or words to search for.

use std::io::Read;
use std::path::PathBuf;
use std::sync::mpsc::Sender;
use std::sync::Arc;
use std::time::{Duration, Instant};

use crate::formats::{self, parse_any};
use crate::job::{Event, Outcome, Progress};
use crate::protein::{Source, Structure};
use anyhow::{anyhow, Context, Result};

#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Input {
    File(PathBuf),
    Sequence(String),
    Pdb(String),
    UniProt(String),
    /// A PubChem compound, typed as `CID 2519`.
    Cid(u64),
    /// Words to search for by name.
    Query(String),
    /// What looks like a path, to a file that is not there.
    Missing(String),
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
            Input::Cid(_) => "PUBCHEM",
            Input::Query(_) => "SEARCH",
            Input::Missing(_) => "FILE",
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
            Input::Cid(cid) => format!("CID {cid}"),
            Input::Query(_) => "PDB · AlphaFold DB · PubChem".into(),
            Input::Missing(_) => "no such file".into(),
            Input::Empty => String::new(),
            Input::Unknown => "not a file, identifier, name or sequence".into(),
        }
    }

    /// The primary button's word in the viewer, where a sequence is looked
    /// up rather than folded.
    pub fn action(&self) -> &'static str {
        match self {
            Input::File(_) => "Open",
            Input::Sequence(_) | Input::Query(_) | Input::Empty | Input::Unknown => "Search",
            Input::Missing(_) => "Open",
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

/// Whether a word ends like a structure file's name.
fn has_extension(text: &str) -> bool {
    let name = text.strip_suffix(".gz").unwrap_or(text);
    name.rsplit_once('.').is_some_and(|(stem, ext)| {
        !stem.is_empty() && formats::Format::from_name(&format!("x.{ext}")).is_some()
    })
}

/// `CID 2519`, `cid:2519` or `CID2519`.
fn as_cid(text: &str) -> Option<u64> {
    let lower = text.trim().to_ascii_lowercase();
    let rest = lower.strip_prefix("cid")?;
    let rest = rest.trim_start_matches([':', ' ', '#']);
    rest.parse().ok().filter(|n| *n > 0)
}

/// Residue letters, if that is all the text is. A FASTA header settles it;
/// otherwise words are taken for a name unless they are long runs of capitals
/// or one long run of letters, the way sequences are written.
pub fn as_sequence(raw: &str, loose: bool) -> Option<String> {
    let trimmed = raw.trim();
    let seq = clean_sequence(trimmed);
    if seq.len() < 12 || !seq.chars().all(|c| AA.contains(c)) {
        return None;
    }
    if loose || trimmed.starts_with('>') {
        return Some(seq);
    }
    let tokens: Vec<&str> = trimmed
        .split_whitespace()
        .filter(|t| t.chars().any(|c| c.is_ascii_alphabetic()))
        .collect();
    let written_as_sequence = match tokens.as_slice() {
        [one] => one.len() >= 16 || one.chars().all(|c| !c.is_ascii_lowercase()),
        // Blocks of ten, as GenBank and many papers print them.
        many => many.iter().all(|t| t.len() >= 10),
    };
    written_as_sequence.then_some(seq)
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
    if !trimmed.contains(char::is_whitespace)
        && (trimmed.contains('/') || trimmed.starts_with('~') || has_extension(trimmed))
    {
        return Input::Missing(trimmed.to_string());
    }
    if let Some(cid) = as_cid(trimmed) {
        return Input::Cid(cid);
    }
    if let Some(seq) = as_sequence(trimmed, false) {
        return Input::Sequence(seq);
    }
    if trimmed.chars().filter(|c| c.is_alphanumeric()).count() >= 2 && !trimmed.contains('\n') {
        return Input::Query(trimmed.split_whitespace().collect::<Vec<_>>().join(" "));
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
    // The model's expected errors, which a plot shows; the model is worth
    // having without them.
    if let Some(url) = entry.get("paeDocUrl").and_then(|v| v.as_str()) {
        let _ = tx.send(Event::Stage("Fetching the predicted aligned error".into()));
        match get(url, Some(tx)).map_err(anyhow::Error::from).and_then(|t| crate::protein::Pae::from_json(&t)) {
            Ok(pae) if pae.n == s.residues.iter().filter(|r| !r.outside_polymer()).count() => s.pae = Some(Arc::new(pae)),
            Ok(_) => {
                let _ = tx.send(Event::Log("The predicted aligned error does not match the model; left out.".into()));
            }
            Err(e) => {
                let _ = tx.send(Event::Log(format!("No predicted aligned error ({e:#}).")));
            }
        }
    }
    Ok(s)
}

/// A compound from PubChem, in 3D where PubChem has a conformer for it and
/// as its flat drawing where it does not (proteins and other large ones).
pub fn fetch_pubchem(cid: u64, tx: &Sender<Event>) -> Result<Structure> {
    let _ = tx.send(Event::Stage(format!("Fetching CID {cid} from PubChem")));
    let base = format!("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}");
    let text = match get(&format!("{base}/SDF?record_type=3d"), Some(tx)) {
        Ok(text) => text,
        Err(_) => {
            let _ = tx.send(Event::Log(format!(
                "PubChem has no 3D conformer for CID {cid}; showing its 2D drawing."
            )));
            get(&format!("{base}/SDF?record_type=2d"), Some(tx))
                .map_err(|e| anyhow!("CID {cid} not found in PubChem ({e})"))?
        }
    };
    let _ = tx.send(Event::Stage("Reading structure".into()));
    let mut s = formats::parse(&text, formats::Format::Sdf, format!("CID {cid}"), Source::PubChem)?;
    if s.atoms.is_empty() {
        anyhow::bail!("no atoms in CID {cid}");
    }
    s.format = Some(formats::Format::Sdf.label());
    if let Some(title) = crate::search::compound_title(cid) {
        s.title = title;
    }
    Ok(s)
}

/// A UniProt entry's sequence, for folding.
#[cfg_attr(not(feature = "simplefold"), allow(dead_code))]
pub fn fetch_uniprot_sequence(acc: &str) -> Result<String> {
    let fasta = get(&format!("https://rest.uniprot.org/uniprotkb/{acc}.fasta"), None)?;
    let seq = clean_sequence(&fasta);
    if seq.is_empty() {
        anyhow::bail!("no sequence for {acc}");
    }
    Ok(seq)
}

/// Something to fetch: an entry in the PDB, an AlphaFold model, or a
/// compound in PubChem.
#[derive(Clone, Debug, PartialEq)]
pub enum Retrieval {
    Pdb(String),
    AlphaFold(String),
    PubChem(u64),
}

impl Retrieval {
    pub fn run(self, tx: &Sender<Event>) -> Result<Outcome> {
        match self {
            Retrieval::Pdb(id) => fetch_pdb(&id, tx),
            Retrieval::AlphaFold(acc) => fetch_alphafold(&acc, tx),
            Retrieval::PubChem(cid) => fetch_pubchem(cid, tx),
        }
        .map(Outcome::Structure)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn words_are_searched_and_residues_are_not() {
        assert_eq!(classify("green fluorescent protein"), Input::Query("green fluorescent protein".into()));
        assert_eq!(classify("caffeine"), Input::Query("caffeine".into()));
        assert_eq!(classify("samples/nothing.pdb"), Input::Missing("samples/nothing.pdb".into()));
        assert_eq!(classify("melt.lammpstrj"), Input::Missing("melt.lammpstrj".into()));
        assert_eq!(classify("GREEN FLUORESCENT PROTEIN"), Input::Query("GREEN FLUORESCENT PROTEIN".into()));
        assert_eq!(classify("1gfl"), Input::Pdb("1GFL".into()));
        assert_eq!(classify("P42212"), Input::UniProt("P42212".into()));
        assert_eq!(classify("CID 2519"), Input::Cid(2519));
        assert_eq!(classify("cid:2244"), Input::Cid(2244));
        let crambin = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN";
        assert_eq!(classify(crambin), Input::Sequence(crambin.into()));
        assert_eq!(classify(&format!(">crambin\n{crambin}")), Input::Sequence(crambin.into()));
        // GenBank's blocks of ten.
        assert!(matches!(classify("ttccpsivar snfnvcrlpg tpeaicatyt"), Input::Sequence(_)));
        // A fold tab takes words that could be residues as residues.
        assert!(as_sequence("green fluorescent protein", true).is_some());
        assert!(as_sequence("green fluorescent protein", false).is_none());
    }
}
