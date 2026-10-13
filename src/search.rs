//! Searching by name, so nobody has to know an identifier by heart.
//!
//! A query goes to three places at once: the PDB for experimental structures,
//! UniProt for proteins with an AlphaFold model, and PubChem for everything
//! that is not a protein. Each database ranks by its own idea of relevance,
//! and the PDB's is loose: a search for green fluorescent protein turns up red
//! ones first. So the hits are scored again here, by how closely their names
//! match what was typed, and the database whose best hit is a match for the
//! name leads.

use std::sync::OnceLock;
use std::time::Duration;

use anyhow::{anyhow, Result};
use serde_json::{json, Value};

use crate::fetch::Retrieval;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    /// An experimental structure in the PDB.
    Pdb,
    /// A UniProt entry, opened as its AlphaFold model.
    AlphaFold,
    /// A compound in PubChem.
    Molecule,
}

impl Kind {
    /// What the group of hits is called in the list.
    pub fn heading(self) -> &'static str {
        match self {
            Kind::Pdb => "Structures · RCSB PDB",
            Kind::AlphaFold => "Predicted models · AlphaFold DB",
            Kind::Molecule => "Molecules · PubChem",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Hit {
    pub kind: Kind,
    /// The PDB ID, UniProt accession, or PubChem CID.
    pub id: String,
    pub title: String,
    /// Organism, method and resolution; or formula and mass.
    pub detail: String,
    score: f32,
    /// Whether its name is the query, give or take the boilerplate.
    exact: bool,
    /// Found as a synonym of the query rather than by its own name.
    alias: bool,
}

impl Hit {
    pub fn retrieval(&self) -> Retrieval {
        match self.kind {
            Kind::Pdb => Retrieval::Pdb(self.id.clone()),
            Kind::AlphaFold => Retrieval::AlphaFold(self.id.clone()),
            Kind::Molecule => Retrieval::PubChem(self.id.parse().unwrap_or(0)),
        }
    }

    #[cfg(test)]
    pub fn for_test(kind: Kind, id: &str, title: &str) -> Hit {
        Hit {
            kind,
            id: id.into(),
            title: title.into(),
            detail: String::new(),
            score: 0.0,
            exact: false,
            alias: false,
        }
    }

    /// The identifier as shown beside the title.
    pub fn label(&self) -> String {
        match self.kind {
            Kind::Molecule => format!("CID {}", self.id),
            _ => self.id.clone(),
        }
    }
}

const PDB_ROWS: usize = 30;
const PDB_EARLIEST: usize = 10;
const PDB_SHOWN: usize = 5;
const UNIPROT_SHOWN: usize = 3;
const MOLECULES_SHOWN: usize = 4;

/// Search all three databases for a name, best match first. A database that
/// fails is left out; only when all of them fail is it an error.
pub fn search(query: &str) -> Result<Vec<Hit>> {
    let query = &american(query.trim());
    let (pdb, uniprot, pubchem) = std::thread::scope(|s| {
        let pdb = s.spawn(|| search_pdb(query));
        let uniprot = s.spawn(|| search_uniprot(query));
        let pubchem = s.spawn(|| search_pubchem(query));
        (
            pdb.join().unwrap_or_else(|_| Err(anyhow!("PDB search failed"))),
            uniprot.join().unwrap_or_else(|_| Err(anyhow!("UniProt search failed"))),
            pubchem.join().unwrap_or_else(|_| Err(anyhow!("PubChem search failed"))),
        )
    });
    if let (Err(e), Err(_), Err(_)) = (&pdb, &uniprot, &pubchem) {
        return Err(anyhow!("search failed: {e:#}"));
    }
    let (pdb, uniprot, mut molecules) = (
        pdb.unwrap_or_default(),
        uniprot.unwrap_or_default(),
        pubchem.unwrap_or_default(),
    );
    // PubChem's synonyms are generous: vitamin C is ascorbic acid, but it
    // also has ubiquitin down as citric acid. A synonym that shares no word
    // with the query is believed only when no protein answers to the name.
    let protein = |hits: &[Hit]| hits.first().is_some_and(|h| h.exact);
    if protein(&pdb) || protein(&uniprot) {
        molecules.retain(|h| !h.alias);
    }
    Ok(arrange(pdb, uniprot, molecules))
}

/// Put the groups in order: the one whose best hit is named exactly what was
/// asked for leads, proteins before molecules when both are. Without an exact
/// match anywhere, structures lead, then models, then molecules.
fn arrange(pdb: Vec<Hit>, uniprot: Vec<Hit>, molecules: Vec<Hit>) -> Vec<Hit> {
    let exact = |hits: &[Hit]| hits.first().is_some_and(|h| h.exact);
    let mut groups = vec![pdb, uniprot, molecules];
    // A stable sort keeps the default order among equals.
    groups.sort_by_key(|g| !exact(g));
    groups.into_iter().flatten().collect()
}

/// Entries in the PDB whose chains match a sequence, closest first.
pub fn search_sequence(sequence: &str) -> Result<Vec<Hit>> {
    let body = json!({
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 0.1,
                "identity_cutoff": 0.3,
                "sequence_type": "protein",
                "value": sequence,
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": { "start": 0, "rows": 8 },
            "results_content_type": ["experimental"],
            "scoring_strategy": "sequence",
        }
    });
    let ids = rcsb_ids(&body)?;
    // The order is the search's: by alignment score.
    Ok(describe_entries(&ids)?
        .into_iter()
        .enumerate()
        .map(|(i, e)| e.hit(-(i as f32), false))
        .collect())
}

// ---------------------------------------------------------------------------
// The PDB
// ---------------------------------------------------------------------------

fn search_pdb(query: &str) -> Result<Vec<Hit>> {
    // Two pools: the text search's own idea of relevance, and the earliest
    // entries that match at all. The first structure of a protein is often
    // the one in the textbooks (1UBQ, 1LYZ, 1EMA), and among a thousand
    // later variants the relevance ranking alone never finds it.
    let request = |rows: usize, sort: Value| {
        json!({
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": { "value": query }
            },
            "return_type": "entry",
            "request_options": {
                "paginate": { "start": 0, "rows": rows },
                "results_content_type": ["experimental"],
                "sort": [sort],
            }
        })
    };
    let relevant = request(PDB_ROWS, json!({ "sort_by": "score", "direction": "desc" }));
    let earliest = request(
        PDB_EARLIEST,
        json!({ "sort_by": "rcsb_accession_info.initial_release_date", "direction": "asc" }),
    );
    let (relevant, earliest) = std::thread::scope(|s| {
        let a = s.spawn(|| rcsb_ids(&relevant));
        let b = s.spawn(|| rcsb_ids(&earliest));
        (
            a.join().unwrap_or_else(|_| Err(anyhow!("PDB search failed"))),
            b.join().unwrap_or_else(|_| Err(anyhow!("PDB search failed"))),
        )
    });
    let relevant = relevant?;
    let earliest = earliest.unwrap_or_default();
    let mut ids = relevant.clone();
    ids.extend(earliest.iter().filter(|id| !relevant.contains(id)).cloned());

    let entries = describe_entries(&ids)?;
    let wanted = words(query);
    let n = relevant.len().max(1) as f32;
    let mut hits: Vec<Hit> = entries
        .into_iter()
        .map(|e| {
            let (title_score, title_exact) = matching(&wanted, &subject_words(&e.title));
            let (name_score, name_exact) = e
                .names
                .iter()
                .map(|d| matching(&wanted, &words(d)))
                .fold((0.0f32, false), |a, b| (a.0.max(b.0), a.1 || b.1));
            let relevance = relevant
                .iter()
                .position(|id| *id == e.id)
                .map_or(0.0, |rank| 0.6 * (1.0 - rank as f32 / n));
            let classic = if earliest.contains(&e.id) { 0.7 } else { 0.0 };
            let score = title_score + name_score + relevance + classic + e.quality
                // Most who ask for a protein by name mean the human one.
                + if e.detail.starts_with("Homo sapiens") { 0.25 } else { 0.0 };
            e.hit(score, title_exact || name_exact)
        })
        .collect();
    hits.sort_by(|a, b| b.score.total_cmp(&a.score));
    hits.truncate(PDB_SHOWN);
    Ok(hits)
}

/// What the PDB says about an entry, enough to list and rank it.
struct Entry {
    id: String,
    title: String,
    /// Organism, method and resolution.
    detail: String,
    /// The names of its chains.
    names: Vec<String>,
    /// Better resolved structures are a little more likely to be the one
    /// wanted; it only ever settles near ties.
    quality: f32,
}

impl Entry {
    fn hit(self, score: f32, exact: bool) -> Hit {
        Hit {
            kind: Kind::Pdb,
            id: self.id,
            title: tidy_title(&self.title),
            detail: self.detail,
            score,
            exact,
            alias: false,
        }
    }
}

/// The identifiers a search returns. No hits comes back as an empty body.
fn rcsb_ids(body: &Value) -> Result<Vec<String>> {
    let reply = post_json("https://search.rcsb.org/rcsbsearch/v2/query", body)?;
    Ok(reply
        .get("result_set")
        .and_then(Value::as_array)
        .map(|rows| {
            rows.iter()
                .filter_map(|r| r.get("identifier").and_then(Value::as_str))
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default())
}

/// Titles, organisms and resolutions for a list of entries, in one request.
fn describe_entries(ids: &[String]) -> Result<Vec<Entry>> {
    if ids.is_empty() {
        return Ok(Vec::new());
    }
    let list = ids
        .iter()
        .map(|id| format!("\"{id}\""))
        .collect::<Vec<_>>()
        .join(",");
    let query = format!(
        "{{entries(entry_ids:[{list}]){{rcsb_id struct{{title}} \
         rcsb_entry_info{{resolution_combined experimental_method}} \
         polymer_entities{{rcsb_polymer_entity{{pdbx_description}} \
         rcsb_entity_source_organism{{scientific_name}}}}}}}}"
    );
    let reply = post_json("https://data.rcsb.org/graphql", &json!({ "query": query }))?;
    let entries = reply
        .pointer("/data/entries")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();

    let mut hits = Vec::new();
    // GraphQL does not keep the order asked for; the search's order matters.
    for id in ids {
        let Some(e) = entries
            .iter()
            .find(|e| e.get("rcsb_id").and_then(Value::as_str) == Some(id))
        else {
            continue;
        };
        let title = e
            .pointer("/struct/title")
            .and_then(Value::as_str)
            .unwrap_or("")
            .to_string();
        let info = e.get("rcsb_entry_info");
        let method = info
            .and_then(|i| i.get("experimental_method"))
            .and_then(Value::as_str)
            .unwrap_or("");
        let resolution = info
            .and_then(|i| i.get("resolution_combined"))
            .and_then(|r| r.get(0))
            .and_then(Value::as_f64);
        let entities = e
            .get("polymer_entities")
            .and_then(Value::as_array)
            .cloned()
            .unwrap_or_default();
        let names: Vec<String> = entities
            .iter()
            .filter_map(|p| p.pointer("/rcsb_polymer_entity/pdbx_description"))
            .filter_map(Value::as_str)
            .map(str::to_string)
            .collect();
        let organism = entities
            .iter()
            .filter_map(|p| p.pointer("/rcsb_entity_source_organism/0/scientific_name"))
            .filter_map(Value::as_str)
            .next()
            .map(sentence_organism);

        let method = match method {
            "X-ray" => "X-ray",
            "EM" => "cryo-EM",
            "NMR" => "NMR",
            "Multiple methods" => "several methods",
            other => other,
        };
        let mut detail = Vec::new();
        detail.extend(organism);
        detail.push(match resolution {
            Some(r) if method != "NMR" => format!("{method} {r:.1} Å"),
            _ => method.to_string(),
        });
        hits.push(Entry {
            id: id.clone(),
            title,
            detail: detail.join(" · "),
            names,
            quality: match resolution {
                Some(r) if r <= 2.0 => 0.3,
                Some(r) if r <= 3.0 => 0.15,
                _ => 0.0,
            },
        });
    }
    Ok(hits)
}

// ---------------------------------------------------------------------------
// UniProt
// ---------------------------------------------------------------------------

fn search_uniprot(query: &str) -> Result<Vec<Hit>> {
    let reply = agent()
        .get("https://rest.uniprot.org/uniprotkb/search")
        .query("query", format!("({query}) AND reviewed:true"))
        .query("fields", "accession,protein_name,organism_name,length")
        .query("size", "6")
        .query("format", "json")
        .header("accept-encoding", "identity")
        .header("user-agent", USER_AGENT)
        .call()?
        .body_mut()
        .with_config()
        .limit(8 * 1024 * 1024)
        .read_to_string()?;
    let json: Value = serde_json::from_str(&reply)?;
    let wanted = words(query);

    let mut hits: Vec<Hit> = json
        .get("results")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default()
        .iter()
        .enumerate()
        .filter_map(|(rank, r)| {
            let id = r.get("primaryAccession")?.as_str()?.to_string();
            let name = r
                .pointer("/proteinDescription/recommendedName/fullName/value")
                .or_else(|| r.pointer("/proteinDescription/submissionNames/0/fullName/value"))
                .and_then(Value::as_str)?
                .to_string();
            let organism = r
                .pointer("/organism/scientificName")
                .and_then(Value::as_str)
                .unwrap_or("");
            let length = r.pointer("/sequence/length").and_then(Value::as_u64);
            let (score, exact) = matching(&wanted, &words(&name));
            let mut detail = vec![organism.to_string()];
            detail.extend(length.map(|n| format!("{n} aa")));
            Some(Hit {
                kind: Kind::AlphaFold,
                id,
                title: name,
                detail: detail.join(" · "),
                score: score + 0.6 * (1.0 - rank as f32 / 6.0),
                exact,
                alias: false,
            })
        })
        .collect();
    hits.sort_by(|a, b| b.score.total_cmp(&a.score));
    hits.truncate(UNIPROT_SHOWN);
    Ok(hits)
}

// ---------------------------------------------------------------------------
// PubChem
// ---------------------------------------------------------------------------

const PUG: &str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug";

/// The compound a name stands for, when PubChem knows it as a synonym, and
/// then the rest of what NCBI's search finds for it. Mixtures (the "aspirin;
/// caffeine" of a cold remedy) are left out: they have no one structure.
fn search_pubchem(query: &str) -> Result<Vec<Hit>> {
    let (named, found) = std::thread::scope(|s| {
        let named = s.spawn(|| {
            properties(&format!("{PUG}/compound/name/{}/property/{PROPS}/JSON", encode(query)))
        });
        let found = s.spawn(|| -> Result<Vec<Hit>> {
            let reply = agent()
                .get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi")
                .query("db", "pccompound")
                .query("term", query)
                .query("retmode", "json")
                .query("retmax", "12")
                .query("sort", "relevance")
                .header("user-agent", USER_AGENT)
                .call()?
                .body_mut()
                .read_to_string()?;
            let json: Value = serde_json::from_str(&reply)?;
            let ids: Vec<&str> = json
                .pointer("/esearchresult/idlist")
                .and_then(Value::as_array)
                .map(|a| a.iter().filter_map(Value::as_str).collect())
                .unwrap_or_default();
            if ids.is_empty() {
                return Ok(Vec::new());
            }
            properties(&format!("{PUG}/compound/cid/{}/property/{PROPS}/JSON", ids.join(",")))
        });
        (
            named.join().unwrap_or_else(|_| Err(anyhow!("PubChem failed"))),
            found.join().unwrap_or_else(|_| Err(anyhow!("PubChem failed"))),
        )
    });
    if let (Err(e), Err(_)) = (&named, &found) {
        // PubChem says 404 for a name it does not know; that is no hits.
        if e.to_string().contains("404") {
            return Ok(Vec::new());
        }
    }
    let wanted = words(query);
    let mut hits = Vec::new();
    // A synonym match is the compound the name means, whatever its title.
    for mut h in named.unwrap_or_default().into_iter().take(1) {
        h.exact = true;
        h.alias = !words(&h.title).iter().any(|w| wanted.contains(w));
        h.score = 10.0;
        hits.push(h);
    }
    let found = found.unwrap_or_default();
    let n = found.len().max(1) as f32;
    let mut rest: Vec<Hit> = found
        .into_iter()
        .enumerate()
        .filter(|(_, h)| !h.title.contains(';') && !hits.iter().any(|k| k.id == h.id))
        .map(|(rank, mut h)| {
            let (score, exact) = matching(&wanted, &words(&h.title));
            h.exact = exact;
            h.score = score + 0.6 * (1.0 - rank as f32 / n);
            h
        })
        // Something whose name has none of the words asked for was found by
        // a synonym or a patent, and is noise in a list like this.
        .filter(|h| h.score > 0.6)
        .collect();
    rest.sort_by(|a, b| b.score.total_cmp(&a.score));
    hits.extend(rest);
    hits.truncate(MOLECULES_SHOWN);
    Ok(hits)
}

const PROPS: &str = "Title,MolecularFormula,MolecularWeight";

fn properties(url: &str) -> Result<Vec<Hit>> {
    let json = get_json(url)?;
    Ok(json
        .pointer("/PropertyTable/Properties")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default()
        .iter()
        .filter_map(|p| {
            let cid = p.get("CID")?.as_u64()?;
            let title = p.get("Title").and_then(Value::as_str).unwrap_or("").to_string();
            let formula = p.get("MolecularFormula").and_then(Value::as_str).unwrap_or("");
            // PubChem gives the weight as a string.
            let weight = p
                .get("MolecularWeight")
                .and_then(|w| w.as_str().map(str::to_string).or_else(|| w.as_f64().map(|f| f.to_string())))
                .unwrap_or_default();
            let mut detail = vec![formula_subscripts(formula)];
            if !weight.is_empty() {
                detail.push(format!("{weight} g/mol"));
            }
            Some(Hit {
                kind: Kind::Molecule,
                id: cid.to_string(),
                title: if title.is_empty() { format!("CID {cid}") } else { title },
                detail: detail.join(" · "),
                score: 0.0,
                exact: false,
                alias: false,
            })
        })
        .collect())
}

/// The name, formula and weight of one compound, for a structure's title.
pub fn compound_title(cid: u64) -> Option<String> {
    properties(&format!("{PUG}/compound/cid/{cid}/property/{PROPS}/JSON"))
        .ok()?
        .into_iter()
        .next()
        .map(|h| h.title)
}

// ---------------------------------------------------------------------------
// Matching names
// ---------------------------------------------------------------------------

/// Words that say nothing about which molecule is meant.
const STOP: [&str; 22] = [
    "a", "an", "the", "of", "in", "from", "with", "and", "for", "by", "to", "at",
    "structure", "structures", "crystal", "crystallographic", "solution", "nmr",
    "x", "ray", "chain", "subunit",
];

/// The databases index American spellings, and the PDB's text search finds
/// little for haemoglobin that it would for hemoglobin.
fn american(query: &str) -> String {
    const SPELLINGS: [(&str, &str); 9] = [
        ("haem", "hem"),
        ("aemia", "emia"),
        ("oestr", "estr"),
        ("oedem", "edem"),
        ("foet", "fet"),
        ("sulph", "sulf"),
        ("caesium", "cesium"),
        ("aluminium", "aluminum"),
        ("ischaem", "ischem"),
    ];
    query
        .split(' ')
        .map(|word| {
            let mut w = word.to_string();
            for (from, to) in SPELLINGS {
                let lower = w.to_lowercase();
                if let Some(at) = lower.find(from) {
                    w.replace_range(at..at + from.len(), to);
                }
            }
            w
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// Lower-case words, the filler left out, and British spellings folded into
/// American so haemoglobin matches hemoglobin.
fn words(text: &str) -> Vec<String> {
    text.split(|c: char| !c.is_alphanumeric())
        .filter(|w| !w.is_empty())
        .map(|w| american(&w.to_lowercase()))
        .filter(|w| !STOP.contains(&w.as_str()))
        .collect()
}

/// The words of a PDB title, after the experiment and the conditions are
/// taken off: "Crystal structure of X in complex with Y at 2 Å" is about X.
fn subject_words(title: &str) -> Vec<String> {
    let lower = title.to_lowercase();
    let cut = [" in complex", " complexed", " bound to", " with ", " at ", ", "]
        .iter()
        .filter_map(|m| lower.find(m))
        .min()
        .unwrap_or(lower.len());
    words(&title[..cut.min(title.len())])
}

/// How well a name answers a query: 3 when it is the query, 2 less a little
/// per extra word when it contains it, and the share of the query's words it
/// has otherwise. The flag says whether it was the query exactly.
fn matching(wanted: &[String], name: &[String]) -> (f32, bool) {
    if wanted.is_empty() || name.is_empty() {
        return (0.0, false);
    }
    // Long words also count inside longer ones: hemoglobin in deoxyhemoglobin.
    let has = |w: &String| name.iter().any(|n| n == w || (w.len() > 5 && n.contains(w.as_str())));
    let found = wanted.iter().filter(|w| has(w)).count();
    if found == wanted.len() {
        // A chain letter after a name ("insulin A") does not make it another
        // molecule; a letter asked for ("vitamin C") still has to be there.
        let extra = name
            .iter()
            .filter(|n| n.chars().count() > 1 && !wanted.contains(n))
            .count();
        if extra == 0 {
            return (3.0, true);
        }
        return (2.0 - 0.08 * extra.min(10) as f32, false);
    }
    (found as f32 / wanted.len() as f32, false)
}

/// PDB titles shout. Take it down to sentence case, keeping identifiers and
/// single letters as they were.
fn tidy_title(title: &str) -> String {
    let shouted = title
        .chars()
        .filter(|c| c.is_alphabetic())
        .all(|c| c.is_uppercase());
    if !shouted {
        return title.to_string();
    }
    let mut out = String::with_capacity(title.len());
    for (i, word) in title.split_whitespace().enumerate() {
        if i > 0 {
            out.push(' ');
        }
        let keep = word.chars().any(|c| c.is_ascii_digit())
            || word.chars().filter(|c| c.is_alphabetic()).count() <= 1;
        if keep {
            out.push_str(word);
        } else {
            out.extend(word.chars().map(|c| c.to_ascii_lowercase()));
        }
    }
    let mut chars = out.chars();
    match chars.next() {
        Some(first) => first.to_uppercase().collect::<String>() + chars.as_str(),
        None => out,
    }
}

/// "HOMO SAPIENS" as "Homo sapiens".
fn sentence_organism(name: &str) -> String {
    if name.chars().filter(|c| c.is_alphabetic()).all(|c| c.is_uppercase()) {
        let lower = name.to_lowercase();
        let mut chars = lower.chars();
        match chars.next() {
            Some(first) => first.to_uppercase().collect::<String>() + chars.as_str(),
            None => lower,
        }
    } else {
        name.to_string()
    }
}

/// C8H10N4O2 as C₈H₁₀N₄O₂.
fn formula_subscripts(formula: &str) -> String {
    formula
        .chars()
        .map(|c| match c.to_digit(10) {
            Some(d) => char::from_u32(0x2080 + d).unwrap_or(c),
            None => c,
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Network
// ---------------------------------------------------------------------------

const USER_AGENT: &str = "profold/0.1 (https://github.com/odexion/profold)";

/// One client for every search, with a time limit: a list that never comes
/// is worse than one that says the search failed.
fn agent() -> &'static ureq::Agent {
    static AGENT: OnceLock<ureq::Agent> = OnceLock::new();
    AGENT.get_or_init(|| {
        ureq::Agent::config_builder()
            .timeout_global(Some(Duration::from_secs(15)))
            .build()
            .new_agent()
    })
}

fn get_json(url: &str) -> Result<Value> {
    let text = agent()
        .get(url)
        .header("user-agent", USER_AGENT)
        .call()?
        .body_mut()
        .with_config()
        .limit(8 * 1024 * 1024)
        .read_to_string()?;
    Ok(serde_json::from_str(&text)?)
}

fn post_json(url: &str, body: &Value) -> Result<Value> {
    let mut resp = agent()
        .post(url)
        .header("content-type", "application/json")
        .header("user-agent", USER_AGENT)
        .send(body.to_string())?;
    // The PDB's search answers no hits with an empty 204.
    if resp.status() == 204 {
        return Ok(Value::Null);
    }
    let text = resp.body_mut().with_config().limit(8 * 1024 * 1024).read_to_string()?;
    if text.trim().is_empty() {
        return Ok(Value::Null);
    }
    Ok(serde_json::from_str(&text)?)
}

/// Percent-encode a path segment.
pub fn encode(text: &str) -> String {
    let mut out = String::with_capacity(text.len() * 3);
    for b in text.bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_' | b'.' | b'~' => out.push(b as char),
            _ => out.push_str(&format!("%{b:02X}")),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn w(s: &str) -> Vec<String> {
        words(s)
    }

    #[test]
    fn names_match_through_the_boilerplate() {
        let q = w("green fluorescent protein");
        assert_eq!(matching(&q, &subject_words("STRUCTURE OF GREEN FLUORESCENT PROTEIN")), (3.0, true));
        assert!(!matching(&q, &subject_words("Extra-superfolder GFP")).1);
        let (contained, exact) = matching(&q, &w("Green fluorescent protein from Aequorea victoria"));
        assert!(!exact && contained > 1.0 && contained < 2.0);
        assert!(matching(&w("haemoglobin"), &w("Hemoglobin")).1);
        assert!(matching(&w("insulin"), &w("Insulin A chain")).1);
        assert!(!matching(&w("vitamin c"), &w("vitamin D binding protein")).1);
        assert_eq!(american("haemoglobin"), "hemoglobin");
    }

    #[test]
    fn exact_groups_lead() {
        let hit = |kind, exact| Hit {
            kind,
            id: String::new(),
            title: String::new(),
            detail: String::new(),
            score: 0.0,
            exact,
            alias: false,
        };
        let order = arrange(
            vec![hit(Kind::Pdb, false)],
            vec![hit(Kind::AlphaFold, false)],
            vec![hit(Kind::Molecule, true)],
        );
        assert_eq!(order[0].kind, Kind::Molecule);
        let order = arrange(
            vec![hit(Kind::Pdb, true)],
            vec![],
            vec![hit(Kind::Molecule, true)],
        );
        assert_eq!(order[0].kind, Kind::Pdb);
    }

    #[test]
    fn encodes_names_for_a_path() {
        assert_eq!(encode("vitamin c"), "vitamin%20c");
        assert_eq!(encode("2,4-D"), "2%2C4-D");
    }

    /// Live: `cargo test search_live -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn search_live() {
        for q in ["green fluorescent protein", "caffeine", "haemoglobin", "hemoglobin", "ubiquitin", "ATP", "heme", "myoglobin", "spike protein", "hiv protease", "kinase", "insulin", "aspirin", "lysozyme", "vitamin c"] {
            let hits = search(q).unwrap();
            println!("== {q}");
            for h in hits.iter() {
                println!("  {:?} {} {} | {} | {:.2} {}", h.kind, h.id, h.title, h.detail, h.score, h.exact);
            }
        }
        let hits = search("green fluorescent protein").unwrap();
        assert_eq!(hits[0].kind, Kind::Pdb);
        assert!(["1GFL", "1EMA"].contains(&hits[0].id.as_str()), "{}", hits[0].id);
        assert_eq!(search("caffeine").unwrap()[0].id, "2519");
        let crambin = search_sequence("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN").unwrap();
        assert_eq!(crambin[0].id, "1CRN");
    }
}
