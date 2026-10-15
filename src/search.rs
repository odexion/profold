//! Searching by name, so nobody has to know an identifier by heart.
//!
//! A query goes to three places at once: the PDB for experimental structures,
//! UniProt for proteins with an AlphaFold model, and PubChem for everything
//! that is not a protein. Each database ranks by its own idea of relevance,
//! and the PDB's is loose: a search for green fluorescent protein turns up red
//! ones first. So the hits are scored again here, by how closely their names
//! match what was typed, and the database whose best hit is a match for the
//! name leads.
//!
//! A gene symbol (TP53, BRCA1, EGFR) names no protein, so it is looked up in
//! UniProt's gene names as well: the protein it codes for comes up as its
//! AlphaFold model, and the PDB entries whose chains are that protein.

use std::sync::{Mutex, OnceLock};
use std::time::{Duration, Instant};

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
    /// The gene it was found by, when it was.
    gene: Option<String>,
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
            gene: None,
        }
    }

    /// Whether this is what was typed, by name or identifier, letter for
    /// letter save case and spacing. Enter before the list is back opens
    /// only such a hit; anything looser is left for the user to choose.
    pub fn is_named(&self, query: &str) -> bool {
        let squash = |t: &str| t.split_whitespace().collect::<Vec<_>>().join(" ").to_lowercase();
        let query = squash(query);
        !query.is_empty()
            && (squash(&self.title) == query
                || squash(&self.label()) == query
                || self.gene.as_deref().is_some_and(|g| squash(g) == query))
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
    // The speller is asked alongside; only a misspelling costs a second round.
    let (found, respelt) = std::thread::scope(|s| {
        let respelt = s.spawn(|| respelt(query));
        (search_all(query), respelt.join().ok().flatten())
    });
    let ((mut pdb, mut uniprot, mut molecules), gene) = found?;
    // A gene symbol is not a misspelling: the speller has BRCA1 as brca.
    if let Some(fixed) = respelt.filter(|_| !gene) {
        if let Ok(((p, u, m), _)) = search_all(&american(&fixed)) {
            pdb = merge(pdb, p, PDB_SHOWN);
            uniprot = merge(uniprot, u, UNIPROT_SHOWN);
            molecules = merge(molecules, m, MOLECULES_SHOWN);
        }
    }
    // PubChem's synonyms are generous: vitamin C is ascorbic acid, but it
    // also has ubiquitin down as citric acid. A synonym that shares no word
    // with the query is believed only when no protein answers to the name.
    let protein = |hits: &[Hit]| hits.first().is_some_and(|h| h.exact);
    if protein(&pdb) || protein(&uniprot) {
        molecules.retain(|h| !h.alias);
    }
    Ok(arrange(pdb, uniprot, molecules))
}

type Groups = (Vec<Hit>, Vec<Hit>, Vec<Hit>);

/// The three databases at once, and whether the query is surely a gene's
/// name. One that fails is left out; only when all of them fail is it an
/// error.
fn search_all(query: &str) -> Result<(Groups, bool)> {
    let (pdb, uniprot, pubchem, gene) = std::thread::scope(|s| {
        let pdb = s.spawn(|| search_pdb(query));
        let uniprot = s.spawn(|| search_uniprot(query));
        let pubchem = s.spawn(|| search_pubchem(query));
        // A gene symbol, looked up alongside; its protein's structures are
        // asked for once UniProt says which protein that is.
        let genes = genes(query);
        let gene = genes.first().is_some_and(|g| g.certain);
        let gene_pdb = genes
            .first()
            .filter(|g| g.certain)
            .cloned()
            .map(|g| s.spawn(move || search_pdb_gene(&g)));
        let mut pdb = pdb.join().unwrap_or_else(|_| Err(anyhow!("PDB search failed")));
        let mut uniprot = uniprot.join().unwrap_or_else(|_| Err(anyhow!("UniProt search failed")));
        // The text search finds whatever mentions the symbol: KRAS's
        // promoter, antibodies against EGFR. The protein's own entries are
        // what a gene symbol asks for.
        if let Some(Ok(Ok(hits))) = gene_pdb.map(|t| t.join()) {
            if !hits.is_empty() {
                pdb = Ok(hits);
            }
        }
        if !genes.is_empty() {
            let hits = genes.into_iter().map(Gene::hit).collect();
            uniprot = Ok(merge(uniprot.unwrap_or_default(), hits, UNIPROT_SHOWN));
        }
        (
            pdb,
            uniprot,
            pubchem.join().unwrap_or_else(|_| Err(anyhow!("PubChem search failed"))),
            gene,
        )
    });
    if let (Err(e), Err(_), Err(_)) = (&pdb, &uniprot, &pubchem) {
        return Err(anyhow!("search failed: {e:#}"));
    }
    Ok((
        (pdb.unwrap_or_default(), uniprot.unwrap_or_default(), pubchem.unwrap_or_default()),
        gene,
    ))
}

/// Two rankings of one database's hits as one, best first.
fn merge(mut hits: Vec<Hit>, more: Vec<Hit>, keep: usize) -> Vec<Hit> {
    for h in more {
        match hits.iter_mut().find(|k| k.id == h.id) {
            Some(k) if h.score > k.score => *k = h,
            Some(_) => {}
            None => hits.push(h),
        }
    }
    hits.sort_by(|a, b| b.score.total_cmp(&a.score));
    hits.truncate(keep);
    hits
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
    let request = |value: &str, rows: usize, sort: Value| {
        json!({
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": { "value": value }
            },
            "return_type": "entry",
            "request_options": {
                "paginate": { "start": 0, "rows": rows },
                "results_content_type": ["experimental"],
                "sort": [sort],
            }
        })
    };
    let pools = |value: &str| {
        let relevant = request(value, PDB_ROWS, json!({ "sort_by": "score", "direction": "desc" }));
        let earliest = request(
            value,
            PDB_EARLIEST,
            json!({ "sort_by": "rcsb_accession_info.initial_release_date", "direction": "asc" }),
        );
        std::thread::scope(|s| {
            let a = s.spawn(|| rcsb_ids(&relevant));
            let b = s.spawn(|| rcsb_ids(&earliest));
            (
                a.join().unwrap_or_else(|_| Err(anyhow!("PDB search failed"))),
                b.join().unwrap_or_else(|_| Err(anyhow!("PDB search failed"))),
            )
        })
    };
    // The words as typed, and again with the last left open in case it is
    // half-written; a database finds something for almost anything, so the
    // second is asked every time and the ranking below sorts the two.
    let tolerant = forgiving(query);
    let ((relevant, earliest), (more, first)) = std::thread::scope(|s| {
        let loose = s.spawn(|| pools(&tolerant));
        let strict = pools(query);
        (strict, loose.join().unwrap_or((Ok(Vec::new()), Ok(Vec::new()))))
    });
    let mut relevant = relevant?;
    let mut earliest = earliest.unwrap_or_default();
    for id in more.unwrap_or_default().into_iter().take(PDB_ROWS * 2 / 3) {
        if !relevant.contains(&id) {
            relevant.push(id);
        }
    }
    for id in first.unwrap_or_default() {
        if !earliest.contains(&id) {
            earliest.push(id);
        }
    }
    let mut ids = relevant.clone();
    ids.extend(earliest.iter().filter(|id| !relevant.contains(id)).cloned());

    let entries = describe_entries(&ids)?;
    let wanted = words(query);
    let n = relevant.len().max(1) as f32;
    let mut hits: Vec<Hit> = entries
        .into_iter()
        // An entry that names none of the words was found by a loose
        // reading of them, and has nothing to do with what was asked.
        .filter(|e| {
            matching(&wanted, &words(&e.title)).0 > 0.0
                || e.names.iter().any(|d| matching(&wanted, &words(d)).0 > 0.0)
        })
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
            let classic = if earliest.contains(&e.id) { 1.0 } else { 0.0 };
            let score = title_score + name_score + relevance + classic + e.quality
                - variation(&e.title)
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
    /// Each chain's UniProt accessions and length in residues.
    chains: Vec<(Vec<String>, u64)>,
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
            gene: None,
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
         rcsb_entity_source_organism{{scientific_name}} \
         entity_poly{{rcsb_sample_sequence_length}} \
         rcsb_polymer_entity_container_identifiers{{uniprot_ids}}}}}}}}"
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
        let chains = entities
            .iter()
            .map(|p| {
                let ids = p
                    .pointer("/rcsb_polymer_entity_container_identifiers/uniprot_ids")
                    .and_then(Value::as_array)
                    .map(|a| a.iter().filter_map(Value::as_str).map(str::to_string).collect())
                    .unwrap_or_default();
                let length = p
                    .pointer("/entity_poly/rcsb_sample_sequence_length")
                    .and_then(Value::as_u64)
                    .unwrap_or(0);
                (ids, length)
            })
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
            chains,
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
    // UniProt has no room for a slip, but the last word may be unfinished.
    let open = open_ended(query);
    let (strict, loose) = std::thread::scope(|s| {
        let loose = s.spawn(|| open.as_deref().map(|open| uniprot(open, query)));
        (uniprot(query, query), loose.join().ok().flatten())
    });
    let mut hits = match (strict, loose) {
        (Err(e), None | Some(Err(_))) => return Err(e),
        (strict, loose) => {
            let mut hits = strict.unwrap_or_default();
            for h in loose.and_then(Result::ok).unwrap_or_default() {
                if !hits.iter().any(|k| k.id == h.id) {
                    hits.push(h);
                }
            }
            hits
        }
    };
    hits.sort_by(|a, b| b.score.total_cmp(&a.score));
    hits.truncate(UNIPROT_SHOWN);
    Ok(hits)
}

/// UniProt's reviewed entries for `asked`, ranked against what was typed.
fn uniprot(asked: &str, query: &str) -> Result<Vec<Hit>> {
    let reply = agent()
        .get("https://rest.uniprot.org/uniprotkb/search")
        .query("query", format!("({asked}) AND reviewed:true"))
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

    let hits: Vec<Hit> = json
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
            // Found for a word in its notes, not its name.
            if score == 0.0 {
                return None;
            }
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
                gene: None,
            })
        })
        .collect();
    Ok(hits)
}

// ---------------------------------------------------------------------------
// Gene symbols
// ---------------------------------------------------------------------------

/// A reviewed UniProt entry found by the name of its gene.
#[derive(Clone, Debug)]
struct Gene {
    accession: String,
    /// The gene's own name, as UniProt writes it: TP53, Tp53 or tp53.
    symbol: String,
    protein: String,
    organism: String,
    length: Option<u64>,
    score: f32,
    /// Surely the gene that was meant: its name, not a synonym, and either
    /// human or written the way it was typed. Lower-case heme is a molecule,
    /// not hemE of E. coli.
    certain: bool,
}

impl Gene {
    fn hit(self) -> Hit {
        let mut detail = vec![self.symbol.clone(), self.organism];
        detail.extend(self.length.map(|n| format!("{n} aa")));
        Hit {
            kind: Kind::AlphaFold,
            id: self.accession,
            title: self.protein,
            detail: detail.join(" · "),
            score: self.score,
            exact: self.certain,
            alias: false,
            gene: Some(self.symbol),
        }
    }
}

/// Whether a query could be a gene symbol: one short word of letters and
/// digits, with perhaps a hyphen (HLA-A).
fn gene_like(query: &str) -> bool {
    let n = query.chars().count();
    (2..=15).contains(&n)
        && query.chars().all(|c| c.is_ascii_alphanumeric() || c == '-')
        && query.chars().any(|c| c.is_ascii_alphabetic())
}

/// The proteins a gene symbol codes for, best first: those whose gene it is
/// by name before those that have it as a synonym, and the human one first
/// of each, as most who type TP53 mean ours. Empty when the query is not a
/// gene, or UniProt cannot be reached.
fn genes(query: &str) -> Vec<Gene> {
    if !gene_like(query) {
        return Vec::new();
    }
    let reply = agent()
        .get("https://rest.uniprot.org/uniprotkb/search")
        .query("query", format!("gene_exact:\"{query}\" AND reviewed:true"))
        .query("fields", "accession,protein_name,organism_name,gene_names,length")
        .query("size", "8")
        .query("format", "json")
        .header("accept-encoding", "identity")
        .header("user-agent", USER_AGENT)
        .call()
        .and_then(|mut r| r.body_mut().with_config().limit(8 * 1024 * 1024).read_to_string());
    let Some(json) = reply.ok().and_then(|t| serde_json::from_str::<Value>(&t).ok()) else {
        return Vec::new();
    };
    let results = json.get("results").and_then(Value::as_array).cloned().unwrap_or_default();
    let n = results.len().max(1) as f32;
    let mut genes: Vec<Gene> = results
        .iter()
        .enumerate()
        .filter_map(|(rank, r)| {
            let accession = r.get("primaryAccession")?.as_str()?.to_string();
            let protein = r
                .pointer("/proteinDescription/recommendedName/fullName/value")
                .or_else(|| r.pointer("/proteinDescription/submissionNames/0/fullName/value"))
                .and_then(Value::as_str)?
                .to_string();
            let named = |v: &Value| v.as_str().is_some_and(|g| g.eq_ignore_ascii_case(query));
            let entries = r.get("genes").and_then(Value::as_array).cloned().unwrap_or_default();
            let primary = entries.iter().find_map(|g| {
                let name = g.pointer("/geneName/value")?;
                named(name).then(|| name.as_str().unwrap_or(query).to_string())
            });
            // gene_exact also finds synonyms, ordered locus names and the like.
            let synonym = entries.iter().any(|g| {
                ["synonyms", "orderedLocusNames", "orfNames"].iter().any(|k| {
                    g.get(*k)
                        .and_then(Value::as_array)
                        .is_some_and(|a| a.iter().any(|s| s.get("value").is_some_and(named)))
                })
            });
            if primary.is_none() && !synonym {
                return None;
            }
            let organism = r
                .pointer("/organism/scientificName")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_string();
            let human = organism == "Homo sapiens";
            let certain = primary.as_deref().is_some_and(|p| human || p == query);
            Some(Gene {
                accession,
                symbol: primary.unwrap_or_else(|| query.to_uppercase()),
                protein,
                length: r.pointer("/sequence/length").and_then(Value::as_u64),
                // An uncertain one is listed, but below proteins named so.
                score: if certain { 4.0 } else { 1.2 }
                    + if human { 0.5 } else { 0.0 }
                    + 0.3 * (1.0 - rank as f32 / n),
                certain,
                organism,
            })
        })
        .collect();
    genes.sort_by(|a, b| b.score.total_cmp(&a.score));
    genes
}

/// The PDB's experimental structures of the protein a gene codes for: the
/// earliest, which are the classics, and the best resolved. Those with more
/// of the protein in them rank higher, so a domain beats a peptide of it,
/// and the protein alone beats it caught with partners.
fn search_pdb_gene(gene: &Gene) -> Result<Vec<Hit>> {
    let request = |rows: usize, sort: &str| {
        json!({
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                            "operator": "exact_match",
                            "value": gene.accession,
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                            "operator": "exact_match",
                            "value": "UniProt",
                        }
                    }
                ]
            },
            "return_type": "entry",
            "request_options": {
                "paginate": { "start": 0, "rows": rows },
                "results_content_type": ["experimental"],
                "sort": [{ "sort_by": sort, "direction": "asc" }],
            }
        })
    };
    let earliest = request(PDB_EARLIEST, "rcsb_accession_info.initial_release_date");
    let sharpest = request(PDB_ROWS, "rcsb_entry_info.resolution_combined");
    let (earliest, sharpest) = std::thread::scope(|s| {
        let b = s.spawn(|| rcsb_ids(&sharpest));
        (rcsb_ids(&earliest), b.join().unwrap_or_else(|_| Err(anyhow!("PDB search failed"))))
    });
    let earliest = earliest?;
    let mut ids = earliest.clone();
    ids.extend(sharpest.unwrap_or_default().into_iter().filter(|id| !earliest.contains(id)));

    let protein = words(&gene.protein);
    let mut hits: Vec<Hit> = describe_entries(&ids)?
        .into_iter()
        .map(|e| {
            let named = e
                .names
                .iter()
                .map(|d| matching(&protein, &words(d)).0)
                .fold(0.0f32, f32::max);
            let classic = if earliest.contains(&e.id) { 0.7 } else { 0.0 };
            let covered = e
                .chains
                .iter()
                .filter(|(ids, _)| ids.contains(&gene.accession))
                .map(|(_, n)| *n)
                .max()
                .unwrap_or(0);
            let coverage = gene.length.map_or(0.0, |l| (covered as f32 / l.max(1) as f32).min(1.0));
            // Every other chain is a partner the protein was caught with.
            let partners = 0.2 * e.names.len().saturating_sub(1).min(3) as f32;
            let score = 2.0 + 0.3 * named + 1.5 * coverage + classic + e.quality
                - variation(&e.title)
                - partners;
            e.hit(score, true)
        })
        .collect();
    hits.sort_by(|a, b| b.score.total_cmp(&a.score));
    hits.truncate(PDB_SHOWN);
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
            let esearch = |term: &str| -> Result<Vec<String>> {
                ncbi_turn();
                let reply = agent()
                    .get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi")
                    .query("db", "pccompound")
                    .query("term", term)
                    .query("retmode", "json")
                    .query("retmax", "12")
                    .query("sort", "relevance")
                    .header("user-agent", USER_AGENT)
                    .call()?
                    .body_mut()
                    .read_to_string()?;
                let json: Value = serde_json::from_str(&reply)?;
                Ok(json
                    .pointer("/esearchresult/idlist")
                    .and_then(Value::as_array)
                    .map(|a| a.iter().filter_map(Value::as_str).map(str::to_string).collect())
                    .unwrap_or_default())
            };
            // As typed, and with the last word left open in case it is
            // unfinished: "caffe" finds caffeine only the second way.
            let open = open_ended(query);
            let (strict, loose) = std::thread::scope(|s| {
                let loose = s.spawn(|| open.as_deref().map(esearch));
                (esearch(query), loose.join().ok().flatten())
            });
            let ids = match (strict, loose) {
                (Err(e), None | Some(Err(_))) => return Err(e),
                (strict, loose) => {
                    let mut ids = strict.unwrap_or_default();
                    for id in loose.and_then(Result::ok).unwrap_or_default() {
                        if !ids.contains(&id) {
                            ids.push(id);
                        }
                    }
                    ids
                }
            };
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
        .filter_map(|(rank, mut h)| {
            let (score, exact) = matching(&wanted, &words(&h.title));
            // Something whose name lacks a word asked for was found by a
            // synonym, a patent or that word alone (Fluo-4 for "green
            // fluo"), and is noise in a list like this.
            if score <= 1.0 {
                return None;
            }
            h.exact = exact;
            h.score = score + 0.6 * (1.0 - rank as f32 / n);
            Some(h)
        })
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
                gene: None,
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
const STOP: [&str; 33] = [
    "a", "an", "the", "of", "in", "from", "with", "and", "for", "by", "to", "at",
    "structure", "structures", "crystal", "crystallographic", "solution", "nmr",
    "x", "ray", "chain", "subunit", "refinement", "refined", "resolution",
    "angstrom", "angstroms", "analysis", "determination", "three", "dimensional",
    "high", "form",
];

/// Words that name a chain or a state of the molecule rather than another
/// molecule: haemoglobin's alpha chain or its deoxy form is still haemoglobin.
const NEUTRAL: [&str; 10] = [
    "alpha", "beta", "gamma", "delta", "deoxy", "oxy", "apo", "holo", "native", "wild",
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

/// How much a title says the entry is a variation on the molecule rather
/// than the molecule: a mutant, a complex, a fragment. Asked for by name,
/// the plain one is wanted.
fn variation(title: &str) -> f32 {
    const MARKS: [&str; 12] = [
        "mutant", "variant", "complex", "complexed", "inhibitor", "bound",
        "fragment", "domain", "analog", "analogue", "fusion", "engineered",
    ];
    let lower = title.to_lowercase();
    let found = lower
        .split(|c: char| !c.is_alphanumeric())
        .filter(|w| MARKS.contains(w))
        .count();
    0.35 * found.min(2) as f32
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
/// per extra word when it has all of it, and the share of the query's words
/// it has otherwise. A word typed half-way or with a slip counts for a little
/// less than the word itself, and a name with every word answered never
/// scores below 1.1. The flag says whether it was the query exactly.
fn matching(wanted: &[String], name: &[String]) -> (f32, bool) {
    if wanted.is_empty() || name.is_empty() {
        return (0.0, false);
    }
    let best = |w: &String| name.iter().map(|n| closeness(w, n)).fold(0.0f32, f32::max);
    let found: Vec<f32> = wanted.iter().map(best).collect();
    if found.iter().any(|c| *c == 0.0) {
        return (found.iter().sum::<f32>() / wanted.len() as f32, false);
    }
    // A chain letter after a name ("insulin A") does not make it another
    // molecule; a letter asked for ("vitamin C") still has to be there.
    let extra = name
        .iter()
        .filter(|n| {
            n.chars().count() > 1
                && !NEUTRAL.contains(&n.as_str())
                && !wanted.iter().any(|w| closeness(w, n) > 0.0)
        })
        .count();
    let slack: f32 = found.iter().map(|c| 1.0 - c).sum();
    let base = if extra == 0 { 3.0 } else { 2.0 - 0.08 * extra.min(10) as f32 };
    ((base - slack).max(1.1), extra == 0 && slack == 0.0)
}

/// How far a name's word answers one word of the query: 1 when it is that
/// word, and less when it holds it (hemoglobin in deoxyhemoglobin), begins
/// with it (fluo, still being typed), or is it with a slip of the keys
/// (flourescent). Nothing otherwise.
fn closeness(w: &str, n: &str) -> f32 {
    if n == w {
        return 1.0;
    }
    let len = w.chars().count();
    if len > 5 && n.contains(w) {
        return 0.9;
    }
    if len >= 3 && n.starts_with(w) {
        return 0.85;
    }
    // One slip in a word of four letters or more, two from eight.
    let slips = match len {
        0..=3 => return 0.0,
        4..=7 => 1,
        _ => 2,
    };
    if slips_between(w, n) <= slips {
        return 0.7;
    }
    // A slip in a word still being typed: flou for fluorescent.
    let head: String = n.chars().take(len).collect();
    if n.chars().count() > len && slips_between(w, &head) <= 1 {
        return 0.6;
    }
    0.0
}

/// Letters typed wrong, missed, added, or swapped with the next one.
fn slips_between(a: &str, b: &str) -> usize {
    let a: Vec<char> = a.chars().collect();
    let b: Vec<char> = b.chars().collect();
    let mut d = vec![vec![0usize; b.len() + 1]; a.len() + 1];
    for (i, row) in d.iter_mut().enumerate() {
        row[0] = i;
    }
    for j in 0..=b.len() {
        d[0][j] = j;
    }
    for i in 1..=a.len() {
        for j in 1..=b.len() {
            let cost = usize::from(a[i - 1] != b[j - 1]);
            d[i][j] = (d[i - 1][j] + 1).min(d[i][j - 1] + 1).min(d[i - 1][j - 1] + cost);
            if i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1] {
                d[i][j] = d[i][j].min(d[i - 2][j - 2] + 1);
            }
        }
    }
    d[a.len()][b.len()]
}

/// The query with its last word left open, for a word still being typed:
/// "green fluo" as "green fluo*".
fn open_ended(query: &str) -> Option<String> {
    let last = query.split_whitespace().last()?;
    (last.chars().count() >= 3 && last.chars().all(char::is_alphanumeric))
        .then(|| format!("{query}*"))
}

/// The query as the PDB's text search takes it with its last word left
/// open, as typed and with each pair of its letters swapped: "green flou"
/// finds green fluorescent protein by way of fluo*. Misspelt whole words
/// are for [`respelt`].
fn forgiving(query: &str) -> String {
    let mut words: Vec<String> = query.split_whitespace().map(str::to_string).collect();
    if let Some(last) = words.last_mut() {
        let chars: Vec<char> = last.chars().collect();
        if chars.len() >= 3 && chars.iter().all(|c| c.is_alphanumeric()) {
            let mut options = vec![format!("{last}*")];
            for k in 1..chars.len().min(10) {
                let mut swapped = chars.clone();
                swapped.swap(k - 1, k);
                let option = format!("{}*", swapped.into_iter().collect::<String>());
                if !options.contains(&option) {
                    options.push(option);
                }
            }
            *last = format!("({})", options.join(" | "));
        }
    }
    words.join(" ")
}

/// The query with its misspellings put right, by NCBI's speller for the
/// biomedical literature: "green flourescent protein", "hemoglobn",
/// "asprin". None when it has nothing to correct, or would change how many
/// words there are ("green f lou" for a word still being typed).
fn respelt(query: &str) -> Option<String> {
    ncbi_turn();
    let reply = agent()
        .get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/espell.fcgi")
        .query("db", "pubmed")
        .query("term", query)
        .header("user-agent", USER_AGENT)
        .call()
        .ok()?
        .body_mut()
        .read_to_string()
        .ok()?;
    let corrected = reply
        .split_once("<CorrectedQuery>")?
        .1
        .split_once("</CorrectedQuery>")?
        .0
        .trim()
        .to_string();
    let count = |t: &str| t.split_whitespace().count();
    (!corrected.is_empty()
        && corrected.to_lowercase() != query.to_lowercase()
        && count(&corrected) == count(query))
    .then_some(corrected)
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

/// Wait for a turn at NCBI's E-utilities, which turn away more than three
/// requests a second from one address; a search asks it up to five times.
fn ncbi_turn() {
    static LAST: Mutex<Option<Instant>> = Mutex::new(None);
    let mut last = LAST.lock().unwrap_or_else(|e| e.into_inner());
    if let Some(at) = *last {
        let next = at + Duration::from_millis(350);
        let now = Instant::now();
        if next > now {
            std::thread::sleep(next - now);
        }
    }
    *last = Some(Instant::now());
}

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
    fn half_typed_and_misspelt_words_still_match() {
        let gfp = w("Green fluorescent protein");
        // Still typing, or a slip of the keys: every word is answered.
        assert!(matching(&w("green fluo"), &gfp).0 > 1.0);
        assert!(matching(&w("green flou"), &gfp).0 > 1.0);
        assert!(matching(&w("green flourescent protein"), &gfp).0 > 1.0);
        assert!(matching(&w("hemoglobn"), &w("Hemoglobin")).0 > 1.0);
        // ...but never as well as the words themselves, nor as exactly.
        let (typo, exact) = matching(&w("green flourescent protein"), &gfp);
        assert!(!exact && typo < matching(&w("green fluorescent protein"), &gfp).0);
        // A name lacking a word asked for is not a match: Fluo-4 has no green.
        assert!(matching(&w("green fluo"), &w("Fluo-4")).0 <= 1.0);
        assert!(matching(&w("green fluo"), &w("Fluo-3")).0 <= 1.0);
        // Short words are not stretched to fit.
        assert_eq!(closeness("cat", "car"), 0.0);
        assert_eq!(slips_between("flourescent", "fluorescent"), 1);
        assert_eq!(forgiving("green flou"), "green (flou* | lfou* | folu* | fluo*)");
        assert_eq!(open_ended("green fluo").as_deref(), Some("green fluo*"));
        assert_eq!(open_ended("vitamin c"), None);
    }

    #[test]
    fn only_the_name_itself_is_named() {
        let caffeine = Hit::for_test(Kind::Molecule, "2519", "Caffeine");
        assert!(caffeine.is_named("caffeine"));
        assert!(caffeine.is_named("  CAFFEINE "));
        assert!(caffeine.is_named("cid 2519"));
        assert!(!caffeine.is_named("caffein"));
        assert!(!caffeine.is_named("caffeine citrate"));
        let gfp = Hit::for_test(Kind::Pdb, "1EMA", "Green fluorescent protein");
        assert!(gfp.is_named("green  fluorescent protein"));
        assert!(!gfp.is_named("green protein"));
        assert!(!gfp.is_named("asdkjh qwe"));
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
            gene: None,
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
    fn gene_symbols_look_like_genes() {
        for g in ["TP53", "brca1", "HLA-A", "Kras"] {
            assert!(gene_like(g), "{g}");
        }
        for q in ["green fluorescent protein", "vitamin c", "53", "2,4-D", "a"] {
            assert!(!gene_like(q), "{q}");
        }
        let mut hit = Hit::for_test(Kind::AlphaFold, "P04637", "Cellular tumor antigen p53");
        hit.gene = Some("TP53".into());
        assert!(hit.is_named("tp53"));
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
        for q in ["green fluorescent protein", "caffeine", "haemoglobin", "hemoglobin", "ubiquitin", "ATP", "heme", "myoglobin", "spike protein", "hiv protease", "kinase", "insulin", "aspirin", "lysozyme", "vitamin c", "green fluo", "green flou", "green flourescent protein", "caffe", "hemoglobn", "ubiquit", "asprin", "TP53", "brca1", "EGFR", "kras", "GFP", "gfp"] {
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
        // Half-typed and misspelt: the protein, and no molecule lacking a word.
        for q in ["green fluo", "green flou", "green flourescent protein"] {
            let hits = search(q).unwrap();
            assert!(hits.iter().any(|h| h.title.to_lowercase().contains("green fluorescent")), "{q}");
            assert!(hits.iter().all(|h| h.kind != Kind::Molecule || h.title.to_lowercase().contains("green")), "{q}");
        }
        // A gene symbol: the human protein, as a model and as structures.
        let p53 = search("TP53").unwrap();
        assert!(p53.iter().any(|h| h.kind == Kind::AlphaFold && h.id == "P04637"));
        assert!(p53.iter().any(|h| h.kind == Kind::Pdb));
        assert!(p53[0].is_named("tp53") || p53[0].kind == Kind::Pdb);
        // Lower-case heme is the molecule, whatever E. coli calls a gene.
        assert_eq!(search("heme").unwrap()[0].kind, Kind::Molecule);
        // Not respelt as brca, whose hits are BRCA2's.
        let brca1 = search("brca1").unwrap();
        assert_eq!(brca1[0].kind, Kind::Pdb);
        assert!(brca1.iter().any(|h| h.id == "P38398"));
        let crambin = search_sequence("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN").unwrap();
        assert_eq!(crambin[0].id, "1CRN");
    }
}
