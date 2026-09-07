//! Structure model plus PDB and mmCIF parsers.
//!
//! Everything downstream (renderer, sequence view, colouring) reads this and
//! nothing else, so the two parsers are the only places that know about file
//! formats.

use std::ops::Range;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Element {
    C,
    N,
    O,
    S,
    P,
    H,
    Other,
}

impl Element {
    fn from_symbol(sym: &str) -> Self {
        match sym.trim().to_ascii_uppercase().as_str() {
            "C" => Element::C,
            "N" => Element::N,
            "O" => Element::O,
            "S" => Element::S,
            "P" => Element::P,
            "H" | "D" => Element::H,
            _ => Element::Other,
        }
    }

    /// Guess from a PDB atom name when the element column is missing.
    fn from_atom_name(name: &str) -> Self {
        let n = name.trim_start_matches(|c: char| c.is_ascii_digit());
        match n.chars().next().map(|c| c.to_ascii_uppercase()) {
            Some('C') => Element::C,
            Some('N') => Element::N,
            Some('O') => Element::O,
            Some('S') => Element::S,
            Some('P') => Element::P,
            Some('H') => Element::H,
            _ => Element::Other,
        }
    }

    /// Van der Waals radius in ångström.
    pub fn vdw(self) -> f32 {
        match self {
            Element::C => 1.70,
            Element::N => 1.55,
            Element::O => 1.52,
            Element::S => 1.80,
            Element::P => 1.80,
            Element::H => 1.10,
            Element::Other => 1.60,
        }
    }

    /// CPK-ish colours, desaturated a little so they sit well on paper white.
    pub fn color(self) -> [f32; 3] {
        match self {
            Element::C => [0.42, 0.44, 0.47],
            Element::N => [0.24, 0.42, 0.78],
            Element::O => [0.83, 0.28, 0.26],
            Element::S => [0.85, 0.68, 0.20],
            Element::P => [0.90, 0.51, 0.20],
            Element::H => [0.85, 0.85, 0.87],
            Element::Other => [0.55, 0.45, 0.70],
        }
    }
}

#[derive(Clone, Copy)]
pub struct Atom {
    pub pos: [f32; 3],
    pub element: Element,
    /// B-factor, or pLDDT for predicted models.
    pub b: f32,
    pub residue: u32,
}

/// Secondary structure, assigned from backbone geometry.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Ss {
    Coil,
    Helix,
    Sheet,
}

#[derive(Clone)]
pub struct Residue {
    pub name: String,
    /// Author numbering, which often does not start at one.
    pub seq: i32,
    pub chain: u16,
    /// Indices into `Structure::atoms` of the backbone atoms.
    pub ca: Option<u32>,
    pub n: Option<u32>,
    pub c: Option<u32>,
    pub o: Option<u32>,
    pub b: f32,
    pub code: char,
    pub ss: Ss,
}

impl Residue {
    pub fn backbone(&self) -> Option<(u32, u32, u32, u32)> {
        Some((self.n?, self.ca?, self.c?, self.o?))
    }
}

#[derive(Clone)]
pub struct Chain {
    pub id: String,
    /// Indices into `Structure::residues`, in file order.
    pub residues: Vec<u32>,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Source {
    Pdb,
    AlphaFold,
    SimpleFold,
}

impl Source {
    pub fn label(self) -> &'static str {
        match self {
            Source::Pdb => "RCSB PDB",
            Source::AlphaFold => "AlphaFold DB",
            Source::SimpleFold => "SimpleFold",
        }
    }
}

#[derive(Clone)]
pub struct Structure {
    pub id: String,
    pub title: String,
    pub source: Source,
    pub atoms: Vec<Atom>,
    pub residues: Vec<Residue>,
    pub chains: Vec<Chain>,
    /// True when the B-factor column holds pLDDT rather than a B-factor.
    pub has_confidence: bool,
    /// Per-residue deviation from a superimposed structure, in angstrom.
    /// Negative where there is no counterpart. Empty unless comparing.
    pub deviation: Vec<f32>,
    /// A flat colour, set when this structure is one half of a comparison.
    pub tint: Option<[f32; 3]>,
}

impl Structure {
    pub fn sequence(&self) -> String {
        self.residues.iter().map(|r| r.code).collect()
    }

    /// The longest protein chain: its residue indices and one-letter sequence.
    /// This is what gets folded when comparing against a prediction, so the
    /// two structures line up residue for residue afterwards.
    pub fn primary_chain(&self) -> (Vec<u32>, String) {
        let best = self
            .chains
            .iter()
            .map(|c| {
                c.residues
                    .iter()
                    .filter(|&&ri| {
                        let r = &self.residues[ri as usize];
                        r.ca.is_some() && r.code != 'X'
                    })
                    .copied()
                    .collect::<Vec<u32>>()
            })
            .max_by_key(|v| v.len())
            .unwrap_or_default();
        let seq = best
            .iter()
            .map(|&ri| self.residues[ri as usize].code)
            .collect();
        (best, seq)
    }

    pub fn mean_confidence(&self) -> Option<f32> {
        if !self.has_confidence || self.residues.is_empty() {
            return None;
        }
        let sum: f32 = self.residues.iter().map(|r| r.b).sum();
        Some(sum / self.residues.len() as f32)
    }

}

// ---------------------------------------------------------------------------
// Shared residue handling
// ---------------------------------------------------------------------------

pub fn three_to_one(name: &str) -> char {
    match name.trim().to_ascii_uppercase().as_str() {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" => 'C',
        "GLN" => 'Q',
        "GLU" => 'E',
        "GLY" => 'G',
        "HIS" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" => 'K',
        "MET" => 'M',
        "PHE" => 'F',
        "PRO" => 'P',
        "SER" => 'S',
        "THR" => 'T',
        "TRP" => 'W',
        "TYR" => 'Y',
        "VAL" => 'V',
        "SEC" => 'U',
        "PYL" => 'O',
        "MSE" => 'M',
        _ => 'X',
    }
}

/// Accumulates atoms into residues and chains regardless of source format.
struct Builder {
    atoms: Vec<Atom>,
    residues: Vec<Residue>,
    chains: Vec<Chain>,
    cur_res_key: Option<(usize, String, i32, String)>,
    b_sum: f32,
    b_n: f32,
}

impl Builder {
    fn new() -> Self {
        Self {
            atoms: Vec::new(),
            residues: Vec::new(),
            chains: Vec::new(),
            cur_res_key: None,
            b_sum: 0.0,
            b_n: 0.0,
        }
    }

    fn flush_residue(&mut self) {
        if let Some(r) = self.residues.last_mut() {
            if self.b_n > 0.0 {
                r.b = self.b_sum / self.b_n;
            }
        }
        self.b_sum = 0.0;
        self.b_n = 0.0;
    }

    #[allow(clippy::too_many_arguments)]
    fn push(
        &mut self,
        chain_id: &str,
        res_name: &str,
        res_seq: i32,
        ins_code: &str,
        atom_name: &str,
        element: Element,
        pos: [f32; 3],
        b: f32,
    ) {
        // Ligand records revisit chains that already appeared, so look the
        // chain up by name rather than assuming one contiguous block per chain.
        let chain = match self.chains.iter().position(|c| c.id == chain_id) {
            Some(i) => i,
            None => {
                self.chains.push(Chain {
                    id: chain_id.to_string(),
                    residues: Vec::new(),
                });
                self.chains.len() - 1
            }
        };

        let key = (chain, res_name.to_string(), res_seq, ins_code.to_string());
        if self.cur_res_key.as_ref() != Some(&key) {
            self.flush_residue();
            self.cur_res_key = Some(key);
            self.chains[chain].residues.push(self.residues.len() as u32);
            self.residues.push(Residue {
                name: res_name.to_string(),
                seq: res_seq,
                chain: chain as u16,
                ca: None,
                n: None,
                c: None,
                o: None,
                b,
                code: three_to_one(res_name),
                ss: Ss::Coil,
            });
        }

        let residue = (self.residues.len() - 1) as u32;
        let idx = self.atoms.len() as u32;
        {
            let r = &mut self.residues[residue as usize];
            match (atom_name.trim(), element) {
                ("CA", Element::C) => r.ca = Some(idx),
                ("N", Element::N) => r.n = Some(idx),
                ("C", Element::C) => r.c = Some(idx),
                ("O", Element::O) => r.o = Some(idx),
                _ => {}
            }
        }
        self.b_sum += b;
        self.b_n += 1.0;
        self.atoms.push(Atom {
            pos,
            element,
            b,
            residue,
        });
    }

    fn finish(mut self, id: String, title: String, source: Source) -> Structure {
        self.flush_residue();
        self.chains.retain(|c| !c.residues.is_empty());

        // A predicted model stores pLDDT in the B-factor column, in 0..=100.
        // It has to actually vary: SimpleFold writes a flat 100.00 for every
        // atom when confidence scoring is switched off, which is not a score.
        let (lo, hi) = self
            .atoms
            .iter()
            .fold((f32::MAX, f32::MIN), |(lo, hi), a| (lo.min(a.b), hi.max(a.b)));
        let has_confidence = matches!(source, Source::AlphaFold | Source::SimpleFold)
            && lo >= 0.0
            && hi <= 100.0
            && hi - lo > 0.5;

        let mut structure = Structure {
            id,
            title,
            source,
            atoms: self.atoms,
            residues: self.residues,
            chains: self.chains,
            has_confidence,
            deviation: Vec::new(),
            tint: None,
        };
        crate::dssp::assign(&mut structure);
        structure
    }
}

// ---------------------------------------------------------------------------
// PDB
// ---------------------------------------------------------------------------

pub fn parse_pdb(text: &str, id: String, source: Source) -> anyhow::Result<Structure> {
    let mut b = Builder::new();
    let mut title = String::new();

    for line in text.lines() {
        if line.starts_with("TITLE ") {
            let part = line.get(10..).unwrap_or("").trim();
            if !title.is_empty() {
                title.push(' ');
            }
            title.push_str(part);
            continue;
        }
        // Only the first model of an ensemble.
        if line.starts_with("ENDMDL") {
            break;
        }
        if !(line.starts_with("ATOM  ") || line.starts_with("HETATM")) {
            continue;
        }
        let col = |r: Range<usize>| line.get(r).unwrap_or("").trim().to_string();

        let res_name = col(17..20);
        if res_name == "HOH" || res_name == "DOD" {
            continue;
        }
        // Keep only the first alternate location.
        let alt = line.as_bytes().get(16).copied().unwrap_or(b' ');
        if alt != b' ' && alt != b'A' {
            continue;
        }

        let atom_name = col(12..16);
        let chain_id = {
            let c = col(21..22);
            if c.is_empty() {
                "A".to_string()
            } else {
                c
            }
        };
        let res_seq: i32 = col(22..26).parse().unwrap_or(0);
        let ins_code = col(26..27);
        let (Ok(x), Ok(y), Ok(z)) = (
            col(30..38).parse::<f32>(),
            col(38..46).parse::<f32>(),
            col(46..54).parse::<f32>(),
        ) else {
            continue;
        };
        let bf: f32 = col(60..66).parse().unwrap_or(0.0);
        let element = {
            let sym = col(76..78);
            if sym.is_empty() {
                Element::from_atom_name(&atom_name)
            } else {
                Element::from_symbol(&sym)
            }
        };

        b.push(
            &chain_id, &res_name, res_seq, &ins_code, &atom_name, element, [x, y, z], bf,
        );
    }

    let s = b.finish(id, title, source);
    if s.atoms.is_empty() {
        anyhow::bail!("no atoms found in PDB file");
    }
    Ok(s)
}

// ---------------------------------------------------------------------------
// mmCIF
// ---------------------------------------------------------------------------

fn unquote(v: &str) -> String {
    let v = v.trim();
    if v.len() >= 2
        && ((v.starts_with('\'') && v.ends_with('\'')) || (v.starts_with('"') && v.ends_with('"')))
    {
        v[1..v.len() - 1].to_string()
    } else {
        v.to_string()
    }
}

/// Split an mmCIF data line into values, honouring single and double quotes.
fn cif_tokens(line: &str) -> Vec<&str> {
    let bytes = line.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        while i < bytes.len() && (bytes[i] == b' ' || bytes[i] == b'\t') {
            i += 1;
        }
        if i >= bytes.len() {
            break;
        }
        let quote = bytes[i];
        if quote == b'\'' || quote == b'"' {
            i += 1;
            let start = i;
            while i < bytes.len() && bytes[i] != quote {
                i += 1;
            }
            out.push(&line[start..i]);
            i += 1;
        } else {
            let start = i;
            while i < bytes.len() && bytes[i] != b' ' && bytes[i] != b'\t' {
                i += 1;
            }
            out.push(&line[start..i]);
        }
    }
    out
}

pub fn parse_mmcif(text: &str, id: String, source: Source) -> anyhow::Result<Structure> {
    let mut b = Builder::new();
    let mut title = String::new();

    let mut lines = text.lines().peekable();
    while let Some(line) = lines.next() {
        let trimmed = line.trim();

        if let Some(rest) = trimmed.strip_prefix("_struct.title") {
            if !rest.is_empty() && !rest.starts_with(char::is_whitespace) {
                continue;
            }
            let rest = rest.trim();
            if !rest.is_empty() {
                title = unquote(rest);
            } else if let Some(next) = lines.next() {
                // The value sits on the following line, either quoted or as a
                // semicolon-delimited block.
                let next = next.trim_end();
                if let Some(first) = next.strip_prefix(';') {
                    title = first.trim().to_string();
                    for line in lines.by_ref() {
                        let l = line.trim_end();
                        if l.trim() == ";" {
                            break;
                        }
                        if !title.is_empty() {
                            title.push(' ');
                        }
                        title.push_str(l.trim());
                    }
                } else {
                    title = unquote(next.trim());
                }
            }
            continue;
        }

        if trimmed != "loop_" {
            continue;
        }

        // Collect this loop's column names.
        let mut cols: Vec<String> = Vec::new();
        while let Some(peek) = lines.peek() {
            let p = peek.trim();
            if p.starts_with('_') {
                cols.push(p.split_whitespace().next().unwrap_or("").to_string());
                lines.next();
            } else {
                break;
            }
        }
        if !cols.first().is_some_and(|c| c.starts_with("_atom_site.")) {
            continue;
        }

        let idx = |name: &str| cols.iter().position(|c| c == name);
        let i_group = idx("_atom_site.group_PDB");
        let i_symbol = idx("_atom_site.type_symbol");
        let i_atom = idx("_atom_site.label_atom_id");
        let i_alt = idx("_atom_site.label_alt_id");
        let i_comp = idx("_atom_site.auth_comp_id").or_else(|| idx("_atom_site.label_comp_id"));
        let i_asym = idx("_atom_site.auth_asym_id").or_else(|| idx("_atom_site.label_asym_id"));
        let i_seq = idx("_atom_site.auth_seq_id").or_else(|| idx("_atom_site.label_seq_id"));
        let i_ins = idx("_atom_site.pdbx_PDB_ins_code");
        let i_x = idx("_atom_site.Cartn_x");
        let i_y = idx("_atom_site.Cartn_y");
        let i_z = idx("_atom_site.Cartn_z");
        let i_b = idx("_atom_site.B_iso_or_equiv");
        let i_model = idx("_atom_site.pdbx_PDB_model_num");

        let (Some(i_x), Some(i_y), Some(i_z), Some(i_comp), Some(i_asym), Some(i_seq)) =
            (i_x, i_y, i_z, i_comp, i_asym, i_seq)
        else {
            continue;
        };

        let mut first_model: Option<String> = None;
        for row in lines.by_ref() {
            let t = row.trim();
            if t.is_empty() || t == "#" || t.starts_with('_') || t == "loop_" {
                break;
            }
            let f = cif_tokens(row);
            if f.len() < cols.len() {
                continue;
            }
            let get = |i: Option<usize>| i.and_then(|i| f.get(i).copied()).unwrap_or("");

            if let Some(model) = i_model.and_then(|i| f.get(i).copied()) {
                match &first_model {
                    None => first_model = Some(model.to_string()),
                    Some(m) if m != model => break,
                    _ => {}
                }
            }
            let group = get(i_group);
            let res_name = f[i_comp];
            if res_name == "HOH" || res_name == "DOD" || group == "HETATM" && res_name == "HOH" {
                continue;
            }
            let alt = get(i_alt);
            if !(alt.is_empty() || alt == "." || alt == "?" || alt == "A") {
                continue;
            }

            let (Ok(x), Ok(y), Ok(z)) = (
                f[i_x].parse::<f32>(),
                f[i_y].parse::<f32>(),
                f[i_z].parse::<f32>(),
            ) else {
                continue;
            };
            let atom_name = get(i_atom);
            let element = {
                let sym = get(i_symbol);
                if sym.is_empty() || sym == "." {
                    Element::from_atom_name(atom_name)
                } else {
                    Element::from_symbol(sym)
                }
            };
            let bf = get(i_b).parse::<f32>().unwrap_or(0.0);
            let seq = f[i_seq].parse::<i32>().unwrap_or(0);
            let ins = {
                let v = get(i_ins);
                if v == "?" || v == "." { "" } else { v }
            };

            b.push(
                f[i_asym], res_name, seq, ins, atom_name, element, [x, y, z], bf,
            );
        }
        break;
    }

    let s = b.finish(id, title, source);
    if s.atoms.is_empty() {
        anyhow::bail!("no atoms found in mmCIF file");
    }
    Ok(s)
}

/// Dispatch on content rather than filename: SimpleFold, RCSB and AlphaFold all
/// hand us slightly different things.
pub fn parse_any(text: &str, id: String, source: Source) -> anyhow::Result<Structure> {
    let looks_cif = text.contains("_atom_site.");
    if looks_cif {
        parse_mmcif(text, id, source)
    } else {
        parse_pdb(text, id, source)
    }
}
