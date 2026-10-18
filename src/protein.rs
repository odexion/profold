//! The structure model.
//!
//! Everything downstream (renderer, sequence view, colouring) reads this and
//! nothing else, so the parsers in `formats` are the only places that know
//! about file formats. A protein and a box of argon are the same thing here:
//! atoms, grouped into residues and chains where the file says so, with bonds,
//! a periodic cell and further frames where it has them.

use std::collections::HashMap;
use std::sync::Arc;

pub use crate::elements::Element;

#[derive(Clone, Copy)]
pub struct Atom {
    pub pos: [f32; 3],
    pub element: Element,
    /// B-factor, or pLDDT for predicted models.
    pub b: f32,
    pub residue: u32,
    /// Index into `Structure::kinds`: the atom type in a simulation file.
    pub kind: u16,
    /// The atom's name in the file, such as `CA` or `OG1`, padded with
    /// zeros. Four bytes, as a PDB file allows; longer names are cut short.
    pub name: [u8; 4],
}

impl Atom {
    pub fn name(&self) -> &str {
        let n = self.name.iter().position(|&b| b == 0).unwrap_or(4);
        std::str::from_utf8(&self.name[..n]).unwrap_or("")
    }
}

/// An atom name packed for `Atom::name`, cut at a character boundary.
pub fn atom_name(name: &str) -> [u8; 4] {
    let mut out = [0; 4];
    let mut n = 0;
    for c in name.trim().chars() {
        let len = c.len_utf8();
        if n + len > 4 {
            break;
        }
        c.encode_utf8(&mut out[n..n + len]);
        n += len;
    }
    out
}

/// A periodic simulation cell: a corner and the three edge vectors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cell {
    pub origin: [f32; 3],
    pub vectors: [[f32; 3]; 3],
}

impl Cell {
    /// An orthorhombic box from its corners.
    pub fn from_bounds(lo: [f32; 3], hi: [f32; 3]) -> Self {
        Cell {
            origin: lo,
            vectors: [
                [hi[0] - lo[0], 0.0, 0.0],
                [0.0, hi[1] - lo[1], 0.0],
                [0.0, 0.0, hi[2] - lo[2]],
            ],
        }
    }

    pub fn lengths(&self) -> [f32; 3] {
        self.vectors
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
    }

    /// Fractional coordinates to Cartesian.
    pub fn to_cartesian(&self, f: [f32; 3]) -> [f32; 3] {
        let v = &self.vectors;
        std::array::from_fn(|k| {
            self.origin[k] + f[0] * v[0][k] + f[1] * v[1][k] + f[2] * v[2][k]
        })
    }

    /// The eight corners.
    pub fn corners(&self) -> [[f32; 3]; 8] {
        std::array::from_fn(|i| {
            self.to_cartesian([
                (i & 1) as f32,
                ((i >> 1) & 1) as f32,
                ((i >> 2) & 1) as f32,
            ])
        })
    }

    /// The twelve edges, as pairs of corners.
    pub fn edges(&self) -> [([f32; 3], [f32; 3]); 12] {
        let c = self.corners();
        const E: [(usize, usize); 12] = [
            (0, 1), (2, 3), (4, 5), (6, 7),
            (0, 2), (1, 3), (4, 6), (5, 7),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        E.map(|(a, b)| (c[a], c[b]))
    }
}

/// One frame of a trajectory: where every atom is, and the cell, which moves
/// too under constant pressure.
#[derive(Clone)]
pub struct Snapshot {
    pub pos: Vec<[f32; 3]>,
    pub cell: Option<Cell>,
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

    /// Neither an amino acid nor a nucleotide: a ligand, an ion or solvent.
    pub fn outside_polymer(&self) -> bool {
        self.ca.is_none()
            && !matches!(
                self.name.as_str(),
                "A" | "C" | "G" | "U" | "T" | "I" | "DA" | "DC" | "DG" | "DT" | "DU" | "DI"
            )
    }
}

#[derive(Clone)]
pub struct Chain {
    pub id: String,
    /// Indices into `Structure::residues`, in file order.
    pub residues: Vec<u32>,
    /// A colour of its own in place of the palette's, when colouring by
    /// chain. Kept with the structure, since chain A of one protein has
    /// nothing to do with chain A of the next.
    pub color: Option<[f32; 3]>,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Source {
    Pdb,
    AlphaFold,
    PubChem,
    #[cfg_attr(not(feature = "simplefold"), allow(dead_code))]
    SimpleFold,
    /// Opened from disk.
    File,
    /// Built into the app.
    Sample,
}

impl Source {
    pub fn label(self) -> &'static str {
        match self {
            Source::Pdb => "RCSB PDB",
            Source::AlphaFold => "AlphaFold DB",
            Source::PubChem => "PubChem",
            Source::SimpleFold => "SimpleFold",
            Source::File => "Local file",
            Source::Sample => "Sample",
        }
    }

    /// Whether the coordinates were measured or predicted.
    pub fn kind(self) -> &'static str {
        match self {
            Source::Pdb => "experimental",
            Source::PubChem => "compound",
            Source::AlphaFold | Source::SimpleFold => "predicted",
            Source::File => "file",
            Source::Sample => "sample",
        }
    }

    /// Read from a file rather than fetched or predicted.
    pub fn is_local(self) -> bool {
        matches!(self, Source::File | Source::Sample)
    }
}

/// How a kind of ligand is painted when given a colour of its own. Without
/// one, its atoms wear their elements' colours.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum LigandPaint {
    /// Carbons in the colour, other atoms in their elements': the usual
    /// way to set a ligand off in a figure while its chemistry still reads.
    Carbons([f32; 3]),
    /// Every atom in the colour.
    Solid([f32; 3]),
}

impl LigandPaint {
    pub fn color(self) -> [f32; 3] {
        match self {
            LigandPaint::Carbons(c) | LigandPaint::Solid(c) => c,
        }
    }

    pub fn with_color(self, c: [f32; 3]) -> Self {
        match self {
            LigandPaint::Carbons(_) => LigandPaint::Carbons(c),
            LigandPaint::Solid(_) => LigandPaint::Solid(c),
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
    /// Colours of their own for kinds of ligand, by residue name. Kept with
    /// the structure, like the chains' colours.
    pub ligand_paint: HashMap<String, LigandPaint>,
    /// Pairs of atom indices, read from the file or worked out from distance.
    pub bonds: Vec<[u32; 2]>,
    pub cell: Option<Cell>,
    /// Every frame of a trajectory, the current one included. Empty for a
    /// single structure. Shared, so changing frame does not copy the lot.
    pub frames: Arc<Vec<Snapshot>>,
    pub frame: usize,
    /// Atom type names, indexed by `Atom::kind`. Empty when the file has no
    /// types beyond elements.
    pub kinds: Vec<String>,
    /// The file format it was read from, for display.
    pub format: Option<&'static str>,
    /// Drawing radius for atoms of unknown element, set from the spacing of
    /// the atoms: a Lennard-Jones fluid in reduced units is one wide, not 1.6.
    pub unknown_radius: f32,
    /// The span of B-factors colouring runs over, worked out when first asked.
    b_span: std::sync::OnceLock<Option<(f32, f32)>>,
}

impl Structure {
    pub fn new(id: String, title: String, source: Source) -> Self {
        Structure {
            id,
            title,
            source,
            atoms: Vec::new(),
            residues: Vec::new(),
            chains: Vec::new(),
            has_confidence: false,
            deviation: Vec::new(),
            tint: None,
            ligand_paint: HashMap::new(),
            bonds: Vec::new(),
            cell: None,
            frames: Arc::new(Vec::new()),
            frame: 0,
            kinds: Vec::new(),
            format: None,
            unknown_radius: Element::UNKNOWN.vdw(),
            b_span: std::sync::OnceLock::new(),
        }
    }

    /// The B-factors of the polymer's residues, from the 5th to the 95th
    /// percentile, so that a floppy tail or two does not wash the rest out.
    /// `None` where there are none worth colouring by: a predicted model's
    /// column is pLDDT, and a flat column says nothing.
    pub fn b_span(&self) -> Option<(f32, f32)> {
        *self.b_span.get_or_init(|| {
            if self.has_confidence {
                return None;
            }
            let mut b: Vec<f32> = self
                .residues
                .iter()
                .filter(|r| !r.outside_polymer())
                .map(|r| r.b)
                .filter(|b| b.is_finite())
                .collect();
            if b.len() < 3 {
                return None;
            }
            b.sort_by(f32::total_cmp);
            let at = |q: f32| b[((b.len() - 1) as f32 * q).round() as usize];
            let (lo, hi) = (at(0.05), at(0.95));
            (hi - lo > 0.5).then_some((lo, hi))
        })
    }

    pub fn sequence(&self) -> String {
        self.residues
            .iter()
            .filter(|r| r.ca.is_some())
            .map(|r| r.code)
            .collect()
    }

    /// Whether there is a protein backbone to draw as a cartoon or a tube.
    /// The database's own page for what was fetched.
    pub fn page(&self) -> Option<String> {
        let id = self.id.trim();
        Some(match self.source {
            Source::Pdb => format!("https://www.rcsb.org/structure/{}", id.to_uppercase()),
            Source::AlphaFold => format!("https://alphafold.ebi.ac.uk/entry/{id}"),
            Source::PubChem => format!(
                "https://pubchem.ncbi.nlm.nih.gov/compound/{}",
                id.trim_start_matches("CID").trim()
            ),
            Source::SimpleFold | Source::File | Source::Sample => return None,
        })
    }

    /// The colour chain `i` is painted when colouring by chain.
    pub fn chain_color(&self, i: usize, palette: crate::render::Palette) -> [f32; 3] {
        self.chains.get(i).and_then(|c| c.color).unwrap_or(palette.chain(i))
    }

    pub fn has_polymer(&self) -> bool {
        self.residues.iter().filter(|r| r.ca.is_some()).take(3).count() >= 3
    }

    /// Each copy of a ligand, by residue name: its residues, in file order.
    pub fn copies(&self, name: &str) -> Vec<usize> {
        (0..self.residues.len())
            .filter(|&r| self.residues[r].name == name && self.residues[r].outside_polymer())
            .collect()
    }

    /// The atoms of a residue.
    pub fn residue_atoms(&self, residue: usize) -> Vec<usize> {
        (0..self.atoms.len()).filter(|&i| self.atoms[i].residue as usize == residue).collect()
    }

    /// The binding site of residue `of`: the polymer residues with an atom
    /// within `cutoff` ångström of one of its atoms, in file order.
    pub fn site(&self, of: usize, cutoff: f32) -> Vec<usize> {
        let own: Vec<[f32; 3]> = self.residue_atoms(of).iter().map(|&i| self.atoms[i].pos).collect();
        if own.is_empty() {
            return Vec::new();
        }
        let mut lo = [f32::MAX; 3];
        let mut hi = [f32::MIN; 3];
        for p in &own {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k] - cutoff);
                hi[k] = hi[k].max(p[k] + cutoff);
            }
        }
        let near = |q: [f32; 3]| {
            (0..3).all(|k| q[k] >= lo[k] && q[k] <= hi[k])
                && own.iter().any(|p| {
                    let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
                    d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= cutoff * cutoff
                })
        };
        let mut hit = vec![false; self.residues.len()];
        for a in &self.atoms {
            let r = a.residue as usize;
            if r != of && !hit[r] && !self.residues[r].outside_polymer() && near(a.pos) {
                hit[r] = true;
            }
        }
        (0..self.residues.len()).filter(|&r| hit[r]).collect()
    }

    /// Where a residue's side chain is, for writing its name beside it: the
    /// middle of the atoms past its backbone, or its alpha carbon for a
    /// glycine.
    pub fn side_chain_centre(&self, residue: usize) -> Option<[f32; 3]> {
        let r = self.residues.get(residue)?;
        let backbone = [r.n, r.c, r.o, r.ca];
        let side: Vec<[f32; 3]> = self
            .residue_atoms(residue)
            .into_iter()
            .filter(|&i| !backbone.contains(&Some(i as u32)))
            .map(|i| self.atoms[i].pos)
            .collect();
        if side.is_empty() {
            return r.ca.map(|ca| self.atoms[ca as usize].pos);
        }
        let n = side.len() as f32;
        let sum = side.iter().fold([0.0; 3], |s, p| [s[0] + p[0], s[1] + p[1], s[2] + p[2]]);
        Some([sum[0] / n, sum[1] / n, sum[2] / n])
    }

    /// A short name for an atom, as a measurement lists it: `SER42 OG`
    /// in a polymer or beside one, else its name or element and number.
    pub fn atom_tag(&self, atom: usize) -> String {
        let Some(a) = self.atoms.get(atom) else {
            return String::new();
        };
        let r = &self.residues[a.residue as usize];
        let name = match (a.name(), a.element.is_known()) {
            ("", true) => a.element.symbol(),
            ("", false) => "atom",
            (n, _) => n,
        };
        if !r.outside_polymer() || self.has_polymer() {
            format!("{}{} {name}", r.name, r.seq)
        } else {
            format!("{name} {}", atom + 1)
        }
    }

    /// What to say about an atom pointed at: a heading, and the lines under
    /// it. With `residue`, the residue it belongs to is what was pointed at,
    /// as on a ribbon, and its atom is not named.
    pub fn describe(&self, atom: usize, residue: bool) -> (String, Vec<String>) {
        let Some(a) = self.atoms.get(atom) else {
            return (String::new(), Vec::new());
        };
        let r = &self.residues[a.residue as usize];
        let chain = self.chains.get(r.chain as usize).map(|c| c.id.trim()).unwrap_or("");
        // A residue is worth naming in a polymer, or beside one; a lone
        // molecule's one residue, or a simulation's, is not.
        let named = !r.outside_polymer() || self.has_polymer();
        let symbol = a.element.is_known().then(|| a.element.symbol());
        let heading = if named {
            let mut h = format!("{} {}", r.name, r.seq);
            if !chain.is_empty() {
                h.push_str(&format!(" · chain {chain}"));
            }
            h
        } else {
            let label = match (a.name(), symbol) {
                ("", Some(e)) => e.to_string(),
                ("", None) => "Atom".to_string(),
                (n, _) => n.to_string(),
            };
            format!("{label} {}", atom + 1)
        };

        let mut lines = Vec::new();
        if residue {
            lines.push(
                match r.ss {
                    Ss::Helix => "helix",
                    Ss::Sheet => "strand",
                    Ss::Coil => "coil",
                }
                .to_string(),
            );
        } else {
            let mut what = Vec::new();
            if named && !a.name().is_empty() {
                what.push(a.name().to_string());
            }
            match symbol {
                Some(e) if what.first().is_none_or(|n| n != e) => what.push(e.to_string()),
                Some(_) => {}
                None => {
                    if let Some(kind) = self.kinds.get(a.kind as usize) {
                        what.push(format!("type {kind}"));
                    }
                }
            }
            if !what.is_empty() {
                lines.push(what.join(" · "));
            }
        }
        let b = if residue { r.b } else { a.b };
        if self.has_confidence {
            lines.push(format!("pLDDT {b:.1}"));
        } else if b != 0.0 && matches!(self.format, Some("PDB" | "mmCIF")) {
            lines.push(format!("B-factor {b:.1}"));
        }
        if !residue {
            let [x, y, z] = a.pos;
            lines.push(format!("{x:.2}, {y:.2}, {z:.2} Å"));
        }
        (heading, lines)
    }

    /// The kinds of ligand and ion beside a polymer, by residue name, with
    /// how many of each, in the order they first appear. Solvent is left out.
    pub fn ligands(&self) -> Vec<(String, usize)> {
        let mut out: Vec<(String, usize)> = Vec::new();
        if !self.has_polymer() {
            return out;
        }
        for r in self.residues.iter().filter(|r| r.outside_polymer() && !is_solvent(&r.name)) {
            match out.iter_mut().find(|(name, _)| *name == r.name) {
                Some((_, n)) => *n += 1,
                None => out.push((r.name.clone(), 1)),
            }
        }
        out
    }

    /// A polymer with molecules beside it, which are coloured by element
    /// whatever the polymer is coloured by.
    pub fn has_ligands(&self) -> bool {
        self.has_polymer() && self.residues.iter().any(|r| r.outside_polymer())
    }

    pub fn frame_count(&self) -> usize {
        self.frames.len().max(1)
    }

    /// Move every atom to where it is in frame `i`. Bonds and secondary
    /// structure are kept from the first frame's topology, but the helices are
    /// reassigned, since they do change over a simulation.
    pub fn set_frame(&mut self, i: usize) {
        let Some(snap) = self.frames.get(i).cloned() else {
            return;
        };
        if snap.pos.len() != self.atoms.len() {
            return;
        }
        for (a, p) in self.atoms.iter_mut().zip(&snap.pos) {
            a.pos = *p;
        }
        if snap.cell.is_some() {
            self.cell = snap.cell;
        }
        self.frame = i;
        if self.has_polymer() {
            crate::dssp::assign(self);
        }
    }

    /// Element counts, most common first, for the info panel.
    pub fn composition(&self) -> Vec<(Element, usize)> {
        let mut counts = [0usize; 256];
        for a in &self.atoms {
            counts[a.element.0 as usize] += 1;
        }
        let mut out: Vec<(Element, usize)> = counts
            .iter()
            .enumerate()
            .filter(|(_, n)| **n > 0)
            .map(|(z, n)| (Element(z as u8), *n))
            .collect();
        out.sort_by(|a, b| b.1.cmp(&a.1).then(a.0 .0.cmp(&b.0 .0)));
        out
    }

    /// Molecular weight in dalton from the atoms themselves.
    pub fn atomic_mass(&self) -> f32 {
        self.atoms.iter().map(|a| a.element.mass()).sum()
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

    /// Molecular weight in kilodalton, from the sequence rather than the
    /// atoms: most files carry no hydrogens, so adding up what is in the
    /// coordinates would come out roughly seven per cent light.
    pub fn mass_kda(&self) -> f32 {
        let residues: f32 = self.residues.iter().map(|r| residue_mass(r.code)).sum();
        // One water is released per peptide bond, so each chain keeps one.
        let water = 18.0153 * self.chains.len().max(1) as f32;
        (residues + water) / 1000.0
    }

    /// The chains worth folding, one entry per distinct sequence: the ids of
    /// the chains that share it, and the sequence itself. A tetramer of two
    /// pairs comes back as two entries, since folding the same sequence twice
    /// would return the same answer twice.
    #[cfg_attr(not(feature = "simplefold"), allow(dead_code))]
    pub fn distinct_chains(&self) -> Vec<(Vec<String>, String)> {
        let mut out: Vec<(Vec<String>, String)> = Vec::new();
        for (i, chain) in self.chains.iter().enumerate() {
            let sequence: String = chain
                .residues
                .iter()
                .map(|&ri| &self.residues[ri as usize])
                .filter(|r| r.ca.is_some() && r.code != 'X')
                .map(|r| r.code)
                .collect();
            if sequence.is_empty() {
                continue;
            }
            let id = if chain.id.trim().is_empty() {
                (i + 1).to_string()
            } else {
                chain.id.clone()
            };
            match out.iter_mut().find(|(_, seq)| *seq == sequence) {
                Some((ids, _)) => ids.push(id),
                None => out.push((vec![id], sequence)),
            }
        }
        out
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

/// Average mass of a residue in the chain, in dalton. Anything unrecognised
/// counts as the mean over the twenty, which is close enough for a total.
fn residue_mass(code: char) -> f32 {
    match code {
        'A' => 71.0788,
        'R' => 156.1875,
        'N' => 114.1038,
        'D' => 115.0886,
        'C' => 103.1388,
        'E' => 129.1155,
        'Q' => 128.1307,
        'G' => 57.0519,
        'H' => 137.1411,
        'I' => 113.1594,
        'L' => 113.1594,
        'K' => 128.1741,
        'M' => 131.1926,
        'F' => 147.1766,
        'P' => 97.1167,
        'S' => 87.0782,
        'T' => 101.1051,
        'W' => 186.2132,
        'Y' => 163.1760,
        'V' => 99.1326,
        _ => 110.0,
    }
}

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

/// One atom as a parser hands it over. Formats without residues or chains
/// leave those fields at their defaults, and everything lands in one residue.
pub struct Record<'a> {
    pub chain: &'a str,
    pub res_name: &'a str,
    pub res_seq: i32,
    pub ins: &'a str,
    pub name: &'a str,
    pub element: Element,
    pub pos: [f32; 3],
    pub b: f32,
    pub kind: u16,
}

impl Default for Record<'_> {
    fn default() -> Self {
        Record {
            chain: "",
            res_name: "",
            res_seq: 0,
            ins: "",
            name: "",
            element: Element::UNKNOWN,
            pos: [0.0; 3],
            b: 0.0,
            kind: 0,
        }
    }
}

/// Accumulates atoms into residues and chains regardless of source format.
#[derive(Default)]
pub struct Builder {
    atoms: Vec<Atom>,
    residues: Vec<Residue>,
    chains: Vec<Chain>,
    cur_chain: usize,
    cur_res_key: Option<(usize, String, i32, String)>,
    b_sum: f32,
    b_n: f32,
    /// Bonds the file lists. When it lists none they are worked out.
    pub bonds: Vec<[u32; 2]>,
    pub cell: Option<Cell>,
    /// Frames after the first; the first is the atoms themselves.
    pub frames: Vec<Snapshot>,
    pub kinds: Vec<String>,
    pub format: Option<&'static str>,
}

impl Builder {
    pub fn new(format: &'static str) -> Self {
        Self {
            format: Some(format),
            ..Default::default()
        }
    }

    pub fn len(&self) -> usize {
        self.atoms.len()
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

    /// The index of a named atom type, adding it if new.
    pub fn kind(&mut self, name: &str) -> u16 {
        match self.kinds.iter().position(|k| k == name) {
            Some(i) => i as u16,
            None => {
                self.kinds.push(name.to_string());
                (self.kinds.len() - 1) as u16
            }
        }
    }

    pub fn push(&mut self, r: Record) {
        // Ligand records revisit chains that already appeared, so look the
        // chain up by name rather than assuming one contiguous block per chain.
        let chain = if self.chains.get(self.cur_chain).is_some_and(|c| c.id == r.chain) {
            self.cur_chain
        } else {
            match self.chains.iter().position(|c| c.id == r.chain) {
                Some(i) => i,
                None => {
                    self.chains.push(Chain {
                        id: r.chain.to_string(),
                        residues: Vec::new(),
                        color: None,
                    });
                    self.chains.len() - 1
                }
            }
        };
        self.cur_chain = chain;

        let same = self.cur_res_key.as_ref().is_some_and(|(c, name, seq, ins)| {
            *c == chain && name == r.res_name && *seq == r.res_seq && ins == r.ins
        });
        if !same {
            self.flush_residue();
            self.cur_res_key = Some((chain, r.res_name.to_string(), r.res_seq, r.ins.to_string()));
            self.chains[chain].residues.push(self.residues.len() as u32);
            self.residues.push(Residue {
                name: r.res_name.to_string(),
                seq: r.res_seq,
                chain: chain as u16,
                ca: None,
                n: None,
                c: None,
                o: None,
                b: r.b,
                code: three_to_one(r.res_name),
                ss: Ss::Coil,
            });
        }

        let residue = (self.residues.len() - 1) as u32;
        let idx = self.atoms.len() as u32;
        {
            let res = &mut self.residues[residue as usize];
            // Backbone atoms only count in an amino acid: a calcium named CA
            // in a ligand is not an alpha carbon.
            if res.code != 'X' {
                match (r.name.trim(), r.element) {
                    ("CA", Element::C) => res.ca = Some(idx),
                    ("N", Element::N) => res.n = Some(idx),
                    ("C", Element::C) => res.c = Some(idx),
                    ("O", Element::O) => res.o = Some(idx),
                    _ => {}
                }
            }
        }
        self.b_sum += r.b;
        self.b_n += 1.0;
        self.atoms.push(Atom {
            pos: r.pos,
            element: r.element,
            b: r.b,
            residue,
            kind: r.kind,
            name: atom_name(r.name),
        });
    }

    pub fn finish(mut self, id: String, title: String, source: Source) -> Structure {
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

        let mut structure = Structure::new(id, title, source);
        structure.has_confidence = has_confidence;
        structure.format = self.format;
        structure.cell = self.cell;
        structure.kinds = self.kinds;
        structure.bonds = if self.bonds.is_empty() {
            crate::bonds::infer(&self.atoms)
        } else {
            self.bonds
        };
        structure.unknown_radius = crate::bonds::spacing_radius(&self.atoms);
        if !self.frames.is_empty() {
            let first = Snapshot {
                pos: self.atoms.iter().map(|a| a.pos).collect(),
                cell: self.cell,
            };
            let n = self.atoms.len();
            let mut frames = vec![first];
            // A trajectory whose atom count changes partway is cut there.
            frames.extend(self.frames.into_iter().take_while(|f| f.pos.len() == n));
            if frames.len() > 1 {
                structure.frames = Arc::new(frames);
            }
        }
        structure.atoms = self.atoms;
        structure.residues = self.residues;
        structure.chains = self.chains;
        if structure.has_polymer() {
            crate::dssp::assign(&mut structure);
        }
        structure
    }
}

/// Residue names that are solvent. Hidden under a cartoon, where a box of
/// water would bury the protein; shown in the atomistic styles.
pub fn is_solvent(name: &str) -> bool {
    matches!(
        name,
        "HOH" | "WAT" | "SOL" | "TIP3" | "TIP3P" | "TIP4" | "TIP4P" | "TIP5P" | "SPC" | "SPCE"
            | "T3P" | "T4P" | "H2O" | "DOD"
    )
}
