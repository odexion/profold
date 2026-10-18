//! Biological assemblies: the molecule as it is in the cell, built from the
//! chains a file deposits and the symmetry operations that copy them.
//!
//! A crystal's asymmetric unit is often a part of the real thing, one chain
//! of a tetramer, say, or two copies of a monomer; an assembly says which
//! chains to copy and how to move each copy.

use std::collections::HashMap;

use crate::protein::{Chain, Structure};

/// A rotation and translation, as three rows of a 3x4 matrix.
pub type Op = [[f32; 4]; 3];

pub const IDENTITY: Op = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]];

#[derive(Clone, Debug, PartialEq)]
pub struct Assembly {
    pub id: String,
    /// What the file calls it, such as "tetrameric", when it says.
    pub name: String,
    pub parts: Vec<Part>,
}

/// Some chains, each copied by each of some operations.
#[derive(Clone, Debug, PartialEq)]
pub struct Part {
    pub chains: Vec<String>,
    pub ops: Vec<Op>,
}

impl Assembly {
    /// How many chains it has, counting every copy.
    pub fn chain_count(&self, s: &Structure) -> usize {
        self.parts
            .iter()
            .map(|p| p.ops.len() * s.chains.iter().filter(|c| p.chains.contains(&c.id)).count())
            .sum()
    }

    /// Whether it is just the deposited chains where they are: then there is
    /// nothing to build.
    pub fn is_deposited(&self, s: &Structure) -> bool {
        let identity = |op: &Op| {
            (0..3).all(|i| (0..4).all(|j| (op[i][j] - IDENTITY[i][j]).abs() < 1e-4))
        };
        let all: Vec<&String> = s.chains.iter().map(|c| &c.id).collect();
        let covered: Vec<&String> = self.parts.iter().flat_map(|p| p.chains.iter()).collect();
        self.parts.iter().all(|p| p.ops.len() == 1 && identity(&p.ops[0]))
            && all.iter().all(|c| covered.contains(c))
            && covered.len() == all.len()
    }

    /// A word for it: the file's own, or how many chains it has.
    pub fn label(&self, s: &Structure) -> String {
        let name = self.name.trim().to_lowercase();
        if !name.is_empty() && name != "?" {
            // "tetrameric" reads better on a button as "tetramer".
            return capitalise(name.strip_suffix("ic").unwrap_or(&name));
        }
        match self.chain_count(s) {
            1 => "Monomer".into(),
            2 => "Dimer".into(),
            3 => "Trimer".into(),
            4 => "Tetramer".into(),
            n => format!("{n} chains"),
        }
    }
}

fn capitalise(w: &str) -> String {
    let mut c = w.chars();
    match c.next() {
        Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
        None => String::new(),
    }
}

/// Apply an operation to a point.
pub fn apply(op: &Op, p: [f32; 3]) -> [f32; 3] {
    std::array::from_fn(|i| op[i][0] * p[0] + op[i][1] * p[1] + op[i][2] * p[2] + op[i][3])
}

/// One operation after another: `a` applied to what `b` gives.
fn compose(a: &Op, b: &Op) -> Op {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let rot = (0..3).map(|k| a[i][k] * b[k][j]).sum::<f32>();
            if j == 3 {
                rot + a[i][3]
            } else {
                rot
            }
        })
    })
}

/// The operations an mmCIF `oper_expression` names: `1`, `1,2,5`, `1-60`,
/// or a product of parenthesised lists, `(1-60)(61)`, each of the first
/// applied after each of the second.
pub fn expand(expr: &str, ops: &HashMap<String, Op>) -> Vec<Op> {
    let list = |text: &str| -> Vec<Op> {
        let mut out = Vec::new();
        for item in text.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            match item.split_once('-') {
                Some((a, b)) => {
                    if let (Ok(a), Ok(b)) = (a.trim().parse::<i64>(), b.trim().parse::<i64>()) {
                        out.extend((a..=b).filter_map(|i| ops.get(&i.to_string()).copied()));
                    }
                }
                None => out.extend(ops.get(item).copied()),
            }
        }
        out
    };
    let groups: Vec<&str> = if expr.contains('(') {
        expr.split(['(', ')']).map(str::trim).filter(|g| !g.is_empty()).collect()
    } else {
        vec![expr]
    };
    let mut product = vec![IDENTITY];
    for g in groups.iter().rev() {
        let these = list(g);
        product = these.iter().flat_map(|a| product.iter().map(move |b| compose(a, b))).collect();
    }
    product
}

/// The assemblies of a PDB file, from its REMARK 350 records.
pub fn from_pdb(text: &str) -> Vec<Assembly> {
    let mut out: Vec<Assembly> = Vec::new();
    let mut rows: Vec<[f32; 4]> = Vec::new();
    for line in text.lines().filter(|l| l.starts_with("REMARK 350")) {
        let body = line.get(10..).unwrap_or("").trim();
        if let Some(id) = body.strip_prefix("BIOMOLECULE:") {
            out.push(Assembly { id: id.trim().to_string(), name: String::new(), parts: Vec::new() });
        } else if let Some(name) = body
            .strip_prefix("AUTHOR DETERMINED BIOLOGICAL UNIT:")
            .or_else(|| body.strip_prefix("SOFTWARE DETERMINED QUATERNARY STRUCTURE:"))
        {
            if let Some(a) = out.last_mut().filter(|a| a.name.is_empty()) {
                a.name = name.trim().to_string();
            }
        } else if let Some(chains) = body.strip_prefix("APPLY THE FOLLOWING TO CHAINS:") {
            if let Some(a) = out.last_mut() {
                a.parts.push(Part { chains: split_chains(chains), ops: Vec::new() });
            }
        } else if let Some(chains) = body.strip_prefix("AND CHAINS:") {
            if let Some(p) = out.last_mut().and_then(|a| a.parts.last_mut()) {
                p.chains.extend(split_chains(chains));
            }
        } else if body.starts_with("BIOMT") {
            let f: Vec<f32> = body.split_whitespace().skip(2).filter_map(|v| v.parse().ok()).collect();
            if f.len() >= 4 {
                rows.push([f[0], f[1], f[2], f[3]]);
            }
            if rows.len() == 3 {
                if let Some(p) = out.last_mut().and_then(|a| a.parts.last_mut()) {
                    p.ops.push([rows[0], rows[1], rows[2]]);
                }
                rows.clear();
            }
        }
    }
    out.retain(|a| a.parts.iter().any(|p| !p.ops.is_empty()));
    out
}

fn split_chains(text: &str) -> Vec<String> {
    text.split(',').map(str::trim).filter(|c| !c.is_empty()).map(str::to_string).collect()
}

/// The structure with assembly `a` built: every chain it names, copied by
/// each of its operations. The first copy keeps a chain's name; the rest
/// are numbered, A, A2, A3.
pub fn build(s: &Structure, a: &Assembly) -> Structure {
    let mut out = s.clone();
    out.atoms.clear();
    out.residues.clear();
    out.chains.clear();
    out.bonds.clear();
    out.frames = std::sync::Arc::new(Vec::new());
    out.frame = 0;
    out.cell = None;
    out.deviation.clear();
    // Built, it is its own whole: nothing further to build from it.
    out.assemblies.clear();
    let mut bonds_of: Vec<Vec<[u32; 2]>> = vec![Vec::new(); s.chains.len()];
    // Bonds are copied with their chain; those between chains are kept
    // only within one copy of a part.
    let chain_of_atom: Vec<usize> = s.atoms.iter().map(|a| s.residues[a.residue as usize].chain as usize).collect();
    let mut between: Vec<[u32; 2]> = Vec::new();
    for &[x, y] in &s.bonds {
        let (cx, cy) = (chain_of_atom[x as usize], chain_of_atom[y as usize]);
        if cx == cy {
            bonds_of[cx].push([x, y]);
        } else {
            between.push([x, y]);
        }
    }
    let mut seen: HashMap<String, usize> = HashMap::new();
    for part in &a.parts {
        let chains: Vec<usize> = (0..s.chains.len()).filter(|&c| part.chains.contains(&s.chains[c].id)).collect();
        for op in &part.ops {
            let mut index: HashMap<u32, u32> = HashMap::new();
            for &c in &chains {
                let chain = &s.chains[c];
                let n = seen.entry(chain.id.clone()).or_insert(0);
                *n += 1;
                let id = if *n == 1 { chain.id.clone() } else { format!("{}{}", chain.id, n) };
                let new_chain = out.chains.len() as u16;
                let mut residues = Vec::with_capacity(chain.residues.len());
                for &r in &chain.residues {
                    let old = &s.residues[r as usize];
                    let new_r = out.residues.len() as u32;
                    let mut res = old.clone();
                    res.chain = new_chain;
                    let atoms: Vec<usize> = s.residue_atoms(r as usize);
                    for &i in &atoms {
                        let mut atom = s.atoms[i];
                        atom.pos = apply(op, atom.pos);
                        atom.residue = new_r;
                        index.insert(i as u32, out.atoms.len() as u32);
                        out.atoms.push(atom);
                    }
                    let remap = |x: Option<u32>| x.and_then(|x| index.get(&x).copied());
                    res.ca = remap(old.ca);
                    res.n = remap(old.n);
                    res.c = remap(old.c);
                    res.o = remap(old.o);
                    res.nucleic = remap(old.nucleic);
                    out.residues.push(res);
                    residues.push(new_r);
                }
                out.chains.push(Chain { id, residues, color: chain.color });
                out.bonds.extend(bonds_of[c].iter().filter_map(|[x, y]| Some([*index.get(x)?, *index.get(y)?])));
            }
            out.bonds.extend(between.iter().filter_map(|[x, y]| Some([*index.get(x)?, *index.get(y)?])));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const DIMER: &str = "\
REMARK 350 BIOMOLECULE: 1
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000       10.00000
REMARK 350   BIOMT2   2  0.000000 -1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00 10.00           N
ATOM      2  CA  GLY A   1       2.400   2.000   3.000  1.00 10.00           C
ATOM      3  C   GLY A   1       3.000   3.300   3.000  1.00 10.00           C
ATOM      4  O   GLY A   1       2.400   4.300   3.000  1.00 10.00           O
";

    #[test]
    fn a_pdb_dimer_is_built_from_its_one_chain() {
        let s = crate::formats::parse(DIMER, crate::formats::Format::Pdb, "test".into(), crate::protein::Source::File).unwrap();
        assert_eq!(s.assemblies.len(), 1);
        let a = &s.assemblies[0];
        assert_eq!(a.label(&s), "Dimer");
        assert!(!a.is_deposited(&s));
        assert_eq!(a.chain_count(&s), 2);
        let built = build(&s, a);
        assert_eq!(built.atoms.len(), 8);
        assert_eq!(built.chains.iter().map(|c| c.id.as_str()).collect::<Vec<_>>(), ["A", "A2"]);
        // The copy is turned about z and moved along x.
        let copy = built.atoms[4].pos;
        assert!((copy[0] - 9.0).abs() < 1e-4 && (copy[1] + 2.0).abs() < 1e-4 && (copy[2] - 3.0).abs() < 1e-4);
        assert_eq!(built.residues[1].ca, Some(5));
        assert!(built.bonds.iter().all(|b| b[0] < 8 && b[1] < 8));
        assert_eq!(built.bonds.len(), s.bonds.len() * 2);
    }

    #[test]
    fn operation_expressions() {
        let mut ops = HashMap::new();
        for i in 1..=4 {
            let mut op = IDENTITY;
            op[0][3] = i as f32;
            ops.insert(i.to_string(), op);
        }
        assert_eq!(expand("1", &ops).len(), 1);
        assert_eq!(expand("1,2,4", &ops).len(), 3);
        assert_eq!(expand("(1-4)", &ops).len(), 4);
        let product = expand("(1,2)(3,4)", &ops);
        assert_eq!(product.len(), 4);
        // Shifts add: 1+3, 1+4, 2+3, 2+4.
        let shifts: Vec<f32> = product.iter().map(|o| o[0][3]).collect();
        assert_eq!(shifts, [4.0, 5.0, 5.0, 6.0]);
    }
}

#[cfg(test)]
mod real {
    use super::*;

    /// Entries downloaded from RCSB into `PROFOLD_DL`, when it is set.
    #[test]
    #[ignore]
    fn downloaded_entries() {
        let Some(dir) = std::env::var_os("PROFOLD_DL") else { return };
        for name in ["1stp.cif", "1stp.pdb", "1aon.cif", "1aon.pdb"] {
            let bytes = std::fs::read(std::path::Path::new(&dir).join(name)).unwrap();
            let s = crate::formats::load_bytes(name, bytes).unwrap();
            for a in &s.assemblies {
                let built = build(&s, a);
                eprintln!(
                    "{name}: {} chains, assembly {} '{}' -> {} chains, {} atoms (from {}), deposited: {}",
                    s.chains.len(), a.id, a.label(&s), built.chains.len(), built.atoms.len(), s.atoms.len(), a.is_deposited(&s)
                );
            }
        }
    }
}

#[cfg(test)]
mod nucleic_real {
    #[test]
    #[ignore]
    fn downloaded_strands() {
        let Some(dir) = std::env::var_os("PROFOLD_DL") else { return };
        for name in ["1bna.pdb", "1ehz.cif", "1tsr.pdb"] {
            let bytes = std::fs::read(std::path::Path::new(&dir).join(name)).unwrap();
            let s = crate::formats::load_bytes(name, bytes).unwrap();
            let traced = s.residues.iter().filter(|r| r.nucleic.is_some()).count();
            let ligands: Vec<String> = s.ligands().into_iter().map(|(n, k)| format!("{n}x{k}")).collect();
            eprintln!("{name}: {traced} nucleotides traced, polymer {}, ligands {ligands:?}", s.has_polymer());
        }
    }
}
