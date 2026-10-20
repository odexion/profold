//! Readers for the structure and trajectory formats that VMD and OVITO open.
//!
//! Each reader turns text into a `Structure` through `protein::Builder` and
//! knows nothing else. The format is chosen from the file name first and from
//! the content when the name says nothing, and gzip is undone on the way in.

mod cif;
mod cube;
mod gro;
mod lammps;
mod mol2;
mod pdb;
mod poscar;
mod sdf;
mod xsf;
mod xyz;
pub mod write;
pub mod traj;
mod bcif;

use std::io::Read;
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};

use crate::protein::{Cell, Element, Source, Structure};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Format {
    Pdb,
    Pqr,
    Cif,
    Xyz,
    Gro,
    Mol2,
    Sdf,
    LammpsData,
    LammpsDump,
    Poscar,
    Cube,
    Xsf,
}

impl Format {
    pub fn label(self) -> &'static str {
        match self {
            Format::Pdb => "PDB",
            Format::Pqr => "PQR",
            Format::Cif => "CIF",
            Format::Xyz => "XYZ",
            Format::Gro => "GROMACS",
            Format::Mol2 => "MOL2",
            Format::Sdf => "SDF",
            Format::LammpsData => "LAMMPS data",
            Format::LammpsDump => "LAMMPS dump",
            Format::Poscar => "VASP",
            Format::Cube => "Gaussian cube",
            Format::Xsf => "XSF",
        }
    }

    /// From the file name alone, ignoring a trailing `.gz`. LAMMPS and VASP
    /// files are as often told apart by their stem (`data.melt`,
    /// `dump.crack`, `CONTCAR`) as by an extension.
    pub fn from_name(name: &str) -> Option<Format> {
        let name = name.to_ascii_lowercase();
        let name = name.strip_suffix(".gz").unwrap_or(&name);
        let ext = name.rsplit_once('.').map(|(_, e)| e).unwrap_or("");
        let format = match ext {
            "pdb" | "ent" | "pdb1" => Format::Pdb,
            "pqr" => Format::Pqr,
            "cif" | "mmcif" | "mcif" => Format::Cif,
            "xyz" | "extxyz" => Format::Xyz,
            "gro" => Format::Gro,
            "mol2" | "ml2" => Format::Mol2,
            "sdf" | "sd" | "mol" => Format::Sdf,
            "data" | "lmp" | "lammps" => Format::LammpsData,
            "dump" | "lammpstrj" | "lammpsdump" => Format::LammpsDump,
            "vasp" | "poscar" | "contcar" => Format::Poscar,
            "cube" | "cub" => Format::Cube,
            "xsf" | "axsf" => Format::Xsf,
            _ => {
                let stem = name.rsplit('/').next().unwrap_or(name);
                if stem.starts_with("data.") {
                    Format::LammpsData
                } else if stem.starts_with("dump.") {
                    Format::LammpsDump
                } else if ["poscar", "contcar", "xdatcar"]
                    .iter()
                    .any(|p| stem.starts_with(p))
                {
                    Format::Poscar
                } else {
                    return None;
                }
            }
        };
        Some(format)
    }

    /// From the content, for files whose name says nothing.
    pub fn sniff(text: &str) -> Option<Format> {
        let head: String = text.chars().take(64 * 1024).collect();
        if head.contains("ITEM: TIMESTEP") {
            return Some(Format::LammpsDump);
        }
        if head.contains("@<TRIPOS>") {
            return Some(Format::Mol2);
        }
        if head.contains("_atom_site") {
            return Some(Format::Cif);
        }
        if head.contains("PRIMCOORD") || head.contains("PRIMVEC") {
            return Some(Format::Xsf);
        }
        if head.contains(" V2000") || head.contains(" V3000") || head.contains("$$$$") {
            return Some(Format::Sdf);
        }
        if head
            .lines()
            .any(|l| l.starts_with("ATOM  ") || l.starts_with("HETATM"))
        {
            return Some(Format::Pdb);
        }
        if head.lines().any(|l| {
            let l = l.trim_end();
            l.ends_with(" atoms") && l.split_whitespace().count() == 2
        }) && head.lines().any(|l| l.trim_start().starts_with("Atoms"))
        {
            return Some(Format::LammpsData);
        }
        let mut lines = head.lines();
        let first = lines.next().unwrap_or("").trim();
        if first.parse::<usize>().is_ok() {
            return Some(Format::Xyz);
        }
        // A GROMACS file has a title, then the atom count.
        if lines.next().is_some_and(|l| l.trim().parse::<usize>().is_ok()) {
            return Some(Format::Gro);
        }
        None
    }
}

/// Every extension the reader takes, for the open dialog's description.
pub const SUPPORTED: &str =
    "PDB, PQR, mmCIF, BinaryCIF, CIF, XYZ, GRO, MOL2, SDF, LAMMPS data and dump, POSCAR, cube, XSF, DCD, XTC";

/// Bytes, unzipped if they are gzipped.
pub(crate) fn gunzip(bytes: Vec<u8>) -> Result<Vec<u8>> {
    if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut out = Vec::new();
        flate2::read::MultiGzDecoder::new(&bytes[..])
            .read_to_end(&mut out)
            .context("could not decompress")?;
        Ok(out)
    } else {
        Ok(bytes)
    }
}

/// Whether a file is a binary trajectory, DCD or XTC, by its name; and if so
/// its atom count and frames.
pub fn trajectory(name: &str, bytes: &[u8]) -> Option<Result<(usize, Vec<crate::protein::Snapshot>)>> {
    let lower = name.to_ascii_lowercase();
    if lower.ends_with(".dcd") {
        Some(traj::dcd(bytes))
    } else if lower.ends_with(".xtc") {
        Some(traj::xtc(bytes))
    } else {
        None
    }
}

/// Read bytes as text, undoing gzip if that is what they are.
pub fn decode(bytes: Vec<u8>) -> Result<String> {
    let bytes = gunzip(bytes)?;
    Ok(match String::from_utf8(bytes) {
        Ok(s) => s,
        Err(e) => String::from_utf8_lossy(e.as_bytes()).into_owned(),
    })
}

/// Open a structure or trajectory from disk.
pub fn load_path(path: &Path) -> Result<Structure> {
    let bytes = std::fs::read(path).with_context(|| format!("could not read {}", path.display()))?;
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    load_bytes(&name, bytes)
}

/// Read a file's contents, given its name to tell the format by.
pub fn load_bytes(name: &str, bytes: Vec<u8>) -> Result<Structure> {
    // BinaryCIF is unpacked to mmCIF and read as that.
    let lower = name.to_ascii_lowercase();
    if lower.ends_with(".bcif") || lower.ends_with(".bcif.gz") {
        let bytes = gunzip(bytes)?;
        let text = bcif::to_mmcif(&bytes).with_context(|| format!("reading {name} as BinaryCIF"))?;
        let stem = lower.trim_end_matches(".gz").trim_end_matches(".bcif");
        let mut s = cif::parse(&text, stem.to_string(), Source::File)?;
        if s.atoms.is_empty() {
            bail!("no atoms in {name}");
        }
        s.format = Some("BinaryCIF");
        return Ok(s);
    }
    // A map has no atoms of its own: it is a box with values in it, to lay
    // over the structure shown or to look at alone.
    if crate::volume::is_map(name) {
        let volume = crate::volume::read_map(name, bytes)?;
        let stem = lower.trim_end_matches(".gz");
        let stem = stem.rsplit_once('.').map_or(stem, |(s, _)| s);
        let mut s = Structure::new(stem.to_string(), String::new(), Source::File);
        s.cell = Some(volume.cell());
        s.volume = Some(std::sync::Arc::new(volume));
        s.format = Some("map");
        return Ok(s);
    }
    if let Some(found) = trajectory(name, &bytes) {
        // On its own, a trajectory's atoms are anonymous beads.
        let (n, frames) = found?;
        return Ok(Structure::beads(name, n, frames));
    }
    let text = decode(bytes)?;
    let format = Format::from_name(name)
        .or_else(|| Format::sniff(&text))
        .ok_or_else(|| anyhow!("{name} is not a supported format ({SUPPORTED})"))?;
    let stem = name.strip_suffix(".gz").unwrap_or(name);
    let id = stem
        .rsplit_once('.')
        .map(|(s, _)| s)
        .filter(|s| !s.is_empty())
        .unwrap_or(stem)
        .to_string();
    let mut s = parse(&text, format, id, Source::File)
        .with_context(|| format!("reading {name} as {}", format.label()))?;
    if s.atoms.is_empty() {
        bail!("no atoms in {name}");
    }
    s.format = Some(format.label());
    Ok(s)
}

pub fn parse(text: &str, format: Format, id: String, source: Source) -> Result<Structure> {
    match format {
        Format::Pdb => pdb::parse(text, id, source),
        Format::Pqr => pdb::parse_pqr(text, id, source),
        Format::Cif => cif::parse(text, id, source),
        Format::Xyz => xyz::parse(text, id, source),
        Format::Gro => gro::parse(text, id, source),
        Format::Mol2 => mol2::parse(text, id, source),
        Format::Sdf => sdf::parse(text, id, source),
        Format::LammpsData => lammps::parse_data(text, id, source),
        Format::LammpsDump => lammps::parse_dump(text, id, source),
        Format::Poscar => poscar::parse(text, id, source),
        Format::Cube => cube::parse(text, id, source),
        Format::Xsf => xsf::parse(text, id, source),
    }
}

/// Dispatch on content rather than filename: SimpleFold, RCSB and AlphaFold all
/// hand us slightly different things.
pub fn parse_any(text: &str, id: String, source: Source) -> Result<Structure> {
    if text.contains("_atom_site.") {
        cif::parse(text, id, source)
    } else {
        pdb::parse(text, id, source)
    }
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

/// A cell from its lengths and angles in degrees, `a` along x and `b` in the
/// xy plane, which is the crystallographers' convention.
pub fn cell_from_parameters(a: f32, b: f32, c: f32, alpha: f32, beta: f32, gamma: f32) -> Cell {
    let (al, be, ga) = (alpha.to_radians(), beta.to_radians(), gamma.to_radians());
    let cx = c * be.cos();
    let cy = c * (al.cos() - be.cos() * ga.cos()) / ga.sin();
    let cz = (c * c - cx * cx - cy * cy).max(0.0).sqrt();
    Cell {
        origin: [0.0; 3],
        vectors: [
            [a, 0.0, 0.0],
            [b * ga.cos(), b * ga.sin(), 0.0],
            [cx, cy, cz],
        ],
    }
}

/// A number from a crystallographic file, which may carry its uncertainty in
/// brackets: "5.4310(2)".
pub fn number(s: &str) -> Option<f32> {
    let s = s.split('(').next().unwrap_or(s);
    s.parse().ok()
}

/// The element of an atom known only by its name and residue, as in GROMACS
/// and PQR files. Single-atom ions are named after themselves, in capitals,
/// which would otherwise read CL as carbon.
pub fn guess_element(atom: &str, residue: &str) -> Element {
    let atom = atom.trim();
    let residue = residue.trim().to_ascii_uppercase();
    let ion = match residue.trim_end_matches(['+', '-', '1', '2', '3']) {
        "NA" | "SOD" => Some("Na"),
        "CL" | "CLA" => Some("Cl"),
        "K" | "POT" => Some("K"),
        "MG" => Some("Mg"),
        "CA" | "CAL" => Some("Ca"),
        "ZN" => Some("Zn"),
        "FE" => Some("Fe"),
        "CU" => Some("Cu"),
        "MN" => Some("Mn"),
        "LI" => Some("Li"),
        "CS" => Some("Cs"),
        "BR" => Some("Br"),
        "RB" => Some("Rb"),
        "IOD" => Some("I"),
        _ => None,
    };
    if let Some(sym) = ion {
        if atom.eq_ignore_ascii_case(residue.trim_end_matches(['+', '-', '1', '2', '3']))
            || atom.to_ascii_uppercase().starts_with(&sym.to_ascii_uppercase())
        {
            return Element::from_symbol(sym);
        }
    }
    Element::from_label(atom)
}

/// The element named by a label, or by an atomic number written in its place.
pub fn element_or_number(label: &str) -> Element {
    match label.trim().parse::<u32>() {
        Ok(z) => Element::from_number(z),
        Err(_) => Element::from_label(label),
    }
}

#[cfg(test)]
mod tests;
