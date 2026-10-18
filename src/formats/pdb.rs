//! PDB, and PQR, its whitespace-separated cousin from electrostatics.
//!
//! Every MODEL after the first becomes a frame, so an NMR ensemble or a
//! trajectory written out by VMD can be stepped through. VMD separates its
//! frames with END rather than MODEL, which is honoured too.

use std::ops::Range;

use anyhow::{bail, Result};

use super::{cell_from_parameters, guess_element};
use crate::protein::{Builder, Element, Record, Snapshot, Source, Structure};

/// Water is dropped from crystal structures, where it is clutter. A
/// simulation's solvent is named differently and is kept.
fn is_crystal_water(res_name: &str) -> bool {
    res_name == "HOH" || res_name == "DOD"
}

pub fn parse(text: &str, id: String, source: Source) -> Result<Structure> {
    let mut b = Builder::new("PDB");
    b.assemblies = crate::assembly::from_pdb(text);
    let mut title = String::new();
    // Which model the lines belong to; past the first, only positions are kept.
    let mut model = 0usize;
    let mut frame: Vec<[f32; 3]> = Vec::new();
    let mut seen_atoms = false;

    for line in text.lines() {
        if line.starts_with("TITLE ") {
            let part = line.get(10..).unwrap_or("").trim();
            if !title.is_empty() {
                title.push(' ');
            }
            title.push_str(part);
            continue;
        }
        if line.starts_with("CRYST1") && model == 0 {
            let col = |r: Range<usize>| line.get(r).unwrap_or("").trim().parse::<f32>().ok();
            if let (Some(a), Some(bb), Some(c), Some(al), Some(be), Some(ga)) = (
                col(6..15),
                col(15..24),
                col(24..33),
                col(33..40),
                col(40..47),
                col(47..54),
            ) {
                // 1 x 1 x 1 is the placeholder NMR and EM entries carry.
                if a > 1.0 && bb > 1.0 && c > 1.0 {
                    b.cell = Some(cell_from_parameters(a, bb, c, al, be, ga));
                }
            }
            continue;
        }
        if line.starts_with("ENDMDL") || line.trim_end() == "END" {
            if seen_atoms {
                if model > 0 {
                    b.frames.push(Snapshot {
                        pos: std::mem::take(&mut frame),
                        cell: None,
                    });
                }
                model += 1;
                seen_atoms = false;
            }
            continue;
        }
        if !(line.starts_with("ATOM  ") || line.starts_with("HETATM")) {
            continue;
        }
        let col = |r: Range<usize>| line.get(r).unwrap_or("").trim();

        let res_name = col(17..20);
        if is_crystal_water(res_name) {
            continue;
        }
        // Keep only the first alternate location.
        let alt = line.as_bytes().get(16).copied().unwrap_or(b' ');
        if alt != b' ' && alt != b'A' {
            continue;
        }
        let (Ok(x), Ok(y), Ok(z)) = (
            col(30..38).parse::<f32>(),
            col(38..46).parse::<f32>(),
            col(46..54).parse::<f32>(),
        ) else {
            continue;
        };
        seen_atoms = true;
        if model > 0 {
            frame.push([x, y, z]);
            continue;
        }

        let atom_name = col(12..16);
        let chain = match col(21..22) {
            "" => "A",
            c => c,
        };
        let element = match col(76..78) {
            "" => guess_element(atom_name, res_name),
            sym => Element::from_symbol(sym),
        };
        b.push(Record {
            chain,
            res_name,
            res_seq: col(22..26).parse().unwrap_or(0),
            ins: col(26..27),
            name: atom_name,
            element,
            pos: [x, y, z],
            b: col(60..66).parse().unwrap_or(0.0),
            kind: 0,
        });
    }
    // A last model with no ENDMDL after it.
    if !frame.is_empty() {
        b.frames.push(Snapshot { pos: frame, cell: None });
    }

    if b.len() == 0 {
        bail!("no atoms found in PDB file");
    }
    Ok(b.finish(id, title, source))
}

/// PQR keeps PDB's record names but separates fields by whitespace, with a
/// charge and a radius where the occupancy and B-factor were. The chain is
/// optional, which is told from the number of fields.
pub fn parse_pqr(text: &str, id: String, source: Source) -> Result<Structure> {
    let mut b = Builder::new("PQR");
    for line in text.lines() {
        if line.starts_with("ENDMDL") {
            break;
        }
        if !(line.starts_with("ATOM") || line.starts_with("HETATM")) {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 10 {
            continue;
        }
        let n = f.len();
        let (Ok(x), Ok(y), Ok(z)) = (
            f[n - 5].parse::<f32>(),
            f[n - 4].parse::<f32>(),
            f[n - 3].parse::<f32>(),
        ) else {
            continue;
        };
        let has_chain = n >= 11;
        let (chain, seq) = if has_chain { (f[4], f[5]) } else { ("A", f[4]) };
        b.push(Record {
            chain,
            res_name: f[3],
            res_seq: seq.parse().unwrap_or(0),
            name: f[2],
            element: guess_element(f[2], f[3]),
            pos: [x, y, z],
            b: f[n - 2].parse().unwrap_or(0.0),
            ..Record::default()
        });
    }
    if b.len() == 0 {
        bail!("no atoms found in PQR file");
    }
    Ok(b.finish(id, String::new(), source))
}
