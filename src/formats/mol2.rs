//! Tripos MOL2: atoms with SYBYL types, explicit bonds, and substructures
//! standing in for residues. A file of several molecules — docking poses,
//! usually — is read as frames when they all have the same atoms.

use std::collections::HashMap;

use anyhow::{bail, Result};

use super::{cell_from_parameters, element_or_number};
use crate::protein::{Builder, Record, Snapshot, Source, Structure};

#[derive(PartialEq)]
enum Section {
    Other,
    Molecule,
    Atom,
    Bond,
    Crysin,
}

pub fn parse(text: &str, id: String, source: Source) -> Result<Structure> {
    let mut b = Builder::new("MOL2");
    let mut title = String::new();
    let mut section = Section::Other;
    let mut molecule = 0usize;
    let mut line_in_section = 0usize;
    let mut ids: HashMap<String, u32> = HashMap::new();
    let mut frame: Vec<[f32; 3]> = Vec::new();

    for line in text.lines() {
        let t = line.trim();
        if let Some(name) = t.strip_prefix("@<TRIPOS>") {
            section = match name.to_ascii_uppercase().as_str() {
                "MOLECULE" => {
                    molecule += 1;
                    if molecule > 2 || (molecule == 2 && !frame.is_empty()) {
                        b.frames.push(Snapshot {
                            step: None,
                            pos: std::mem::take(&mut frame),
                            cell: None,
                        });
                    }
                    Section::Molecule
                }
                "ATOM" => Section::Atom,
                "BOND" => Section::Bond,
                "CRYSIN" => Section::Crysin,
                _ => Section::Other,
            };
            line_in_section = 0;
            continue;
        }
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        line_in_section += 1;
        let f: Vec<&str> = t.split_whitespace().collect();
        match section {
            Section::Molecule if molecule == 1 && line_in_section == 1 => {
                title = t.to_string();
            }
            Section::Atom => {
                if f.len() < 6 {
                    continue;
                }
                let (Ok(x), Ok(y), Ok(z)) = (
                    f[2].parse::<f32>(),
                    f[3].parse::<f32>(),
                    f[4].parse::<f32>(),
                ) else {
                    continue;
                };
                if molecule > 1 {
                    frame.push([x, y, z]);
                    continue;
                }
                // "C.ar", "N.am", "Fe": the element is before the dot.
                let sybyl = f[5];
                let symbol = sybyl.split('.').next().unwrap_or(sybyl);
                let mut element = element_or_number(symbol);
                if !element.is_known() {
                    element = element_or_number(f[1]);
                }
                // "ALA12" is residue ALA, number 12.
                let subst = f.get(7).copied().unwrap_or("");
                let res_name = subst.trim_end_matches(|c: char| c.is_ascii_digit());
                let res_seq = f.get(6).and_then(|v| v.parse().ok()).unwrap_or(0);
                let kind = b.kind(sybyl);
                ids.insert(f[0].to_string(), b.len() as u32);
                b.push(Record {
                    res_name,
                    res_seq,
                    name: f[1],
                    element,
                    pos: [x, y, z],
                    kind,
                    ..Record::default()
                });
            }
            Section::Bond if molecule == 1 && f.len() >= 3 => {
                if let (Some(&a), Some(&c)) = (ids.get(f[1]), ids.get(f[2])) {
                    b.bonds.push([a, c]);
                }
            }
            Section::Crysin if molecule == 1 && f.len() >= 6 => {
                let v: Vec<f32> = f[..6].iter().filter_map(|x| x.parse().ok()).collect();
                if v.len() == 6 {
                    b.cell = Some(cell_from_parameters(v[0], v[1], v[2], v[3], v[4], v[5]));
                }
            }
            _ => {}
        }
    }
    if !frame.is_empty() {
        b.frames.push(Snapshot {
            step: None,
            pos: frame,
            cell: None,
        });
    }

    if b.len() == 0 {
        bail!("no atoms found in MOL2 file");
    }
    Ok(b.finish(id, title, source))
}
