//! LAMMPS data files and text dumps.
//!
//! Neither says what element an atom is; both give a numeric type. A data
//! file's Masses section usually settles it, and a dump with an `element` or
//! `mass` column does too. Otherwise atoms stay typed but unknown, coloured by
//! type and sized from their spacing, as OVITO does.

use std::collections::HashMap;

use anyhow::{bail, Result};

use crate::protein::{Builder, Cell, Element, Record, Snapshot, Source, Structure};

fn strip_comment(line: &str) -> (&str, &str) {
    match line.split_once('#') {
        Some((data, comment)) => (data.trim(), comment.trim()),
        None => (line.trim(), ""),
    }
}

/// The box of a data file or a dump, with its tilt factors.
fn triclinic(lo: [f32; 3], hi: [f32; 3], tilt: [f32; 3]) -> Cell {
    let [xy, xz, yz] = tilt;
    Cell {
        origin: lo,
        vectors: [
            [hi[0] - lo[0], 0.0, 0.0],
            [xy, hi[1] - lo[1], 0.0],
            [xz, yz, hi[2] - lo[2]],
        ],
    }
}

// ---------------------------------------------------------------------------
// Data files
// ---------------------------------------------------------------------------

/// Which columns of an Atoms line hold the type, molecule and position.
struct Style {
    kind: usize,
    mol: Option<usize>,
    x: usize,
}

fn style(name: &str, sample: &[&str]) -> Style {
    let s = |kind, mol, x| Style { kind, mol, x };
    match name {
        "atomic" => s(1, None, 2),
        "charge" => s(1, None, 3),
        "molecular" | "bond" | "angle" => s(2, Some(1), 3),
        "full" => s(2, Some(1), 4),
        "sphere" => s(1, None, 4),
        _ => {
            // No style named: go by the number of columns, with three image
            // flags possibly on the end.
            let n = sample.len();
            match n {
                5 | 8 => s(1, None, 2),
                6 | 9 => {
                    // charge has a float third; molecular an integer type.
                    if sample[2].parse::<i64>().is_ok() {
                        s(2, Some(1), 3)
                    } else {
                        s(1, None, 3)
                    }
                }
                _ => s(2, Some(1), 4),
            }
        }
    }
}

pub fn parse_data(text: &str, id: String, source: Source) -> Result<Structure> {
    let mut lines = text.lines();
    let title = lines.next().unwrap_or("").trim().to_string();

    let mut lo = [0.0f32; 3];
    let mut hi = [0.0f32; 3];
    let mut tilt = [0.0f32; 3];
    let mut has_box = false;
    let mut masses: Vec<(String, f32, String)> = Vec::new();
    let mut atom_lines: Vec<(Vec<String>, String)> = Vec::new();
    let mut bond_lines: Vec<[String; 2]> = Vec::new();
    let mut section = String::new();
    let mut atom_style = String::new();

    for raw in lines {
        let (line, comment) = strip_comment(raw);
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        // A section header is a word or two with no numbers in front.
        if f[0].parse::<f64>().is_err() {
            section = f[0].to_string();
            if section == "Atoms" {
                atom_style = comment.to_string();
            }
            continue;
        }
        match section.as_str() {
            "" => match f.as_slice() {
                [a, b, "xlo", "xhi"] => {
                    lo[0] = a.parse()?;
                    hi[0] = b.parse()?;
                    has_box = true;
                }
                [a, b, "ylo", "yhi"] => {
                    lo[1] = a.parse()?;
                    hi[1] = b.parse()?;
                }
                [a, b, "zlo", "zhi"] => {
                    lo[2] = a.parse()?;
                    hi[2] = b.parse()?;
                }
                [a, b, c, "xy", "xz", "yz"] => {
                    tilt = [a.parse()?, b.parse()?, c.parse()?];
                }
                _ => {}
            },
            "Masses" if f.len() >= 2 => {
                masses.push((f[0].to_string(), f[1].parse().unwrap_or(0.0), comment.to_string()));
            }
            "Atoms" => {
                atom_lines.push((f.iter().map(|s| s.to_string()).collect(), comment.to_string()));
            }
            "Bonds" if f.len() >= 4 => {
                bond_lines.push([f[2].to_string(), f[3].to_string()]);
            }
            _ => {}
        }
    }
    if atom_lines.is_empty() {
        bail!("no Atoms section in LAMMPS data file");
    }

    let sample: Vec<&str> = atom_lines[0].0.iter().map(String::as_str).collect();
    let st = style(atom_style.split_whitespace().next().unwrap_or(""), &sample);

    // A type is an element when the Masses comment names one ("# Cu") or
    // when its mass is an element's.
    let element_of: HashMap<String, Element> = masses
        .iter()
        .map(|(t, m, comment)| {
            let named = Element::from_symbol(comment.split_whitespace().next().unwrap_or(""));
            let e = if named.is_known() { named } else { Element::from_mass(*m) };
            (t.clone(), e)
        })
        .collect();

    // Atom ids are not in order, and bonds refer to them by id. Sorting
    // keeps molecules together and the frame order stable.
    let mut rows: Vec<(i64, Vec<String>)> = atom_lines
        .into_iter()
        .filter_map(|(f, _)| Some((f[0].parse().ok()?, f)))
        .collect();
    rows.sort_by_key(|(id, _)| *id);

    let mut b = Builder::new("LAMMPS data");
    if has_box {
        b.cell = Some(triclinic(lo, hi, tilt));
    }
    let mut index: HashMap<i64, u32> = HashMap::with_capacity(rows.len());
    for (atom_id, f) in &rows {
        let (Some(x), Some(y), Some(z)) = (
            f.get(st.x).and_then(|v| v.parse::<f32>().ok()),
            f.get(st.x + 1).and_then(|v| v.parse::<f32>().ok()),
            f.get(st.x + 2).and_then(|v| v.parse::<f32>().ok()),
        ) else {
            continue;
        };
        let type_name = f.get(st.kind).map(String::as_str).unwrap_or("1");
        let kind = b.kind(type_name);
        let mol = st
            .mol
            .and_then(|c| f.get(c))
            .and_then(|v| v.parse().ok())
            .unwrap_or(0);
        index.insert(*atom_id, b.len() as u32);
        b.push(Record {
            res_seq: mol,
            name: type_name,
            element: element_of.get(type_name).copied().unwrap_or_default(),
            pos: [x, y, z],
            kind,
            ..Record::default()
        });
    }
    b.bonds = bond_lines
        .iter()
        .filter_map(|[a, c]| {
            Some([
                *index.get(&a.parse().ok()?)?,
                *index.get(&c.parse().ok()?)?,
            ])
        })
        .collect();
    Ok(b.finish(id, title, source))
}

// ---------------------------------------------------------------------------
// Dumps
// ---------------------------------------------------------------------------

struct DumpFrame {
    step: Option<i64>,
    cell: Cell,
    /// Sorted by atom id.
    atoms: Vec<DumpAtom>,
}

struct DumpAtom {
    id: i64,
    kind: String,
    element: Element,
    mol: i32,
    pos: [f32; 3],
}

/// Columns of interest in an `ITEM: ATOMS` header.
struct Columns {
    id: Option<usize>,
    kind: Option<usize>,
    element: Option<usize>,
    mass: Option<usize>,
    mol: Option<usize>,
    pos: [usize; 3],
    scaled: bool,
}

fn columns(names: &[&str]) -> Result<Columns> {
    let at = |n: &str| names.iter().position(|c| *c == n);
    let pick = |sets: &[[&str; 3]]| {
        sets.iter().find_map(|s| Some([at(s[0])?, at(s[1])?, at(s[2])?]))
    };
    let (pos, scaled) = if let Some(p) = pick(&[["x", "y", "z"], ["xu", "yu", "zu"]]) {
        (p, false)
    } else if let Some(p) = pick(&[["xs", "ys", "zs"], ["xsu", "ysu", "zsu"]]) {
        (p, true)
    } else {
        bail!("the dump has no position columns");
    };
    Ok(Columns {
        id: at("id"),
        kind: at("type"),
        element: at("element"),
        mass: at("mass"),
        mol: at("mol"),
        pos,
        scaled,
    })
}

fn read_frames(text: &str) -> Result<Vec<DumpFrame>> {
    let mut frames = Vec::new();
    let mut lines = text.lines();
    let mut natoms = 0usize;
    let mut step = None;
    let mut cell = Cell::from_bounds([0.0; 3], [1.0; 3]);

    while let Some(line) = lines.next() {
        let Some(item) = line.trim().strip_prefix("ITEM:") else {
            continue;
        };
        let item = item.trim();
        if item.starts_with("TIMESTEP") {
            step = lines.next().and_then(|l| l.trim().parse().ok());
        } else if item.starts_with("NUMBER OF ATOMS") {
            natoms = lines.next().unwrap_or("").trim().parse().unwrap_or(0);
        } else if item.starts_with("BOX BOUNDS") {
            let tilted = item.contains("xy");
            let mut b = [[0.0f32; 3]; 3];
            for row in &mut b {
                let v: Vec<f32> = lines
                    .next()
                    .unwrap_or("")
                    .split_whitespace()
                    .filter_map(|x| x.parse().ok())
                    .collect();
                for (k, x) in v.iter().take(3).enumerate() {
                    row[k] = *x;
                }
            }
            cell = if tilted {
                // The bounds given are those of the tilted box's bounding box.
                let (xy, xz, yz) = (b[0][2], b[1][2], b[2][2]);
                let xlo = b[0][0] - 0.0f32.min(xy).min(xz).min(xy + xz);
                let xhi = b[0][1] - 0.0f32.max(xy).max(xz).max(xy + xz);
                let ylo = b[1][0] - 0.0f32.min(yz);
                let yhi = b[1][1] - 0.0f32.max(yz);
                triclinic([xlo, ylo, b[2][0]], [xhi, yhi, b[2][1]], [xy, xz, yz])
            } else {
                Cell::from_bounds([b[0][0], b[1][0], b[2][0]], [b[0][1], b[1][1], b[2][1]])
            };
        } else if let Some(header) = item.strip_prefix("ATOMS") {
            let names: Vec<&str> = header.split_whitespace().collect();
            let cols = columns(&names)?;
            // Only the first frame's types matter; later frames are positions.
            let detail = frames.is_empty();
            let mut atoms = Vec::with_capacity(natoms);
            for (i, line) in lines.by_ref().take(natoms).enumerate() {
                let f: Vec<&str> = line.split_whitespace().collect();
                let num = |c: usize| f.get(c).and_then(|v| v.parse::<f32>().ok());
                let (Some(x), Some(y), Some(z)) =
                    (num(cols.pos[0]), num(cols.pos[1]), num(cols.pos[2]))
                else {
                    continue;
                };
                let pos = if cols.scaled {
                    cell.to_cartesian([x, y, z])
                } else {
                    [x, y, z]
                };
                if !detail {
                    atoms.push(DumpAtom {
                        id: cols
                            .id
                            .and_then(|c| f.get(c))
                            .and_then(|v| v.parse().ok())
                            .unwrap_or(i as i64),
                        kind: String::new(),
                        element: Element::UNKNOWN,
                        mol: 0,
                        pos,
                    });
                    continue;
                }
                let element = cols
                    .element
                    .and_then(|c| f.get(c))
                    .map(|e| Element::from_symbol(e))
                    .or_else(|| cols.mass.and_then(num).map(Element::from_mass))
                    .unwrap_or_default();
                atoms.push(DumpAtom {
                    id: cols
                        .id
                        .and_then(|c| f.get(c))
                        .and_then(|v| v.parse().ok())
                        .unwrap_or(i as i64),
                    kind: cols
                        .kind
                        .and_then(|c| f.get(c))
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| "1".into()),
                    element,
                    mol: cols
                        .mol
                        .and_then(|c| f.get(c))
                        .and_then(|v| v.parse().ok())
                        .unwrap_or(0),
                    pos,
                });
            }
            atoms.sort_by_key(|a| a.id);
            frames.push(DumpFrame { step: step.take(), cell, atoms });
        }
    }
    Ok(frames)
}

pub fn parse_dump(text: &str, id: String, source: Source) -> Result<Structure> {
    let mut frames = read_frames(text)?.into_iter();
    let Some(first) = frames.next() else {
        bail!("no frames in LAMMPS dump");
    };
    let mut b = Builder::new("LAMMPS dump");
    b.cell = Some(first.cell);
    b.first_step = first.step;
    for a in &first.atoms {
        let kind = b.kind(&a.kind);
        b.push(Record {
            res_seq: a.mol,
            name: &a.kind,
            element: a.element,
            pos: a.pos,
            kind,
            ..Record::default()
        });
    }
    let ids: Vec<i64> = first.atoms.iter().map(|a| a.id).collect();
    for f in frames {
        // Atoms lost or gained partway end the trajectory there.
        if f.atoms.len() != ids.len() || f.atoms.iter().zip(&ids).any(|(a, id)| a.id != *id) {
            break;
        }
        b.frames.push(Snapshot {
            step: f.step,
            pos: f.atoms.iter().map(|a| a.pos).collect(),
            cell: Some(f.cell),
        });
    }
    if b.len() == 0 {
        bail!("no atoms in LAMMPS dump");
    }
    Ok(b.finish(id, String::new(), source))
}
