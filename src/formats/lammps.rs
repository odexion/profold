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
    /// In file order, which LAMMPS does not keep from frame to frame.
    ids: Vec<i64>,
    pos: Vec<[f32; 3]>,
    /// Types, elements and molecules, read in the first frame only.
    detail: Vec<Detail>,
}

struct Detail {
    kind: String,
    element: Element,
    mol: i32,
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

impl Columns {
    /// The last column read, so the rest of a line can be skipped.
    fn last(&self) -> usize {
        [self.id, self.kind, self.element, self.mass, self.mol]
            .into_iter()
            .flatten()
            .chain(self.pos)
            .max()
            .unwrap_or(0)
    }
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

/// The dump cut where each frame begins, so frames can be read in parallel.
/// Text before the first `ITEM: TIMESTEP` stays with the first piece.
fn split_frames(text: &str) -> Vec<&str> {
    const MARK: &str = "ITEM: TIMESTEP";
    let mut starts: Vec<usize> = text
        .match_indices(MARK)
        .map(|(i, _)| i)
        .filter(|&i| i == 0 || text.as_bytes()[i - 1] == b'\n')
        .collect();
    if starts.first() != Some(&0) {
        starts.insert(0, 0);
    }
    starts.push(text.len());
    starts.windows(2).map(|w| &text[w[0]..w[1]]).collect()
}

/// The frames in a piece of a dump; `detail` reads types for the first.
fn read_frames(text: &str, mut detail: bool) -> Result<Vec<DumpFrame>> {
    let mut frames = Vec::new();
    let mut lines = text.lines();
    let mut natoms = 0usize;
    let mut step = None;
    let mut cell = Cell::from_bounds([0.0; 3], [1.0; 3]);
    let mut f: Vec<&str> = Vec::new();

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
            let last = cols.last();
            let mut frame = DumpFrame {
                step: step.take(),
                cell,
                ids: Vec::with_capacity(natoms),
                pos: Vec::with_capacity(natoms),
                detail: Vec::new(),
            };
            for (i, line) in lines.by_ref().take(natoms).enumerate() {
                f.clear();
                f.extend(line.split_whitespace().take(last + 1));
                let num = |c: usize| f.get(c).and_then(|v| v.parse::<f32>().ok());
                let (Some(x), Some(y), Some(z)) =
                    (num(cols.pos[0]), num(cols.pos[1]), num(cols.pos[2]))
                else {
                    continue;
                };
                frame.pos.push(if cols.scaled {
                    cell.to_cartesian([x, y, z])
                } else {
                    [x, y, z]
                });
                frame.ids.push(
                    cols.id
                        .and_then(|c| f.get(c))
                        .and_then(|v| v.parse().ok())
                        .unwrap_or(i as i64),
                );
                if !detail {
                    continue;
                }
                let element = cols
                    .element
                    .and_then(|c| f.get(c))
                    .map(|e| Element::from_symbol(e))
                    .or_else(|| cols.mass.and_then(num).map(Element::from_mass))
                    .unwrap_or_default();
                frame.detail.push(Detail {
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
                });
            }
            // Only the first frame's types matter; later frames are positions.
            detail = false;
            frames.push(frame);
        }
    }
    Ok(frames)
}

/// The order that sorts atoms by id.
fn by_id(ids: &[i64]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..ids.len()).collect();
    order.sort_by_key(|&i| ids[i]);
    order
}

pub fn parse_dump(text: &str, id: String, source: Source) -> Result<Structure> {
    use rayon::prelude::*;

    let pieces = split_frames(text);
    let mut head = read_frames(pieces[0], true)?.into_iter();
    let Some(first) = head.next() else {
        bail!("no frames in LAMMPS dump");
    };
    // Atoms go in by id: that keeps molecules together, and is the one
    // order every frame can be put back into.
    let order = by_id(&first.ids);
    let sorted: Vec<i64> = order.iter().map(|&i| first.ids[i]).collect();
    // A frame in the first frame's file order, as LAMMPS usually writes
    // them, needs no sorting of its own. Atoms lost or gained partway end
    // the trajectory there.
    let snapshot = |f: DumpFrame| -> Option<Snapshot> {
        let own;
        let order = if f.ids == first.ids {
            &order
        } else {
            own = by_id(&f.ids);
            if own.len() != sorted.len() || own.iter().zip(&sorted).any(|(&i, id)| f.ids[i] != *id) {
                return None;
            }
            &own
        };
        Some(Snapshot {
            step: f.step,
            pos: order.iter().map(|&i| f.pos[i]).collect(),
            cell: Some(f.cell),
        })
    };
    let rest: Vec<Vec<Option<Snapshot>>> = pieces[1..]
        .par_iter()
        .map(|piece| Ok(read_frames(piece, false)?.into_iter().map(snapshot).collect()))
        .collect::<Result<_>>()?;
    let frames: Vec<Snapshot> = head
        .map(snapshot)
        .chain(rest.into_iter().flatten())
        .map_while(|s| s)
        .collect();

    let mut b = Builder::new("LAMMPS dump");
    b.cell = Some(first.cell);
    b.first_step = first.step;
    for &i in &order {
        let (a, pos) = (&first.detail[i], first.pos[i]);
        let kind = b.kind(&a.kind);
        b.push(Record {
            res_seq: a.mol,
            name: &a.kind,
            element: a.element,
            pos,
            kind,
            ..Record::default()
        });
    }
    b.frames = frames;
    if b.len() == 0 {
        bail!("no atoms in LAMMPS dump");
    }
    Ok(b.finish(id, String::new(), source))
}
