//! XYZ, plain and extended.
//!
//! Each frame is a count, a comment, and one line per atom. The extended
//! dialect (ASE, OVITO, QUIP) puts `key=value` pairs in the comment: the cell
//! as `Lattice="..."` and the column layout as `Properties=...`.

use anyhow::{bail, Result};

use super::element_or_number;
use crate::protein::{Builder, Cell, Record, Snapshot, Source, Structure};

/// `key=value` pairs from an extended XYZ comment line, with quoted values.
fn pairs(comment: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut chars = comment.chars().peekable();
    loop {
        while chars.peek().is_some_and(|c| c.is_whitespace()) {
            chars.next();
        }
        let key: String = chars
            .by_ref()
            .take_while(|c| *c != '=' && !c.is_whitespace())
            .collect();
        if key.is_empty() {
            if chars.peek().is_none() {
                break;
            }
            continue;
        }
        let value = match chars.peek() {
            Some('"') | Some('\'') => {
                let q = chars.next().unwrap();
                chars.by_ref().take_while(|c| *c != q).collect()
            }
            _ => chars.by_ref().take_while(|c| !c.is_whitespace()).collect(),
        };
        out.push((key.to_ascii_lowercase(), value));
        if chars.peek().is_none() {
            break;
        }
    }
    out
}

/// Where the species and position columns are, from `Properties=`.
/// "species:S:1:pos:R:3:forces:R:3" puts species in column 0 and positions in
/// 1 to 3.
fn layout(properties: &str) -> (usize, usize) {
    let f: Vec<&str> = properties.split(':').collect();
    let mut col = 0;
    let mut species = 0;
    let mut pos = 1;
    for chunk in f.chunks(3) {
        let [name, _, width] = chunk else { break };
        let width: usize = width.parse().unwrap_or(1);
        match name.to_ascii_lowercase().as_str() {
            "species" | "element" | "type_name" => species = col,
            "pos" | "positions" | "position" => pos = col,
            _ => {}
        }
        col += width;
    }
    (species, pos)
}

fn lattice(value: &str) -> Option<Cell> {
    let v: Vec<f32> = value
        .split_whitespace()
        .filter_map(|x| x.parse().ok())
        .collect();
    (v.len() == 9).then(|| Cell {
        origin: [0.0; 3],
        vectors: [[v[0], v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]],
    })
}

pub fn parse(text: &str, id: String, source: Source) -> Result<Structure> {
    let mut b = Builder::new("XYZ");
    let mut title = String::new();
    let mut lines = text.lines();
    let mut first = true;
    let mut expected = 0usize;

    loop {
        let Some(count) = lines.by_ref().find(|l| !l.trim().is_empty()) else {
            break;
        };
        let Ok(n) = count.trim().parse::<usize>() else {
            if first {
                bail!("expected an atom count, found {count:?}");
            }
            break;
        };
        let comment = lines.next().unwrap_or("");
        let kv = pairs(comment);
        let get = |k: &str| kv.iter().find(|(key, _)| key == k).map(|(_, v)| v.as_str());
        let (species, pos) = get("properties").map(layout).unwrap_or((0, 1));
        let mut cell = get("lattice").and_then(lattice);
        if let (Some(c), Some(o)) = (cell.as_mut(), get("origin")) {
            let o: Vec<f32> = o.split_whitespace().filter_map(|x| x.parse().ok()).collect();
            if o.len() == 3 {
                c.origin = [o[0], o[1], o[2]];
            }
        }
        if first {
            if get("lattice").is_none() && get("properties").is_none() {
                title = comment.trim().to_string();
            }
            b.cell = cell;
        }

        let mut frame = Vec::with_capacity(if first { 0 } else { n });
        for _ in 0..n {
            let Some(line) = lines.next() else { break };
            let f: Vec<&str> = line.split_whitespace().collect();
            let (Some(label), Some(x), Some(y), Some(z)) = (
                f.get(species),
                f.get(pos).and_then(|v| v.parse::<f32>().ok()),
                f.get(pos + 1).and_then(|v| v.parse::<f32>().ok()),
                f.get(pos + 2).and_then(|v| v.parse::<f32>().ok()),
            ) else {
                continue;
            };
            if first {
                let kind = b.kind(label);
                b.push(Record {
                    name: label,
                    element: element_or_number(label),
                    pos: [x, y, z],
                    kind,
                    ..Record::default()
                });
            } else {
                frame.push([x, y, z]);
            }
        }
        if first {
            expected = b.len();
            first = false;
        } else {
            if frame.len() != expected {
                break;
            }
            b.frames.push(Snapshot { step: None, pos: frame, cell });
        }
    }

    if b.len() == 0 {
        bail!("no atoms found in XYZ file");
    }
    Ok(b.finish(id, title, source))
}
