//! XCrySDen XSF and AXSF: a molecule under ATOMS, or a crystal under
//! PRIMVEC and PRIMCOORD, with ANIMSTEPS numbering each block of an animation.

use anyhow::{bail, Result};

use super::element_or_number;
use crate::protein::{Builder, Cell, Element, Record, Snapshot, Source, Structure};

fn floats(line: &str) -> Vec<f32> {
    line.split_whitespace()
        .filter_map(|x| x.parse().ok())
        .collect()
}

pub fn parse(text: &str, id: String, source: Source) -> Result<Structure> {
    let lines: Vec<&str> = text
        .lines()
        .map(|l| l.split('#').next().unwrap_or("").trim())
        .filter(|l| !l.is_empty())
        .collect();

    let mut cell: Option<Cell> = None;
    let mut blocks: Vec<(Vec<(Element, [f32; 3])>, Option<Cell>)> = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        let key = lines[i]
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_ascii_uppercase();
        match key.as_str() {
            "PRIMVEC" => {
                let mut v = [[0.0f32; 3]; 3];
                for (k, row) in v.iter_mut().enumerate() {
                    let f = floats(lines.get(i + 1 + k).unwrap_or(&""));
                    if f.len() >= 3 {
                        *row = [f[0], f[1], f[2]];
                    }
                }
                cell = Some(Cell {
                    origin: [0.0; 3],
                    vectors: v,
                });
                i += 4;
            }
            "PRIMCOORD" => {
                let n = floats(lines.get(i + 1).unwrap_or(&""))
                    .first()
                    .copied()
                    .unwrap_or(0.0) as usize;
                let atoms = lines
                    .iter()
                    .skip(i + 2)
                    .take(n)
                    .filter_map(|l| {
                        let f: Vec<&str> = l.split_whitespace().collect();
                        let p: Vec<f32> = f
                            .iter()
                            .skip(1)
                            .take(3)
                            .filter_map(|x| x.parse().ok())
                            .collect();
                        (p.len() == 3).then(|| (element_or_number(f[0]), [p[0], p[1], p[2]]))
                    })
                    .collect();
                blocks.push((atoms, cell));
                i += 2 + n;
            }
            "ATOMS" => {
                let mut atoms = Vec::new();
                i += 1;
                while i < lines.len() {
                    let f: Vec<&str> = lines[i].split_whitespace().collect();
                    let p: Vec<f32> = f
                        .iter()
                        .skip(1)
                        .take(3)
                        .filter_map(|x| x.parse().ok())
                        .collect();
                    if p.len() != 3 {
                        break;
                    }
                    atoms.push((element_or_number(f[0]), [p[0], p[1], p[2]]));
                    i += 1;
                }
                blocks.push((atoms, cell));
            }
            _ => i += 1,
        }
    }

    let mut blocks = blocks.into_iter();
    let Some((first, first_cell)) = blocks.next() else {
        bail!("no atoms found in XSF file");
    };
    let mut b = Builder::new("XSF");
    b.cell = first_cell;
    for (element, pos) in &first {
        let kind = b.kind(element.symbol());
        b.push(Record {
            name: element.symbol(),
            element: *element,
            pos: *pos,
            kind,
            ..Record::default()
        });
    }
    for (atoms, cell) in blocks {
        if atoms.len() != first.len() {
            break;
        }
        b.frames.push(Snapshot {
            step: None,
            pos: atoms.into_iter().map(|a| a.1).collect(),
            cell,
        });
    }
    if b.len() == 0 {
        bail!("no atoms found in XSF file");
    }
    Ok(b.finish(id, String::new(), source))
}
