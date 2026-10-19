//! GROMACS .gro: fixed columns in nanometres, a box line at the end of each
//! frame, and frames simply following one another.

use anyhow::{bail, Result};

use super::guess_element;
use crate::protein::{Builder, Cell, Record, Snapshot, Source, Structure};

const NM: f32 = 10.0;

/// The box line: three numbers for a rectangular box, nine for a triclinic
/// one, in GROMACS's own order (v1x v2y v3z v1y v1z v2x v2z v3x v3y).
fn box_line(line: &str) -> Option<Cell> {
    let v: Vec<f32> = line
        .split_whitespace()
        .filter_map(|x| x.parse::<f32>().ok())
        .map(|x| x * NM)
        .collect();
    let vectors = match v.len() {
        3 => [[v[0], 0.0, 0.0], [0.0, v[1], 0.0], [0.0, 0.0, v[2]]],
        9 => [[v[0], v[3], v[4]], [v[5], v[1], v[6]], [v[7], v[8], v[2]]],
        _ => return None,
    };
    (vectors.iter().any(|r| r.iter().any(|x| *x != 0.0))).then_some(Cell {
        origin: [0.0; 3],
        vectors,
    })
}

/// Coordinates are written with a fixed width that depends on the precision
/// they were saved at: 8 columns for the usual three decimals. The distance
/// between two decimal points gives it away.
fn coordinate_width(line: &str) -> usize {
    let rest = line.get(20..).unwrap_or("");
    let dots: Vec<usize> = rest.match_indices('.').map(|(i, _)| i).take(2).collect();
    match dots.as_slice() {
        [a, b] if b > a => b - a,
        _ => 8,
    }
}

pub fn parse(text: &str, id: String, source: Source) -> Result<Structure> {
    let mut b = Builder::new("GROMACS");
    let mut lines = text.lines();
    let mut title = String::new();
    let mut first = true;

    while let Some(header) = lines.next() {
        if header.trim().is_empty() {
            continue;
        }
        let Some(n) = lines.next().and_then(|l| l.trim().parse::<usize>().ok()) else {
            if first {
                bail!("expected an atom count on the second line");
            }
            break;
        };
        if first {
            // "Protein in water t= 0.00000 step= 0"; the time is not a name.
            title = header.split(" t=").next().unwrap_or("").trim().to_string();
        }
        let mut frame = Vec::with_capacity(if first { 0 } else { n });
        let mut width = 8;
        for i in 0..n {
            let Some(line) = lines.next() else { break };
            if i == 0 {
                width = coordinate_width(line);
            }
            let col = |a: usize, b: usize| line.get(a..b.min(line.len())).unwrap_or("").trim();
            let at = |k: usize| col(20 + k * width, 20 + (k + 1) * width).parse::<f32>().ok();
            let (Some(x), Some(y), Some(z)) = (at(0), at(1), at(2)) else {
                continue;
            };
            let pos = [x * NM, y * NM, z * NM];
            if first {
                let res_name = col(5, 10);
                let name = col(10, 15);
                b.push(Record {
                    res_name,
                    res_seq: col(0, 5).parse().unwrap_or(0),
                    name,
                    element: guess_element(name, res_name),
                    pos,
                    ..Record::default()
                });
            } else {
                frame.push(pos);
            }
        }
        let cell = lines.next().and_then(box_line);
        if first {
            b.cell = cell;
            first = false;
        } else {
            b.frames.push(Snapshot { step: None, pos: frame, cell });
        }
    }

    if b.len() == 0 {
        bail!("no atoms found in GROMACS file");
    }
    Ok(b.finish(id, title, source))
}
