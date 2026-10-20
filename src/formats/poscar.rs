//! VASP POSCAR, CONTCAR and XDATCAR.
//!
//! A comment, a scale, three lattice vectors, the species (VASP 5 names them;
//! VASP 4 leaves them to the comment line), their counts, then positions in
//! direct or Cartesian coordinates. XDATCAR repeats the position block once
//! per ionic step.

use anyhow::{anyhow, bail, Result};

use crate::protein::{Builder, Cell, Element, Record, Snapshot, Source, Structure};

fn floats(line: &str) -> Vec<f32> {
    line.split_whitespace()
        .map_while(|x| x.parse().ok())
        .collect()
}

pub fn parse(text: &str, id: String, source: Source) -> Result<Structure> {
    let mut lines = text.lines();
    let comment = lines.next().unwrap_or("").trim().to_string();
    let scale = floats(lines.next().unwrap_or(""));
    let scale = *scale
        .first()
        .ok_or_else(|| anyhow!("missing scale factor"))?;

    let mut vectors = [[0.0f32; 3]; 3];
    for row in &mut vectors {
        let v = floats(lines.next().unwrap_or(""));
        if v.len() < 3 {
            bail!("missing lattice vectors");
        }
        *row = [v[0], v[1], v[2]];
    }
    // A negative scale is the cell volume rather than a factor.
    let scale = if scale < 0.0 {
        let [a, b, c] = vectors;
        let volume = (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            .abs();
        (-scale / volume).cbrt()
    } else {
        scale
    };
    for row in &mut vectors {
        for x in row.iter_mut() {
            *x *= scale;
        }
    }
    let cell = Cell {
        origin: [0.0; 3],
        vectors,
    };

    let mut line = lines.next().unwrap_or("");
    let names: Vec<String> = if floats(line).is_empty() {
        let n = line.split_whitespace().map(|s| s.to_string()).collect();
        line = lines.next().unwrap_or("");
        n
    } else {
        comment.split_whitespace().map(|s| s.to_string()).collect()
    };
    let counts: Vec<usize> = line
        .split_whitespace()
        .map_while(|x| x.parse().ok())
        .collect();
    if counts.is_empty() {
        bail!("missing atom counts");
    }
    let total: usize = counts.iter().sum();
    let species: Vec<(Element, String)> = counts
        .iter()
        .enumerate()
        .flat_map(|(i, n)| {
            // "Fe_pv" and "O_s" are POTCAR names; the element comes first.
            let name = names
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("type {}", i + 1));
            let sym = name.split(['_', '/']).next().unwrap_or("").to_string();
            std::iter::repeat_n((Element::from_symbol(&sym), name), *n)
        })
        .collect();

    let mut b = Builder::new("VASP");
    b.cell = Some(cell);
    let mut first = true;
    loop {
        // "Selective dynamics", then "Direct"/"Cartesian", or in an XDATCAR
        // "Direct configuration= 12".
        let mut mode = loop {
            match lines.next() {
                None => break None,
                Some(l) if l.trim().is_empty() => continue,
                Some(l) => {
                    let c = l.trim().chars().next().unwrap_or(' ').to_ascii_lowercase();
                    if c == 's' {
                        continue;
                    }
                    break Some(c);
                }
            }
        };
        let Some(m) = mode.take() else { break };
        let cartesian = m == 'c' || m == 'k';
        let mut frame = Vec::with_capacity(total);
        for i in 0..total {
            let Some(l) = lines.next() else { break };
            let v = floats(l);
            if v.len() < 3 {
                break;
            }
            let p = if cartesian {
                [v[0] * scale, v[1] * scale, v[2] * scale]
            } else {
                cell.to_cartesian([v[0], v[1], v[2]])
            };
            if first {
                let (element, name) = &species[i];
                let kind = b.kind(name);
                b.push(Record {
                    name,
                    element: *element,
                    pos: p,
                    kind,
                    ..Record::default()
                });
            } else {
                frame.push(p);
            }
        }
        if first {
            first = false;
        } else {
            if frame.len() != total {
                break;
            }
            b.frames.push(Snapshot {
                step: None,
                pos: frame,
                cell: None,
            });
        }
    }

    if b.len() == 0 {
        bail!("no atoms found in POSCAR");
    }
    Ok(b.finish(id, comment, source))
}
