//! MDL molfiles and SD files, in both the V2000 fixed-column table and the
//! V3000 extended form. Records after the first become frames when they hold
//! the same atoms, which is what a file of conformers does.

use anyhow::{bail, Result};

use crate::protein::{Builder, Element, Record, Snapshot, Source, Structure};

struct Molecule {
    name: String,
    atoms: Vec<(Element, [f32; 3])>,
    /// One-based atom indices, as the file writes them.
    bonds: Vec<[u32; 2]>,
}

fn field(line: &str, a: usize, b: usize) -> &str {
    line.get(a..b.min(line.len())).unwrap_or("").trim()
}

fn read(record: &str) -> Option<Molecule> {
    let lines: Vec<&str> = record.lines().collect();
    let name = lines.first()?.trim().to_string();
    let counts = lines.get(3)?;
    if counts.contains("V3000") {
        return read_v3000(name, &lines[4..]);
    }
    let natoms: usize = field(counts, 0, 3).parse().ok()?;
    let nbonds: usize = field(counts, 3, 6).parse().unwrap_or(0);
    let mut atoms = Vec::with_capacity(natoms);
    for line in lines.get(4..4 + natoms)? {
        let (Ok(x), Ok(y), Ok(z)) = (
            field(line, 0, 10).parse::<f32>(),
            field(line, 10, 20).parse::<f32>(),
            field(line, 20, 30).parse::<f32>(),
        ) else {
            return None;
        };
        atoms.push((Element::from_symbol(field(line, 31, 34)), [x, y, z]));
    }
    let bonds = lines
        .iter()
        .skip(4 + natoms)
        .take(nbonds)
        .filter_map(|l| Some([field(l, 0, 3).parse().ok()?, field(l, 3, 6).parse().ok()?]))
        .collect();
    Some(Molecule { name, atoms, bonds })
}

fn read_v3000(name: String, lines: &[&str]) -> Option<Molecule> {
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut block = "";
    for line in lines {
        let Some(rest) = line.trim().strip_prefix("M  V30 ") else {
            continue;
        };
        let f: Vec<&str> = rest.split_whitespace().collect();
        match f.as_slice() {
            ["BEGIN", what, ..] => block = what,
            ["END", ..] => block = "",
            [_, sym, x, y, z, ..] if block == "ATOM" => {
                atoms.push((
                    Element::from_symbol(sym),
                    [x.parse().ok()?, y.parse().ok()?, z.parse().ok()?],
                ));
            }
            [_, _, a, b, ..] if block == "BOND" => {
                bonds.push([a.parse().ok()?, b.parse().ok()?]);
            }
            _ => {}
        }
    }
    Some(Molecule { name, atoms, bonds })
}

pub fn parse(text: &str, id: String, source: Source) -> Result<Structure> {
    let text = text.replace("\r\n", "\n");
    let mut records = text
        .split("$$$$")
        .map(|r| r.trim_start_matches('\n'))
        .filter(|r| !r.trim().is_empty());
    let Some(first) = records.next().and_then(read) else {
        bail!("no molecule found in SD file");
    };

    let mut b = Builder::new("SDF");
    for (element, pos) in &first.atoms {
        let kind = b.kind(element.symbol());
        b.push(Record {
            name: element.symbol(),
            element: *element,
            pos: *pos,
            kind,
            ..Record::default()
        });
    }
    let n = first.atoms.len() as u32;
    b.bonds = first
        .bonds
        .iter()
        .filter(|[a, c]| (1..=n).contains(a) && (1..=n).contains(c))
        .map(|[a, c]| [a - 1, c - 1])
        .collect();

    for m in records.filter_map(read) {
        let same = m.atoms.len() == first.atoms.len()
            && m.atoms.iter().zip(&first.atoms).all(|(a, f)| a.0 == f.0);
        if !same {
            break;
        }
        b.frames.push(Snapshot {
            step: None,
            pos: m.atoms.iter().map(|a| a.1).collect(),
            cell: None,
        });
    }

    if b.len() == 0 {
        bail!("no atoms found in SD file");
    }
    Ok(b.finish(id, first.name, source))
}
