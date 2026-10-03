//! Gaussian cube files. Only the molecule and the grid's box are read; the
//! volumetric data is left for another day.

use anyhow::{anyhow, bail, Result};

use crate::protein::{Builder, Cell, Element, Record, Source, Structure};

const BOHR: f32 = 0.529_177_2;

fn floats(line: Option<&str>) -> Result<Vec<f32>> {
    let line = line.ok_or_else(|| anyhow!("cube file ends early"))?;
    Ok(line
        .split_whitespace()
        .map_while(|x| x.parse().ok())
        .collect())
}

pub fn parse(text: &str, id: String, source: Source) -> Result<Structure> {
    let mut lines = text.lines();
    let title = lines.next().unwrap_or("").trim().to_string();
    lines.next();
    let head = floats(lines.next())?;
    if head.len() < 4 {
        bail!("missing atom count and origin");
    }
    let natoms = head[0].abs() as usize;

    let mut vectors = [[0.0f32; 3]; 3];
    // A negative voxel count means the vectors are already in ångström.
    let mut unit = BOHR;
    for row in &mut vectors {
        let v = floats(lines.next())?;
        if v.len() < 4 {
            bail!("missing voxel axes");
        }
        if v[0] < 0.0 {
            unit = 1.0;
        }
        let n = v[0].abs();
        *row = [v[1] * n, v[2] * n, v[3] * n];
    }
    for row in &mut vectors {
        for x in row.iter_mut() {
            *x *= unit;
        }
    }

    let mut b = Builder::new("Gaussian cube");
    b.cell = Some(Cell {
        origin: [head[1] * unit, head[2] * unit, head[3] * unit],
        vectors,
    });
    for _ in 0..natoms {
        let v = floats(lines.next())?;
        if v.len() < 5 {
            bail!("short atom line");
        }
        let element = Element::from_number(v[0] as u32);
        let kind = b.kind(element.symbol());
        b.push(Record {
            name: element.symbol(),
            element,
            pos: [v[2] * unit, v[3] * unit, v[4] * unit],
            kind,
            ..Record::default()
        });
    }
    Ok(b.finish(id, title, source))
}
