//! Gaussian cube files: a molecule, and an orbital or a density on a grid
//! through it.

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
    let mut steps = [[0.0f32; 3]; 3];
    let mut counts = [0usize; 3];
    // A negative voxel count means the vectors are already in ångström.
    let mut unit = BOHR;
    for ((row, step), count) in vectors.iter_mut().zip(&mut steps).zip(&mut counts) {
        let v = floats(lines.next())?;
        if v.len() < 4 {
            bail!("missing voxel axes");
        }
        if v[0] < 0.0 {
            unit = 1.0;
        }
        let n = v[0].abs();
        *count = n as usize;
        *step = [v[1], v[2], v[3]];
        *row = [v[1] * n, v[2] * n, v[3] * n];
    }
    for row in vectors.iter_mut().chain(&mut steps) {
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
    // A negative atom count says a line of orbital numbers follows the
    // atoms, and that each point holds a value per orbital: the first is
    // drawn.
    let per_point = if head[0] < 0.0 {
        let ids = floats(lines.next())?;
        ids.first().map_or(1, |n| (*n as usize).max(1))
    } else {
        1
    };
    let total = counts.iter().product::<usize>();
    let mut values = Vec::with_capacity(total);
    'read: for line in lines {
        for x in line.split_whitespace() {
            let Ok(x) = x.parse::<f32>() else { continue };
            values.push(x);
            if values.len() >= total * per_point {
                break 'read;
            }
        }
    }
    let mut s = b.finish(id, title, source);
    if total > 0 && values.len() >= total * per_point {
        // Points run with z fastest, so the cube's last axis is the grid's first.
        let data = values.into_iter().step_by(per_point).collect::<Vec<f32>>();
        let origin = [head[1] * unit, head[2] * unit, head[3] * unit];
        s.volume = Some(std::sync::Arc::new(crate::volume::Volume {
            origin,
            axes: [steps[2], steps[1], steps[0]],
            dims: [counts[2], counts[1], counts[0]],
            data,
        }));
    }
    Ok(s)
}
