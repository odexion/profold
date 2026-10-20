//! Bonds from geometry, for files that do not list them.
//!
//! Two atoms are bonded when they sit closer than their covalent radii added
//! up, plus a tolerance — the rule VMD and most viewers use. Neighbours are
//! found on a grid, so this is linear in the number of atoms.

use std::collections::HashMap;

use crate::protein::Atom;

/// Slack on top of the summed covalent radii, in ångström.
const TOLERANCE: f32 = 0.4;
/// Closer than this is two atoms on top of each other, not a bond.
const TOO_CLOSE: f32 = 0.4;
/// Past this many atoms the bonds are not worth the wait; nothing that size
/// is looked at as ball and stick.
const LIMIT: usize = 3_000_000;

fn cell_of(p: [f32; 3], size: f32) -> (i32, i32, i32) {
    (
        (p[0] / size).floor() as i32,
        (p[1] / size).floor() as i32,
        (p[2] / size).floor() as i32,
    )
}

fn grid(atoms: &[Atom], size: f32) -> HashMap<(i32, i32, i32), Vec<u32>> {
    let mut grid: HashMap<(i32, i32, i32), Vec<u32>> = HashMap::new();
    for (i, a) in atoms.iter().enumerate() {
        grid.entry(cell_of(a.pos, size)).or_default().push(i as u32);
    }
    grid
}

pub fn infer(atoms: &[Atom]) -> Vec<[u32; 2]> {
    if atoms.len() > LIMIT {
        return Vec::new();
    }
    let largest = atoms
        .iter()
        .map(|a| a.element.covalent())
        .fold(0.0f32, f32::max);
    if largest <= 0.0 {
        return Vec::new();
    }
    let size = 2.0 * largest + TOLERANCE;
    let grid = grid(atoms, size);

    let mut bonds = Vec::new();
    for (i, a) in atoms.iter().enumerate() {
        let ra = a.element.covalent();
        if ra <= 0.0 {
            continue;
        }
        let (cx, cy, cz) = cell_of(a.pos, size);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(cell) = grid.get(&(cx + dx, cy + dy, cz + dz)) else {
                        continue;
                    };
                    for &j in cell {
                        // Each pair once.
                        if j as usize <= i {
                            continue;
                        }
                        let b = &atoms[j as usize];
                        let rb = b.element.covalent();
                        if rb <= 0.0 || (a.element.is_metal() && b.element.is_metal()) {
                            continue;
                        }
                        // Hydrogens bond to one thing; two of them this close
                        // are a water's pair, not a molecule of H2.
                        if a.element == crate::elements::Element::H
                            && b.element == crate::elements::Element::H
                        {
                            continue;
                        }
                        let d = [
                            a.pos[0] - b.pos[0],
                            a.pos[1] - b.pos[1],
                            a.pos[2] - b.pos[2],
                        ];
                        let d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                        let reach = ra + rb + TOLERANCE;
                        if d2 < reach * reach && d2 > TOO_CLOSE * TOO_CLOSE {
                            bonds.push([i as u32, j]);
                        }
                    }
                }
            }
        }
    }
    bonds
}

/// Half the typical spacing between atoms, which is how large an atom of
/// unknown element is drawn. Taken from the density over the bounding box, so
/// it is right for a liquid or a crystal in whatever units the file uses.
pub fn spacing_radius(atoms: &[Atom]) -> f32 {
    let default = crate::elements::Element::UNKNOWN.vdw();
    if atoms.len() < 2 || atoms.iter().all(|a| a.element.is_known()) {
        return default;
    }
    let mut lo = [f32::MAX; 3];
    let mut hi = [f32::MIN; 3];
    for a in atoms {
        for k in 0..3 {
            lo[k] = lo[k].min(a.pos[k]);
            hi[k] = hi[k].max(a.pos[k]);
        }
    }
    // A flat or linear arrangement has no volume; measure what it has.
    let extents: Vec<f32> = (0..3)
        .map(|k| hi[k] - lo[k])
        .filter(|e| *e > 1e-3)
        .collect();
    if extents.is_empty() {
        return default;
    }
    let measure: f32 = extents.iter().product();
    let spacing = (measure / atoms.len() as f32).powf(1.0 / extents.len() as f32);
    if spacing.is_finite() && spacing > 0.0 {
        spacing * 0.5
    } else {
        default
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elements::Element;

    fn atom(element: Element, pos: [f32; 3]) -> Atom {
        Atom {
            pos,
            element,
            b: 0.0,
            residue: 0,
            kind: 0,
            name: [0; 4],
        }
    }

    #[test]
    fn water() {
        let atoms = [
            atom(Element::O, [0.0, 0.0, 0.0]),
            atom(Element::H, [0.96, 0.0, 0.0]),
            atom(Element::H, [-0.24, 0.93, 0.0]),
        ];
        let mut bonds = infer(&atoms);
        bonds.sort();
        assert_eq!(bonds, vec![[0, 1], [0, 2]]);
    }

    #[test]
    fn metals_do_not_bond() {
        let fe = Element::from_symbol("Fe");
        let atoms = [atom(fe, [0.0; 3]), atom(fe, [2.48, 0.0, 0.0])];
        assert!(infer(&atoms).is_empty());
    }
}
