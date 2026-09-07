//! Secondary structure from backbone geometry.
//!
//! A compact Kabsch–Sander assignment: place the amide hydrogens, score every
//! candidate N–H···O=C pair with the DSSP electrostatic energy, then read
//! helices off the n-turn pattern and sheets off the bridge pattern. Doing it
//! from geometry rather than from `HELIX`/`_struct_conf` records means it works
//! the same for a SimpleFold prediction, which ships no annotation at all.

use std::collections::{HashMap, HashSet};

use crate::protein::{Ss, Structure};

/// 0.084 · 332, the DSSP charge product in kcal·Å/mol.
const Q: f32 = 27.888;
/// An energy below this counts as a hydrogen bond.
const HBOND: f32 = -0.5;
/// Neighbour search cell, comfortably above any bond distance.
const CELL: f32 = 9.0;

fn sub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn norm(v: [f32; 3]) -> f32 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn unit(v: [f32; 3]) -> [f32; 3] {
    let n = norm(v);
    if n < 1e-6 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

fn dist(a: [f32; 3], b: [f32; 3]) -> f32 {
    norm(sub(a, b))
}

pub fn assign(s: &mut Structure) {
    let n = s.residues.len();
    if n < 5 {
        return;
    }

    let pos = |idx: u32| s.atoms[idx as usize].pos;

    // Backbone positions, and whether each residue continues the previous one.
    let mut bb: Vec<Option<([f32; 3], [f32; 3], [f32; 3], [f32; 3])>> = vec![None; n];
    for (i, r) in s.residues.iter().enumerate() {
        if let Some((nn, ca, c, o)) = r.backbone() {
            bb[i] = Some((pos(nn), pos(ca), pos(c), pos(o)));
        }
    }

    let mut follows = vec![false; n]; // residue i continues residue i-1
    for i in 1..n {
        let (Some(prev), Some(cur)) = (bb[i - 1], bb[i]) else {
            continue;
        };
        follows[i] = s.residues[i].chain == s.residues[i - 1].chain
            // A peptide bond is 1.33 A; allow slack, reject a chain break.
            && dist(prev.2, cur.0) < 2.0;
    }

    // Amide hydrogen, one angstrom off N along the previous carbonyl.
    let mut h: Vec<Option<[f32; 3]>> = vec![None; n];
    for i in 1..n {
        if !follows[i] || s.residues[i].name == "PRO" {
            continue;
        }
        let (Some(prev), Some(cur)) = (bb[i - 1], bb[i]) else {
            continue;
        };
        let d = unit(sub(prev.2, prev.3)); // C - O of the preceding residue
        h[i] = Some([cur.0[0] + d[0], cur.0[1] + d[1], cur.0[2] + d[2]]);
    }

    // Bucket the alpha carbons so the pair search stays linear.
    let key = |p: [f32; 3]| {
        (
            (p[0] / CELL).floor() as i32,
            (p[1] / CELL).floor() as i32,
            (p[2] / CELL).floor() as i32,
        )
    };
    let mut grid: HashMap<(i32, i32, i32), Vec<u32>> = HashMap::new();
    for (i, b) in bb.iter().enumerate() {
        if let Some(b) = b {
            grid.entry(key(b.1)).or_default().push(i as u32);
        }
    }

    // `bonds` holds (carbonyl residue, amide residue) pairs. `pairs` keeps
    // every spatially close residue pair, because a bridge partner need not be
    // hydrogen bonded to its opposite number — in the parallel pattern the
    // bonds run to i-1 and i+1 instead.
    let mut bonds: HashSet<(u32, u32)> = HashSet::new();
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    let energy = |co: usize, nh: usize| -> f32 {
        let (Some(a), Some(d), Some(hp)) = (bb[co], bb[nh], h[nh]) else {
            return 0.0;
        };
        let (_, _, c, o) = a;
        let nn = d.0;
        Q * (1.0 / dist(o, nn) + 1.0 / dist(c, hp) - 1.0 / dist(o, hp) - 1.0 / dist(c, nn))
    };

    for (cell, members) in &grid {
        for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let other = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
                    let Some(others) = grid.get(&other) else {
                        continue;
                    };
                    for &i in members {
                        for &j in others {
                            if i >= j {
                                continue;
                            }
                            let (i, j) = (i as usize, j as usize);
                            if j - i < 2 {
                                continue;
                            }
                            if dist(bb[i].unwrap().1, bb[j].unwrap().1) > CELL {
                                continue;
                            }
                            pairs.push((i, j));
                            if h[j].is_some() && energy(i, j) < HBOND {
                                bonds.insert((i as u32, j as u32));
                            }
                            if h[i].is_some() && energy(j, i) < HBOND {
                                bonds.insert((j as u32, i as u32));
                            }
                        }
                    }
                }
            }
        }
    }

    let bonded = |co: usize, nh: usize| bonds.contains(&(co as u32, nh as u32));
    // Are residues a..=b one uninterrupted stretch of chain?
    let run = |a: usize, b: usize| b < n && (a + 1..=b).all(|k| follows[k]);

    let mut ss = vec![Ss::Coil; n];

    // Helices: two consecutive n-turns of the same size.
    for i in 0..n {
        for (turn, span) in [(4usize, 4usize), (3, 3), (5, 5)] {
            if i + turn + 1 >= n || !run(i, i + turn + 1) {
                continue;
            }
            if bonded(i, i + turn) && bonded(i + 1, i + turn + 1) {
                for k in i + 1..=i + span {
                    ss[k] = Ss::Helix;
                }
            }
        }
    }

    // Bridges, in both the parallel and antiparallel patterns.
    let mut bridge = vec![false; n];
    for &(i, j) in &pairs {
        if i == 0 || j == 0 || i + 1 >= n || j + 1 >= n {
            continue;
        }
        if i.abs_diff(j) < 3 || !run(i - 1, i + 1) || !run(j - 1, j + 1) {
            continue;
        }
        let parallel =
            (bonded(i - 1, j) && bonded(j, i + 1)) || (bonded(j - 1, i) && bonded(i, j + 1));
        let antiparallel =
            (bonded(i, j) && bonded(j, i)) || (bonded(i - 1, j + 1) && bonded(j - 1, i + 1));
        if parallel || antiparallel {
            bridge[i] = true;
            bridge[j] = true;
        }
    }

    // A lone bridge is noise; a strand needs a neighbour in the same ladder.
    for i in 0..n {
        if !bridge[i] || ss[i] == Ss::Helix {
            continue;
        }
        let before = i > 0 && bridge[i - 1] && follows[i];
        let after = i + 1 < n && bridge[i + 1] && follows[i + 1];
        if before || after {
            ss[i] = Ss::Sheet;
        }
    }

    tidy(&mut ss, &follows);
    for (r, ss) in s.residues.iter_mut().zip(ss) {
        r.ss = ss;
    }
}

/// Close single-residue gaps and drop stretches too short to mean anything.
fn tidy(ss: &mut [Ss], follows: &[bool]) {
    let n = ss.len();
    for i in 1..n.saturating_sub(1) {
        if ss[i] == Ss::Coil && ss[i - 1] == ss[i + 1] && ss[i - 1] != Ss::Coil {
            ss[i] = ss[i - 1];
        }
    }
    let mut i = 0;
    while i < n {
        if ss[i] == Ss::Coil {
            i += 1;
            continue;
        }
        let mut j = i + 1;
        while j < n && ss[j] == ss[i] && follows[j] {
            j += 1;
        }
        let least = if ss[i] == Ss::Helix { 4 } else { 2 };
        if j - i < least {
            for k in i..j {
                ss[k] = Ss::Coil;
            }
        }
        i = j;
    }
}
