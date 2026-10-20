//! Periodic boundaries: what a simulation box is for. Atoms put back inside
//! the cell, molecules split across its walls made whole again, a
//! trajectory followed through the walls rather than jumping back, and the
//! cell repeated as its neighbours would be. Every frame is changed alike.

use std::sync::Arc;

use anyhow::{bail, Result};

use crate::protein::{Cell, Structure};

type V = [f32; 3];

/// Fractional coordinates in a cell, and back.
struct Frame {
    cell: Cell,
    /// Inverse of the cell's vectors as columns.
    inv: [[f32; 3]; 3],
}

impl Frame {
    fn new(cell: Cell) -> Option<Self> {
        let m: [[f32; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| cell.vectors[j][i]));
        Some(Frame {
            cell,
            inv: inverse(m)?,
        })
    }

    fn fractional(&self, p: V) -> V {
        let d = std::array::from_fn::<f32, 3, _>(|k| p[k] - self.cell.origin[k]);
        std::array::from_fn(|i| {
            self.inv[i][0] * d[0] + self.inv[i][1] * d[1] + self.inv[i][2] * d[2]
        })
    }

    /// `f` cell vectors' worth of displacement.
    fn shift(&self, f: V) -> V {
        let v = &self.cell.vectors;
        std::array::from_fn(|k| f[0] * v[0][k] + f[1] * v[1][k] + f[2] * v[2][k])
    }

    /// The shortest way from `a` to `b` through the walls.
    fn nearest(&self, a: V, b: V) -> V {
        let d: V = std::array::from_fn(|k| b[k] - a[k]);
        let f: V = std::array::from_fn(|i| {
            self.inv[i][0] * d[0] + self.inv[i][1] * d[1] + self.inv[i][2] * d[2]
        });
        let s = self.shift(f.map(|x| x.round()));
        std::array::from_fn(|k| d[k] - s[k])
    }
}

fn inverse(m: [[f32; 3]; 3]) -> Option<[[f32; 3]; 3]> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-12 {
        return None;
    }
    let d = 1.0 / det;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * d,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * d,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * d,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * d,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * d,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * d,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * d,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * d,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * d,
        ],
    ])
}

/// Run `f` over the positions of every frame with that frame's cell, then
/// show the current frame again.
fn each_frame(s: &mut Structure, mut f: impl FnMut(&mut [V], &Frame)) -> Result<()> {
    let Some(base) = s.cell else {
        bail!("there is no cell")
    };
    if s.frames.is_empty() {
        let frame = Frame::new(base).ok_or_else(|| anyhow::anyhow!("the cell is flat"))?;
        let mut pos: Vec<V> = s.atoms.iter().map(|a| a.pos).collect();
        f(&mut pos, &frame);
        for (a, p) in s.atoms.iter_mut().zip(pos) {
            a.pos = p;
        }
        return Ok(());
    }
    for snap in Arc::make_mut(&mut s.frames).iter_mut() {
        if let Some(frame) = Frame::new(snap.cell.unwrap_or(base)) {
            f(&mut snap.pos, &frame);
        }
    }
    let at = s.frame;
    s.set_frame(at);
    Ok(())
}

/// Every atom moved by whole cells to lie inside it.
pub fn wrap(s: &mut Structure) -> Result<()> {
    each_frame(s, |pos, frame| {
        for p in pos.iter_mut() {
            let f = frame.fractional(*p);
            let s = frame.shift(f.map(|x| -x.floor()));
            *p = std::array::from_fn(|k| p[k] + s[k]);
        }
    })
}

/// The runs of atoms that belong together, each in an order where every
/// atom after the first is near one before it: a polymer's chains residue
/// by residue, and otherwise what bonds join.
fn molecules(s: &Structure) -> Vec<Vec<(usize, usize)>> {
    let n = s.atoms.len();
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &[a, b] in &s.bonds {
        neighbours[a as usize].push(b as usize);
        neighbours[b as usize].push(a as usize);
    }
    // A polymer chain's atoms, residue after residue, are each placed after
    // the one before them, as a chain runs.
    let mut placed = vec![false; n];
    let mut out = Vec::new();
    if s.has_polymer() {
        let by_residue = s.atoms_by_residue();
        for chain in &s.chains {
            let polymer = chain
                .residues
                .iter()
                .any(|&r| !s.residues[r as usize].outside_polymer());
            if !polymer {
                continue;
            }
            let atoms: Vec<usize> = chain
                .residues
                .iter()
                .flat_map(|&r| by_residue.of(r as usize).iter().map(|&a| a as usize))
                .collect();
            let run: Vec<(usize, usize)> = atoms.windows(2).map(|w| (w[0], w[1])).collect();
            for &a in &atoms {
                placed[a] = true;
            }
            out.push(run);
        }
    }
    // Everything else by its bonds, breadth first from each unplaced atom.
    for root in 0..n {
        if placed[root] {
            continue;
        }
        placed[root] = true;
        let mut run = Vec::new();
        let mut queue = std::collections::VecDeque::from([root]);
        while let Some(a) = queue.pop_front() {
            for &b in &neighbours[a] {
                if !placed[b] {
                    placed[b] = true;
                    run.push((a, b));
                    queue.push_back(b);
                }
            }
        }
        if !run.is_empty() {
            out.push(run);
        }
    }
    out
}

/// Molecules cut by the walls put back together, each atom placed next to
/// the one it follows, wherever the first of them lies.
pub fn make_whole(s: &mut Structure) -> Result<()> {
    let runs = molecules(s);
    each_frame(s, |pos, frame| {
        for run in &runs {
            for &(from, to) in run {
                let d = frame.nearest(pos[from], pos[to]);
                pos[to] = std::array::from_fn(|k| pos[from][k] + d[k]);
            }
        }
    })
}

/// A trajectory followed through the walls: each atom goes on from where it
/// was in the frame before, by the shortest way, rather than jumping back
/// into the cell when it leaves.
pub fn unwrap_in_time(s: &mut Structure) -> Result<()> {
    if s.frames.len() < 2 {
        bail!("there is one frame");
    }
    let Some(base) = s.cell else {
        bail!("there is no cell")
    };
    let frames = Arc::make_mut(&mut s.frames);
    let mut last_raw = frames[0].pos.clone();
    for t in 1..frames.len() {
        let Some(frame) = Frame::new(frames[t].cell.unwrap_or(base)) else {
            continue;
        };
        let raw = frames[t].pos.clone();
        if raw.len() != last_raw.len() {
            break;
        }
        let before = frames[t - 1].pos.clone();
        for (i, p) in frames[t].pos.iter_mut().enumerate() {
            let d = frame.nearest(last_raw[i], raw[i]);
            *p = std::array::from_fn(|k| before[i][k] + d[k]);
        }
        last_raw = raw;
    }
    let at = s.frame;
    s.set_frame(at);
    Ok(())
}

/// The cell and everything in it repeated `n[0]` by `n[1]` by `n[2]`
/// times, its copies laid out along its own vectors, every frame alike.
/// The cell grows to hold them all.
pub fn replicate(s: &Structure, n: [usize; 3]) -> Result<Structure> {
    let Some(base) = s.cell else {
        bail!("there is no cell")
    };
    if n.iter().any(|&k| k == 0) || n.iter().product::<usize>() > 1000 {
        bail!("between 1 and 1000 copies, in all");
    }
    if s.chains.len() * n.iter().product::<usize>() > u16::MAX as usize {
        bail!("too many chains to copy");
    }
    let images: Vec<[f32; 3]> = (0..n[0])
        .flat_map(|a| {
            (0..n[1]).flat_map(move |b| (0..n[2]).map(move |c| [a as f32, b as f32, c as f32]))
        })
        .collect();
    let count = s.atoms.len();
    let residues = s.residues.len();
    let shift = |cell: &Cell, f: [f32; 3]| -> V {
        let v = &cell.vectors;
        std::array::from_fn(|k| f[0] * v[0][k] + f[1] * v[1][k] + f[2] * v[2][k])
    };
    let grow = |cell: &Cell| Cell {
        origin: cell.origin,
        vectors: std::array::from_fn(|i| cell.vectors[i].map(|x| x * n[i] as f32)),
    };
    let copy_pos = |pos: &[V], cell: &Cell| -> Vec<V> {
        images
            .iter()
            .flat_map(|&f| {
                let d = shift(cell, f);
                pos.iter()
                    .map(move |p| [p[0] + d[0], p[1] + d[1], p[2] + d[2]])
            })
            .collect()
    };

    let mut out = s.clone();
    out.atoms = images
        .iter()
        .enumerate()
        .flat_map(|(k, &f)| {
            let d = shift(&base, f);
            s.atoms.iter().map(move |a| {
                let mut a = *a;
                a.pos = [a.pos[0] + d[0], a.pos[1] + d[1], a.pos[2] + d[2]];
                a.residue += (k * residues) as u32;
                a
            })
        })
        .collect();
    let offset = |x: Option<u32>, k: usize| x.map(|x| x + (k * count) as u32);
    out.residues = (0..images.len())
        .flat_map(|k| {
            s.residues.iter().map(move |r| {
                let mut r = r.clone();
                r.chain += (k * s.chains.len()) as u16;
                (r.ca, r.n, r.c, r.o) = (
                    offset(r.ca, k),
                    offset(r.n, k),
                    offset(r.c, k),
                    offset(r.o, k),
                );
                (r.nucleic, r.base) = (offset(r.nucleic, k), offset(r.base, k));
                r
            })
        })
        .collect();
    out.chains = (0..images.len())
        .flat_map(|k| {
            s.chains.iter().map(move |c| {
                let mut c = c.clone();
                if k > 0 {
                    c.id = format!("{}{}", c.id, k + 1);
                }
                c.residues = c
                    .residues
                    .iter()
                    .map(|r| r + (k * residues) as u32)
                    .collect();
                c
            })
        })
        .collect();
    out.bonds = (0..images.len())
        .flat_map(|k| {
            s.bonds
                .iter()
                .map(move |&[a, b]| [a + (k * count) as u32, b + (k * count) as u32])
        })
        .collect();
    out.frames = Arc::new(
        s.frames
            .iter()
            .map(|f| {
                let cell = f.cell.unwrap_or(base);
                let mut f = f.clone();
                f.pos = copy_pos(&f.pos, &cell);
                f.cell = Some(grow(&cell));
                f
            })
            .collect(),
    );
    out.cell = Some(grow(&base));
    let times = images.len();
    for p in &mut out.props {
        p.values = p.values.repeat(times);
        p.frames = Arc::new(p.frames.iter().map(|v| v.repeat(times)).collect());
    }
    out.forget_derived();
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const DUMP: &str = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 10.0
0.0 10.0
ITEM: ATOMS id type x y z c_pe
1 1 9.5 5.0 5.0 -1
2 1 5.0 5.0 5.0 -2
ITEM: TIMESTEP
1
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 10.0
0.0 10.0
ITEM: ATOMS id type x y z c_pe
1 1 0.5 5.0 5.0 -1
2 1 5.0 5.0 5.0 -2
";

    fn dump() -> Structure {
        crate::formats::load_bytes("t.lammpstrj", DUMP.as_bytes().to_vec()).unwrap()
    }

    fn close(a: V, b: V) -> bool {
        (0..3).all(|k| (a[k] - b[k]).abs() < 1e-4)
    }

    #[test]
    fn an_atom_that_leaves_the_cell_goes_on_rather_than_jumping_back() {
        let mut s = dump();
        unwrap_in_time(&mut s).unwrap();
        assert!(
            close(s.frames[1].pos[0], [10.5, 5.0, 5.0]),
            "{:?}",
            s.frames[1].pos[0]
        );
        // And wrapped, it is back in the cell.
        wrap(&mut s).unwrap();
        assert!(close(s.frames[1].pos[0], [0.5, 5.0, 5.0]));
    }

    #[test]
    fn a_molecule_cut_by_a_wall_is_made_whole() {
        let mut s = dump();
        s.bonds = vec![[0, 1]];
        s.frames = Arc::new(s.frames[..1].to_vec());
        // Atom 2 one cell over from atom 1's neighbour.
        Arc::make_mut(&mut s.frames)[0].pos[1] = [0.2, 5.0, 5.0];
        s.set_frame(0);
        make_whole(&mut s).unwrap();
        assert!(
            close(s.atoms[1].pos, [10.2, 5.0, 5.0]),
            "{:?}",
            s.atoms[1].pos
        );
    }

    #[test]
    fn replicas_fill_a_bigger_cell_and_keep_their_properties() {
        let s = dump();
        let r = replicate(&s, [2, 1, 3]).unwrap();
        assert_eq!(r.atoms.len(), 12);
        assert_eq!(r.cell.unwrap().lengths(), [20.0, 10.0, 30.0]);
        assert_eq!(r.frames.len(), 2);
        assert!(r.frames.iter().all(|f| f.pos.len() == 12));
        assert_eq!(r.props[0].values.len(), 12);
        // The last copy is one cell along x and two along z.
        assert!(
            close(r.atoms[11].pos, [15.0, 5.0, 25.0]),
            "{:?}",
            r.atoms[11].pos
        );
        assert!(r
            .atoms
            .iter()
            .all(|a| (a.residue as usize) < r.residues.len()));
        assert!(replicate(&s, [0, 1, 1]).is_err());
    }
}
