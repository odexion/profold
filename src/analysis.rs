//! What a trajectory does over time, plotted: how far it drifts from its
//! first frame (RMSD, after fitting each frame onto it), how compact it
//! stays (radius of gyration), how much each atom moves (RMSF), how much of
//! it water could touch (SASA), how many polar contacts hold it together,
//! and its helices and strands frame by frame. Each is worked out only when
//! asked for, over the atoms chosen.

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;

use rayon::prelude::*;

use crate::protein::{Element, Ss, Structure};

/// What can be measured.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Measure {
    Rmsd,
    Rg,
    Rmsf,
    Sasa,
    Contacts,
    Secondary,
}

impl Measure {
    pub const ALL: [Measure; 6] = [
        Measure::Rmsd,
        Measure::Rg,
        Measure::Rmsf,
        Measure::Sasa,
        Measure::Contacts,
        Measure::Secondary,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Measure::Rmsd => "RMSD",
            Measure::Rg => "Rg",
            Measure::Rmsf => "RMSF",
            Measure::Sasa => "SASA",
            Measure::Contacts => "Polar contacts",
            Measure::Secondary => "Secondary",
        }
    }

    /// The name on its pill, where the full one is long.
    pub fn short(self) -> &'static str {
        match self {
            Measure::Contacts => "Contacts",
            m => m.label(),
        }
    }

    pub fn unit(self) -> &'static str {
        match self {
            Measure::Rmsd | Measure::Rg | Measure::Rmsf => "Å",
            Measure::Sasa => "Å²",
            Measure::Contacts | Measure::Secondary => "",
        }
    }

    /// What it is, said in full on hover.
    pub fn about(self) -> &'static str {
        match self {
            Measure::Rmsd => "How far the chosen atoms are from where they were in the first frame, each frame fitted onto it",
            Measure::Rg => "How spread out the chosen atoms are about their centre, frame by frame",
            Measure::Rmsf => "How far each chosen atom wanders from its mean place, over the frames fitted onto the first",
            Measure::Sasa => "The area of the chosen atoms a water molecule could touch, every atom shading the rest (Shrake–Rupley)",
            Measure::Contacts => "Nitrogens and oxygens among the chosen atoms, of different residues, closer than 3.5 Å",
            Measure::Secondary => "Each residue's helix, strand or coil, frame by frame",
        }
    }

    /// Whether it can be taken of these atoms: helices need a protein,
    /// contacts nitrogens or oxygens, and a fit three atoms.
    pub fn applies(self, s: &Structure, atoms: &[usize]) -> bool {
        match self {
            Measure::Secondary => {
                s.has_polymer() && s.residues.iter().filter(|r| r.ca.is_some()).count() >= 5
            }
            Measure::Contacts => s.residues.len() > 1 && atoms.iter().any(|&i| polar(s, i)),
            Measure::Sasa => !atoms.is_empty() && s.atoms.iter().any(|a| a.element.is_known()),
            Measure::Rmsd | Measure::Rmsf => atoms.len() >= 3,
            Measure::Rg => !atoms.is_empty(),
        }
    }
}

/// What a measure comes to.
pub enum Plot {
    /// A value against frame number: the frame, and the value there.
    Time(Vec<(usize, f32)>),
    /// A value for each chosen atom, in order.
    Along(Vec<f32>),
    /// Helix, strand or coil of each residue in `residues`, for the frames
    /// sampled.
    Secondary {
        frames: Vec<(usize, Vec<Ss>)>,
        residues: Vec<usize>,
    },
}

/// The most frames the costlier measures are taken on, spread evenly.
const SAMPLED: usize = 200;
/// Polar atoms closer than this, in Å, are in contact.
const CONTACT: f32 = 3.5;
/// The radius of a water molecule, rolled over the atoms for SASA.
const PROBE: f32 = 1.4;

fn polar(s: &Structure, i: usize) -> bool {
    matches!(s.atoms[i].element, Element::N | Element::O)
}

/// The atoms frames are fitted on to take out tumbling: a polymer's alpha
/// carbons, or every atom.
fn fit_atoms(s: &Structure) -> Vec<usize> {
    let ca: Vec<usize> = s
        .residues
        .iter()
        .filter_map(|r| r.ca.map(|a| a as usize))
        .collect();
    if ca.len() >= 3 {
        ca
    } else {
        (0..s.atoms.len()).collect()
    }
}

/// A frame's positions, read where they are rather than copied.
fn positions(s: &Structure, frame: usize) -> std::borrow::Cow<'_, [[f32; 3]]> {
    match s.frames.get(frame) {
        Some(f) if f.pos.len() == s.atoms.len() => std::borrow::Cow::Borrowed(&f.pos),
        _ => std::borrow::Cow::Owned(s.atoms.iter().map(|a| a.pos).collect()),
    }
}

/// How far each atom strays from its mean place, gathered a frame at a
/// time by Welford's method and merged between threads by Chan's, so the
/// result is the two-pass one without holding the frames.
struct Spread {
    frames: f64,
    mean: Vec<[f64; 3]>,
    /// Summed squared distance from the mean, over the three axes.
    m2: Vec<f64>,
}

impl Spread {
    fn new(atoms: usize) -> Self {
        Spread {
            frames: 0.0,
            mean: vec![[0.0; 3]; atoms],
            m2: vec![0.0; atoms],
        }
    }

    fn add(&mut self, fit: &[[f64; 3]]) {
        self.frames += 1.0;
        for ((mean, m2), p) in self.mean.iter_mut().zip(&mut self.m2).zip(fit) {
            for k in 0..3 {
                let d = p[k] - mean[k];
                mean[k] += d / self.frames;
                *m2 += d * (p[k] - mean[k]);
            }
        }
    }

    fn merge(mut self, other: Spread) -> Spread {
        if other.frames == 0.0 {
            return self;
        }
        if self.frames == 0.0 {
            return other;
        }
        let n = self.frames + other.frames;
        let w = self.frames * other.frames / n;
        for (((mean, m2), om), om2) in self
            .mean
            .iter_mut()
            .zip(&mut self.m2)
            .zip(&other.mean)
            .zip(&other.m2)
        {
            for k in 0..3 {
                let d = om[k] - mean[k];
                *m2 += d * d * w;
                mean[k] += d * other.frames / n;
            }
            *m2 += om2;
        }
        self.frames = n;
        self
    }

    /// Root-mean-square fluctuation of each atom, in Å.
    fn rmsf(&self) -> Vec<f32> {
        self.m2
            .iter()
            .map(|v| (v / self.frames.max(1.0)).sqrt() as f32)
            .collect()
    }
}

/// `pos`'s chosen atoms laid over `reference` by Kabsch: the fitted places,
/// in the order chosen.
fn fitted(pos: &[[f32; 3]], reference: &[[f64; 3]], atoms: &[usize]) -> Vec<[f64; 3]> {
    let p: Vec<[f64; 3]> = atoms.iter().map(|&i| pos[i].map(f64::from)).collect();
    let (r, cp, cq) = crate::align::kabsch(&p, reference);
    p.iter()
        .map(|v| {
            let d = [v[0] - cp[0], v[1] - cp[1], v[2] - cp[2]];
            std::array::from_fn(|i| r[i][0] * d[0] + r[i][1] * d[1] + r[i][2] * d[2] + cq[i])
        })
        .collect()
}

/// Every frame laid over the first by a polymer's alpha carbons, or every
/// atom, so the whole turns and drifts no more and only its motions within
/// are left.
pub fn fit_frames(s: &mut Structure) -> anyhow::Result<()> {
    if s.frames.len() < 2 {
        anyhow::bail!("there is one frame");
    }
    let atoms = fit_atoms(s);
    let first = positions(s, 0);
    let reference: Vec<[f64; 3]> = atoms.iter().map(|&i| first[i].map(f64::from)).collect();
    let frames = Arc::make_mut(&mut s.frames);
    frames.par_iter_mut().skip(1).for_each(|f| {
        let p: Vec<[f64; 3]> = atoms.iter().map(|&i| f.pos[i].map(f64::from)).collect();
        let (r, cp, cq) = crate::align::kabsch(&p, &reference);
        for v in f.pos.iter_mut() {
            let d = [
                f64::from(v[0]) - cp[0],
                f64::from(v[1]) - cp[1],
                f64::from(v[2]) - cp[2],
            ];
            *v = std::array::from_fn(|i| {
                (r[i][0] * d[0] + r[i][1] * d[1] + r[i][2] * d[2] + cq[i]) as f32
            });
        }
        // A cell turned with the atoms would no longer be the box's; it goes.
        f.cell = None;
    });
    s.cell = None;
    let at = s.frame;
    s.set_frame(at);
    Ok(())
}

/// What `done` counts up to for `m`: one for each frame it looks at.
pub fn total(s: &Structure, m: Measure) -> usize {
    let n = s.frame_count().max(1);
    match m {
        Measure::Sasa | Measure::Contacts | Measure::Secondary => {
            n.div_ceil(n.div_ceil(SAMPLED).max(1))
        }
        _ => n,
    }
}

/// Take measure `m` of `atoms` across the trajectory, counting frames into
/// `done` as they are finished. Nothing if `cancel` is set partway, or if
/// it does not apply.
pub fn measure(
    s: &Structure,
    m: Measure,
    atoms: &[usize],
    cancel: &AtomicBool,
    done: &AtomicUsize,
) -> Option<Plot> {
    if !m.applies(s, atoms) {
        return None;
    }
    let n = s.frame_count().max(1);
    let stop = || cancel.load(Ordering::Relaxed);
    let tick = || {
        done.fetch_add(1, Ordering::Relaxed);
    };
    // The costlier ones on a sample of the frames.
    let stride = n.div_ceil(SAMPLED).max(1);
    let sampled: Vec<usize> = (0..n).step_by(stride).collect();
    let plot = match m {
        Measure::Rmsd | Measure::Rmsf => {
            let first = positions(s, 0);
            let reference: Vec<[f64; 3]> = atoms.iter().map(|&i| first[i].map(f64::from)).collect();
            // Each frame is fitted, used and let go: nothing of the length
            // of the trajectory is held.
            if m == Measure::Rmsd {
                let points: Vec<(usize, f32)> = (0..n)
                    .into_par_iter()
                    .map(|t| {
                        if stop() {
                            return (t, 0.0);
                        }
                        let fit = fitted(&positions(s, t), &reference, atoms);
                        let msd = fit
                            .iter()
                            .zip(&reference)
                            .map(|(a, b)| (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>())
                            .sum::<f64>()
                            / atoms.len() as f64;
                        tick();
                        (t, msd.sqrt() as f32)
                    })
                    .collect();
                Plot::Time(points)
            } else {
                let spread = (0..n)
                    .into_par_iter()
                    .fold(
                        || Spread::new(atoms.len()),
                        |mut acc, t| {
                            if !stop() {
                                acc.add(&fitted(&positions(s, t), &reference, atoms));
                                tick();
                            }
                            acc
                        },
                    )
                    .reduce(|| Spread::new(atoms.len()), Spread::merge);
                Plot::Along(spread.rmsf())
            }
        }
        Measure::Rg => Plot::Time(
            (0..n)
                .into_par_iter()
                .map(|t| {
                    let pos = positions(s, t);
                    let k = atoms.len() as f64;
                    let c: [f64; 3] = std::array::from_fn(|d| {
                        atoms.iter().map(|&i| f64::from(pos[i][d])).sum::<f64>() / k
                    });
                    let rg = (atoms
                        .iter()
                        .map(|&i| {
                            (0..3)
                                .map(|d| (f64::from(pos[i][d]) - c[d]).powi(2))
                                .sum::<f64>()
                        })
                        .sum::<f64>()
                        / k)
                        .sqrt();
                    tick();
                    (t, rg as f32)
                })
                .collect(),
        ),
        Measure::Contacts => {
            let polar: Vec<usize> = atoms.iter().copied().filter(|&i| polar(s, i)).collect();
            Plot::Time(
                sampled
                    .par_iter()
                    .map(|&t| {
                        let c = contacts(s, &positions(s, t), &polar) as f32;
                        tick();
                        (t, c)
                    })
                    .collect(),
            )
        }
        Measure::Sasa => Plot::Time(
            sampled
                .par_iter()
                .map(|&t| {
                    let a = sasa(s, &positions(s, t), atoms);
                    tick();
                    (t, a)
                })
                .collect(),
        ),
        Measure::Secondary => {
            let residues: Vec<usize> = (0..s.residues.len())
                .filter(|&r| s.residues[r].ca.is_some())
                .collect();
            let frames = sampled
                .par_iter()
                .map(|&t| {
                    let mut c = s.clone();
                    c.set_frame(t);
                    tick();
                    (t, residues.iter().map(|&r| c.residues[r].ss).collect())
                })
                .collect();
            Plot::Secondary { frames, residues }
        }
    };
    (!stop()).then_some(plot)
}

/// Nitrogens and oxygens of different residues closer than `CONTACT`.
fn contacts(s: &Structure, pos: &[[f32; 3]], polar: &[usize]) -> usize {
    let at: Vec<[f32; 3]> = polar.iter().map(|&i| pos[i]).collect();
    let cells = CellList::new(&at, CONTACT / 2.0);
    (0..polar.len())
        .into_par_iter()
        .map(|a| {
            cells
                .around(at[a])
                .filter(|&&b| {
                    let (i, j) = (polar[a], polar[b]);
                    j > i
                        && s.atoms[i].residue != s.atoms[j].residue
                        && (0..3).map(|k| (at[a][k] - at[b][k]).powi(2)).sum::<f32>()
                            < CONTACT * CONTACT
                })
                .count()
        })
        .sum()
}

/// Atoms sorted into cubes `h` wide on a dense grid over their box, so
/// every atom within twice `h` of a point is found by reading 25 runs of
/// neighbouring cubes, each a contiguous slice: no hashing, no lists to
/// allocate.
struct CellList {
    lo: [f32; 3],
    h: f32,
    dims: [usize; 3],
    /// Where each cube's atoms start in `atoms`, and one past the last.
    start: Vec<u32>,
    atoms: Vec<usize>,
}

impl CellList {
    fn new(pos: &[[f32; 3]], reach: f32) -> Self {
        let mut lo = [f32::MAX; 3];
        let mut hi = [f32::MIN; 3];
        for p in pos {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        if pos.is_empty() {
            (lo, hi) = ([0.0; 3], [0.0; 3]);
        }
        // Cubes no smaller than the reach, and no more of them than a few
        // per atom, however spread out the atoms are.
        let span =
            |h: f32| -> [usize; 3] { std::array::from_fn(|k| ((hi[k] - lo[k]) / h) as usize + 1) };
        let mut h = reach.max(1e-3);
        while span(h).iter().product::<usize>() > 8 * pos.len().max(64) {
            h *= 1.5;
        }
        let dims = span(h);
        let index = |p: [f32; 3]| -> usize {
            let c: [usize; 3] =
                std::array::from_fn(|k| (((p[k] - lo[k]) / h) as usize).min(dims[k] - 1));
            (c[2] * dims[1] + c[1]) * dims[0] + c[0]
        };
        let mut start = vec![0u32; dims.iter().product::<usize>() + 1];
        for p in pos {
            start[index(*p) + 1] += 1;
        }
        for c in 1..start.len() {
            start[c] += start[c - 1];
        }
        let mut fill = start.clone();
        let mut atoms = vec![0; pos.len()];
        for (i, p) in pos.iter().enumerate() {
            let c = index(*p);
            atoms[fill[c] as usize] = i;
            fill[c] += 1;
        }
        CellList {
            lo,
            h,
            dims,
            start,
            atoms,
        }
    }

    /// Every atom in the cubes two either side of the one `p` is in: all
    /// within `2h` of it, and some further.
    fn around(&self, p: [f32; 3]) -> impl Iterator<Item = &usize> + '_ {
        let c: [i64; 3] = std::array::from_fn(|k| ((p[k] - self.lo[k]) / self.h).floor() as i64);
        let [nx, ny, nz] = self.dims.map(|d| d as i64);
        let (x0, x1) = ((c[0] - 2).max(0), (c[0] + 2).min(nx - 1));
        (c[2] - 2..=c[2] + 2)
            .filter(move |z| (0..nz).contains(z))
            .flat_map(move |z| {
                (c[1] - 2..=c[1] + 2)
                    .filter(move |y| (0..ny).contains(y))
                    .map(move |y| (y, z))
            })
            .flat_map(move |(y, z)| {
                let row = ((z * ny + y) * nx) as usize;
                let (a, b) = if x0 <= x1 {
                    (
                        self.start[row + x0 as usize] as usize,
                        self.start[row + x1 as usize + 1] as usize,
                    )
                } else {
                    (0, 0)
                };
                self.atoms[a..b].iter()
            })
    }
}

/// Solvent-accessible surface area by Shrake and Rupley: points over each
/// atom's sphere grown by a water's radius, counted where no other such
/// sphere covers them. Every atom shades the others; the area is that of
/// the atoms in `of`.
pub fn sasa(s: &Structure, pos: &[[f32; 3]], of: &[usize]) -> f32 {
    const POINTS: usize = 64;
    let dots: Vec<[f32; 3]> = (0..POINTS)
        .map(|i| {
            // Evenly over the sphere along a golden spiral.
            let y = 1.0 - (i as f32 + 0.5) * 2.0 / POINTS as f32;
            let r = (1.0 - y * y).sqrt();
            let phi = i as f32 * std::f32::consts::PI * (3.0 - 5f32.sqrt());
            [r * phi.cos(), y, r * phi.sin()]
        })
        .collect();
    let radius: Vec<f32> = s
        .atoms
        .iter()
        .map(|a| if a.element.is_known() { a.element.vdw() } else { s.unknown_radius } + PROBE)
        .collect();
    let cells = CellList::new(pos, radius.iter().cloned().fold(0.0, f32::max));
    of.par_iter()
        .map(|&i| {
            let near = cells
                .around(pos[i])
                .copied()
                .filter(|&j| j != i)
                .map(|j| {
                    (
                        j,
                        (0..3).map(|k| (pos[i][k] - pos[j][k]).powi(2)).sum::<f32>(),
                    )
                })
                .filter(|&(j, d2)| d2 < (radius[i] + radius[j]).powi(2))
                .map(|(j, d2)| (d2, pos[j], radius[j] * radius[j]))
                .collect::<Vec<_>>();
            // Nearest first, and whichever hid the last dot tried first for
            // the next, as dots side by side are mostly hidden by the same
            // neighbour: the same count, sooner.
            let mut near = near;
            near.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
            let mut last = 0;
            let open = dots
                .iter()
                .filter(|d| {
                    let p: [f32; 3] = std::array::from_fn(|k| pos[i][k] + d[k] * radius[i]);
                    let hides = |&(_, c, r2): &(f32, [f32; 3], f32)| {
                        (0..3).map(|k| (p[k] - c[k]).powi(2)).sum::<f32>() < r2
                    };
                    if near.get(last).is_some_and(hides) {
                        return false;
                    }
                    match near.iter().position(hides) {
                        Some(j) => {
                            last = j;
                            false
                        }
                        None => true,
                    }
                })
                .count();
            4.0 * std::f32::consts::PI * radius[i] * radius[i] * open as f32 / POINTS as f32
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn crambin() -> Structure {
        crate::formats::load_bytes(
            "1crn.pdb.gz",
            include_bytes!("../samples/1crn.pdb.gz").to_vec(),
        )
        .unwrap()
    }

    fn take(s: &Structure, m: Measure, atoms: &[usize]) -> Plot {
        measure(s, m, atoms, &AtomicBool::new(false), &AtomicUsize::new(0)).unwrap()
    }

    fn two_frames(s: &mut Structure, second: Vec<[f32; 3]>) {
        let first: Vec<[f32; 3]> = s.atoms.iter().map(|a| a.pos).collect();
        s.frames = Arc::new(vec![
            crate::protein::Snapshot {
                step: None,
                pos: first,
                cell: None,
            },
            crate::protein::Snapshot {
                step: None,
                pos: second,
                cell: None,
            },
        ]);
    }

    #[test]
    fn a_lone_atom_exposes_its_whole_grown_sphere() {
        let s = crate::formats::load_bytes("one.xyz", b"1\n\nC 0 0 0\n".to_vec()).unwrap();
        let r = Element::from_symbol("C").vdw() + PROBE;
        let area = sasa(&s, &[[0.0; 3]], &[0]);
        assert!(
            (area - 4.0 * std::f32::consts::PI * r * r).abs() < 0.01 * area,
            "{area}"
        );
    }

    #[test]
    fn rg_of_a_square_is_its_half_diagonal() {
        let s = crate::formats::load_bytes(
            "sq.xyz",
            b"4\n\nC 0 0 0\nC 2 0 0\nC 2 2 0\nC 0 2 0\n".to_vec(),
        )
        .unwrap();
        let Plot::Time(rg) = take(&s, Measure::Rg, &[0, 1, 2, 3]) else {
            panic!()
        };
        assert!((rg[0].1 - 2f32.sqrt()).abs() < 1e-5);
        // Of two corners only, half the side.
        let Plot::Time(rg) = take(&s, Measure::Rg, &[0, 1]) else {
            panic!()
        };
        assert!((rg[0].1 - 1.0).abs() < 1e-5);
    }

    #[test]
    fn a_turned_and_shifted_copy_has_no_rmsd_and_fitting_undoes_the_turn() {
        let mut s = crambin();
        let first: Vec<[f32; 3]> = s.atoms.iter().map(|a| a.pos).collect();
        let (c, sn) = (0.6f32.cos(), 0.6f32.sin());
        let moved = first
            .iter()
            .map(|p| [c * p[0] - sn * p[1] + 5.0, sn * p[0] + c * p[1], p[2] - 3.0])
            .collect();
        two_frames(&mut s, moved);
        let ca = fit_atoms(&s);
        let Plot::Time(rmsd) = take(&s, Measure::Rmsd, &ca) else {
            panic!()
        };
        assert!(rmsd[1].1 < 1e-3, "{}", rmsd[1].1);
        let Plot::Along(rmsf) = take(&s, Measure::Rmsf, &ca) else {
            panic!()
        };
        assert_eq!(rmsf.len(), ca.len());
        assert!(rmsf.iter().all(|v| *v < 1e-3));
        // Turning does not change how compact it is.
        let Plot::Time(rg) = take(&s, Measure::Rg, &ca) else {
            panic!()
        };
        assert!((rg[0].1 - rg[1].1).abs() < 1e-3);
        assert!(matches!(
            take(&s, Measure::Secondary, &ca),
            Plot::Secondary { .. }
        ));
        // Alpha carbons alone have no polar atoms to make contacts.
        assert!(!Measure::Contacts.applies(&s, &ca));
        let all: Vec<usize> = (0..s.atoms.len()).collect();
        assert!(Measure::Contacts.applies(&s, &all));
        fit_frames(&mut s).unwrap();
        let d = (0..3)
            .map(|k| (s.frames[1].pos[10][k] - first[10][k]).powi(2))
            .sum::<f32>()
            .sqrt();
        assert!(d < 1e-3, "{d}");
    }

    #[test]
    fn the_spread_gathered_in_pieces_is_the_spread_of_the_whole() {
        let frames: Vec<Vec<[f64; 3]>> = (0..7)
            .map(|t| {
                (0..3)
                    .map(|i| {
                        [
                            (t * t + i) as f64 * 0.37,
                            (t as f64).sin() * 2.0 + 100.0,
                            i as f64 - t as f64 * 0.1,
                        ]
                    })
                    .collect()
            })
            .collect();
        let mut whole = Spread::new(3);
        frames.iter().for_each(|f| whole.add(f));
        let (mut a, mut b) = (Spread::new(3), Spread::new(3));
        frames[..3].iter().for_each(|f| a.add(f));
        frames[3..].iter().for_each(|f| b.add(f));
        let merged = a.merge(b).merge(Spread::new(3));
        // Against the two passes written out.
        let direct: Vec<f32> = (0..3)
            .map(|i| {
                let mean: [f64; 3] =
                    std::array::from_fn(|k| frames.iter().map(|f| f[i][k]).sum::<f64>() / 7.0);
                let v = frames
                    .iter()
                    .map(|f| (0..3).map(|k| (f[i][k] - mean[k]).powi(2)).sum::<f64>())
                    .sum::<f64>()
                    / 7.0;
                v.sqrt() as f32
            })
            .collect();
        for (x, y) in [whole.rmsf(), merged.rmsf()]
            .iter()
            .flat_map(|r| r.iter().zip(&direct))
        {
            assert!((x - y).abs() < 1e-5, "{x} against {y}");
        }
    }

    #[test]
    fn sasa_of_a_protein_is_the_brute_force_count() {
        // Every dot against every other atom, as Shrake and Rupley wrote it.
        let s = crambin();
        let pos: Vec<[f32; 3]> = s.atoms.iter().map(|a| a.pos).collect();
        let all: Vec<usize> = (0..pos.len()).collect();
        const POINTS: usize = 64;
        let r: Vec<f32> = s.atoms.iter().map(|a| a.element.vdw() + PROBE).collect();
        let mut brute = 0.0;
        for i in 0..pos.len() {
            let mut open = 0;
            for n in 0..POINTS {
                let y = 1.0 - (n as f32 + 0.5) * 2.0 / POINTS as f32;
                let rr = (1.0 - y * y).sqrt();
                let phi = n as f32 * std::f32::consts::PI * (3.0 - 5f32.sqrt());
                let d = [rr * phi.cos(), y, rr * phi.sin()];
                let p: [f32; 3] = std::array::from_fn(|k| pos[i][k] + d[k] * r[i]);
                if (0..pos.len()).all(|j| {
                    j == i || (0..3).map(|k| (p[k] - pos[j][k]).powi(2)).sum::<f32>() >= r[j] * r[j]
                }) {
                    open += 1;
                }
            }
            brute += 4.0 * std::f32::consts::PI * r[i] * r[i] * open as f32 / POINTS as f32;
        }
        let fast = sasa(&s, &pos, &all);
        assert!((fast - brute).abs() < 1e-2, "{fast} against {brute}");
    }

    #[test]
    fn contacts_are_the_brute_force_count() {
        let s = crambin();
        let pos: Vec<[f32; 3]> = s.atoms.iter().map(|a| a.pos).collect();
        let polar: Vec<usize> = (0..pos.len()).filter(|&i| polar(&s, i)).collect();
        let mut brute = 0;
        for (n, &i) in polar.iter().enumerate() {
            for &j in &polar[n + 1..] {
                let d2 = (0..3).map(|k| (pos[i][k] - pos[j][k]).powi(2)).sum::<f32>();
                if s.atoms[i].residue != s.atoms[j].residue && d2 < CONTACT * CONTACT {
                    brute += 1;
                }
            }
        }
        assert!(brute > 20);
        assert_eq!(contacts(&s, &pos, &polar), brute);
    }

    #[test]
    fn rmsd_and_rmsf_of_one_atom_moved_are_what_geometry_says() {
        // Five atoms; the last moves 1 Å straight out of the plane of the
        // others, which a fit cannot take away but only share out.
        let mut s = crate::formats::load_bytes(
            "five.xyz",
            b"5\n\nC 0 0 0\nC 6 0 0\nC 6 6 0\nC 0 6 0\nC 3 3 0\n".to_vec(),
        )
        .unwrap();
        let mut second: Vec<[f32; 3]> = s.atoms.iter().map(|a| a.pos).collect();
        second[4][2] += 1.0;
        two_frames(&mut s, second);
        let all = [0, 1, 2, 3, 4];
        // The fit moves everything by the centroid's shift, 0.2 Å up: the
        // four corners are 0.2 Å off and the middle 0.8 Å.
        let Plot::Time(rmsd) = take(&s, Measure::Rmsd, &all) else {
            panic!()
        };
        let want = ((4.0 * 0.04 + 0.64) / 5.0f32).sqrt();
        assert!(
            (rmsd[1].1 - want).abs() < 1e-4,
            "{} against {want}",
            rmsd[1].1
        );
        // Over two frames each atom swings half its shift either side.
        let Plot::Along(rmsf) = take(&s, Measure::Rmsf, &all) else {
            panic!()
        };
        assert!(
            (rmsf[4] - 0.4).abs() < 1e-4 && (rmsf[0] - 0.1).abs() < 1e-4,
            "{rmsf:?}"
        );
    }
}
