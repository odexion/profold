//! What a trajectory does over time, worked out once and plotted: how far
//! it drifts from its first frame (RMSD, after fitting each frame onto it),
//! how compact it stays (radius of gyration), how much of it water could
//! touch (SASA), how many polar contacts hold it together, how much each
//! residue moves (RMSF), and its helices and strands frame by frame.

use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;

use rayon::prelude::*;

use crate::protein::{Element, Ss, Structure};

/// A value against frame number.
pub struct Series {
    pub name: &'static str,
    pub unit: &'static str,
    /// Frame index, and the value there.
    pub points: Vec<(usize, f32)>,
}

pub struct Analysis {
    pub series: Vec<Series>,
    /// How far each fitted atom wanders from its mean place, in Å: by
    /// residue for a polymer (its alpha carbon), by atom otherwise.
    pub rmsf: Vec<f32>,
    pub rmsf_by_residue: bool,
    /// Secondary structure of each polymer residue, in `ss_residues`' order,
    /// for the frames sampled.
    pub ss: Vec<(usize, Vec<Ss>)>,
    pub ss_residues: Vec<usize>,
}

/// The most frames the costlier measures are taken on, spread evenly.
const SAMPLED: usize = 200;
/// Polar atoms closer than this, in Å, are in contact.
const CONTACT: f32 = 3.5;
/// The radius of a water molecule, rolled over the atoms for SASA.
const PROBE: f32 = 1.4;

/// The atoms frames are fitted on: a polymer's alpha carbons, or every atom.
fn fit_atoms(s: &Structure) -> (Vec<usize>, bool) {
    let ca: Vec<usize> = s
        .residues
        .iter()
        .filter_map(|r| r.ca.map(|a| a as usize))
        .collect();
    if ca.len() >= 3 {
        (ca, true)
    } else {
        ((0..s.atoms.len()).collect(), false)
    }
}

fn positions(s: &Structure, frame: usize) -> Vec<[f32; 3]> {
    match s.frames.get(frame) {
        Some(f) if f.pos.len() == s.atoms.len() => f.pos.clone(),
        _ => s.atoms.iter().map(|a| a.pos).collect(),
    }
}

/// `pos`'s chosen atoms laid over `reference`'s by Kabsch: the fitted
/// places, in the order chosen.
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

/// Every frame laid over the first by the fitted atoms, so the whole turns
/// and drifts no more and only its motions within are left.
pub fn fit_frames(s: &mut Structure) -> anyhow::Result<()> {
    if s.frames.len() < 2 {
        anyhow::bail!("there is one frame");
    }
    let (atoms, _) = fit_atoms(s);
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

/// Work it all out, counting frames into `done` as they are finished, out
/// of `total(s)`. Nothing if `cancel` is set partway.
pub fn analyse(s: &Structure, cancel: &AtomicBool, done: &AtomicUsize) -> Option<Analysis> {
    let n = s.frame_count().max(1);
    let (atoms, by_residue) = fit_atoms(s);
    let first = positions(s, 0);
    let reference: Vec<[f64; 3]> = atoms.iter().map(|&i| first[i].map(f64::from)).collect();
    let stop = || cancel.load(Ordering::Relaxed);

    // RMSD, Rg and the sums for RMSF, frame by frame.
    let per_frame: Vec<(f32, f32, Vec<[f64; 3]>)> = (0..n)
        .into_par_iter()
        .map(|t| {
            if stop() {
                return (0.0, 0.0, Vec::new());
            }
            let pos = positions(s, t);
            let fit = fitted(&pos, &reference, &atoms);
            let msd = fit
                .iter()
                .zip(&reference)
                .map(|(a, b)| (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>())
                .sum::<f64>()
                / atoms.len().max(1) as f64;
            let c: [f64; 3] = std::array::from_fn(|k| {
                pos.iter().map(|p| f64::from(p[k])).sum::<f64>() / pos.len().max(1) as f64
            });
            let rg = (pos
                .iter()
                .map(|p| {
                    (0..3)
                        .map(|k| (f64::from(p[k]) - c[k]).powi(2))
                        .sum::<f64>()
                })
                .sum::<f64>()
                / pos.len().max(1) as f64)
                .sqrt();
            done.fetch_add(1, Ordering::Relaxed);
            (msd.sqrt() as f32, rg as f32, fit)
        })
        .collect();
    if stop() {
        return None;
    }
    let mut mean = vec![[0.0f64; 3]; atoms.len()];
    for (_, _, fit) in &per_frame {
        for (m, p) in mean.iter_mut().zip(fit) {
            for k in 0..3 {
                m[k] += p[k] / n as f64;
            }
        }
    }
    let mut spread = vec![0.0f64; atoms.len()];
    for (_, _, fit) in &per_frame {
        for ((v, m), p) in spread.iter_mut().zip(&mean).zip(fit) {
            *v += (0..3).map(|k| (p[k] - m[k]).powi(2)).sum::<f64>() / n as f64;
        }
    }
    let rmsf: Vec<f32> = spread.iter().map(|v| v.sqrt() as f32).collect();
    let mut series = vec![
        Series {
            name: "RMSD",
            unit: "Å",
            points: per_frame
                .iter()
                .enumerate()
                .map(|(t, f)| (t, f.0))
                .collect(),
        },
        Series {
            name: "Rg",
            unit: "Å",
            points: per_frame
                .iter()
                .enumerate()
                .map(|(t, f)| (t, f.1))
                .collect(),
        },
    ];
    drop(per_frame);

    // The costlier ones on a sample of the frames.
    let stride = n.div_ceil(SAMPLED).max(1);
    let sampled: Vec<usize> = (0..n).step_by(stride).collect();
    let polar: Vec<usize> = (0..s.atoms.len())
        .filter(|&i| matches!(s.atoms[i].element, Element::N | Element::O))
        .collect();
    if !polar.is_empty() && s.residues.len() > 1 {
        let points = sampled
            .par_iter()
            .map(|&t| (t, contacts(s, &positions(s, t), &polar) as f32))
            .collect();
        series.push(Series {
            name: "Polar contacts",
            unit: "",
            points,
        });
    }
    if stop() {
        return None;
    }
    if s.atoms.iter().any(|a| a.element.is_known()) && s.atoms.len() <= 50_000 {
        let points = sampled
            .par_iter()
            .map(|&t| (t, sasa(s, &positions(s, t))))
            .collect();
        series.push(Series {
            name: "SASA",
            unit: "Å²",
            points,
        });
    }
    if stop() {
        return None;
    }

    // Helices and strands, from each sampled frame's own backbone.
    let ss_residues: Vec<usize> = (0..s.residues.len())
        .filter(|&r| s.residues[r].ca.is_some())
        .collect();
    let ss = if s.has_polymer() && ss_residues.len() >= 5 {
        sampled
            .par_iter()
            .map(|&t| {
                let mut c = s.clone();
                c.set_frame(t);
                (t, ss_residues.iter().map(|&r| c.residues[r].ss).collect())
            })
            .collect()
    } else {
        Vec::new()
    };
    done.store(total(s), Ordering::Relaxed);
    Some(Analysis {
        series,
        rmsf,
        rmsf_by_residue: by_residue,
        ss,
        ss_residues,
    })
}

/// What `done` counts up to.
pub fn total(s: &Structure) -> usize {
    s.frame_count().max(1)
}

/// Nitrogens and oxygens of different residues closer than `CONTACT`.
fn contacts(s: &Structure, pos: &[[f32; 3]], polar: &[usize]) -> usize {
    let key = |p: [f32; 3]| p.map(|x| (x / CONTACT).floor() as i32);
    let mut grid: HashMap<[i32; 3], Vec<usize>> = HashMap::new();
    for &i in polar {
        grid.entry(key(pos[i])).or_default().push(i);
    }
    let mut count = 0;
    for &i in polar {
        let c = key(pos[i]);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(list) = grid.get(&[c[0] + dx, c[1] + dy, c[2] + dz]) else {
                        continue;
                    };
                    for &j in list {
                        if j <= i || s.atoms[i].residue == s.atoms[j].residue {
                            continue;
                        }
                        let d2 = (0..3).map(|k| (pos[i][k] - pos[j][k]).powi(2)).sum::<f32>();
                        if d2 < CONTACT * CONTACT {
                            count += 1;
                        }
                    }
                }
            }
        }
    }
    count
}

/// Solvent-accessible surface area by Shrake and Rupley: points over each
/// atom's sphere grown by a water's radius, counted where no other such
/// sphere covers them.
pub fn sasa(s: &Structure, pos: &[[f32; 3]]) -> f32 {
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
    let cell = radius.iter().cloned().fold(0.0, f32::max) * 2.0;
    let key = |p: [f32; 3]| p.map(|x| (x / cell).floor() as i32);
    let mut grid: HashMap<[i32; 3], Vec<usize>> = HashMap::new();
    for (i, p) in pos.iter().enumerate() {
        grid.entry(key(*p)).or_default().push(i);
    }
    (0..pos.len())
        .into_par_iter()
        .map(|i| {
            let c = key(pos[i]);
            let near: Vec<usize> = (-1..=1)
                .flat_map(|dx| {
                    (-1..=1).flat_map(move |dy| {
                        (-1..=1).map(move |dz| [c[0] + dx, c[1] + dy, c[2] + dz])
                    })
                })
                .filter_map(|k| grid.get(&k))
                .flatten()
                .copied()
                .filter(|&j| j != i)
                .filter(|&j| {
                    (0..3).map(|k| (pos[i][k] - pos[j][k]).powi(2)).sum::<f32>()
                        < (radius[i] + radius[j]).powi(2)
                })
                .collect();
            let open = dots
                .iter()
                .filter(|d| {
                    let p: [f32; 3] = std::array::from_fn(|k| pos[i][k] + d[k] * radius[i]);
                    near.iter().all(|&j| {
                        (0..3).map(|k| (p[k] - pos[j][k]).powi(2)).sum::<f32>()
                            >= radius[j] * radius[j]
                    })
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

    #[test]
    fn a_lone_atom_exposes_its_whole_grown_sphere() {
        let s = crate::formats::load_bytes("one.xyz", b"1\n\nC 0 0 0\n".to_vec()).unwrap();
        let r = Element::from_symbol("C").vdw() + PROBE;
        let area = sasa(&s, &[[0.0; 3]]);
        assert!(
            (area - 4.0 * std::f32::consts::PI * r * r).abs() < 0.01 * area,
            "{area}"
        );
    }

    #[test]
    fn a_turned_and_shifted_copy_has_no_rmsd_and_fitting_undoes_the_turn() {
        let mut s = crambin();
        let first: Vec<[f32; 3]> = s.atoms.iter().map(|a| a.pos).collect();
        let (c, sn) = (0.6f32.cos(), 0.6f32.sin());
        let moved: Vec<[f32; 3]> = first
            .iter()
            .map(|p| [c * p[0] - sn * p[1] + 5.0, sn * p[0] + c * p[1], p[2] - 3.0])
            .collect();
        s.frames = Arc::new(vec![
            crate::protein::Snapshot {
                step: None,
                pos: first.clone(),
                cell: None,
            },
            crate::protein::Snapshot {
                step: None,
                pos: moved,
                cell: None,
            },
        ]);
        let a = analyse(&s, &AtomicBool::new(false), &AtomicUsize::new(0)).unwrap();
        let rmsd = &a.series[0];
        assert_eq!(rmsd.name, "RMSD");
        assert!(rmsd.points[1].1 < 1e-3, "{}", rmsd.points[1].1);
        // Turning does not change how compact it is.
        assert!((a.series[1].points[0].1 - a.series[1].points[1].1).abs() < 1e-3);
        assert!(a.rmsf_by_residue && a.rmsf.iter().all(|v| *v < 1e-3));
        assert!(
            a.series.iter().any(|s| s.name == "SASA")
                && a.series.iter().any(|s| s.name == "Polar contacts")
        );
        assert_eq!(a.ss.len(), 2);
        fit_frames(&mut s).unwrap();
        let d = (0..3)
            .map(|k| (s.frames[1].pos[10][k] - first[10][k]).powi(2))
            .sum::<f32>()
            .sqrt();
        assert!(d < 1e-3, "{d}");
    }
}
