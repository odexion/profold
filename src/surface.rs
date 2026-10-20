//! A molecular surface: the level set of a sum of Gaussians, one per atom,
//! as VMD's QuickSurf and Mol* draw it. Atoms on their own stand at their
//! radius; neighbours flow together into one smooth skin, much like the
//! surface a probe of water would roll over.
//!
//! The density is sampled on a grid and the surface pulled out of it with
//! marching tetrahedra, which needs no large case table and leaves no holes.

use rayon::prelude::*;

/// One triangle of the surface, with a normal per corner and the atom whose
/// patch it lies in.
pub struct Facet {
    pub p: [[f32; 3]; 3],
    pub n: [[f32; 3]; 3],
    pub atom: u32,
}

/// How sharply each atom's Gaussian falls away: larger is tighter round the
/// atom, smaller blobbier.
const SHARPNESS: f32 = 2.0;
/// Each atom's Gaussian is widened by this much, in ångström, so the skin
/// stands off the atoms the way a water probe would and fills the crevices
/// between them.
pub const SWELL: f32 = 0.4;
/// The density the surface is drawn at. An atom alone has exactly this at
/// its swollen radius.
const LEVEL: f32 = 1.0;
/// Past this many radii an atom adds nothing worth counting.
const REACH: f32 = 2.2;
/// The most grid points sampled, which sets the spacing on big structures.
const MAX_POINTS: f32 = 6.0e6;

/// The surface over atoms at `points` with radii `radii`, sampled every
/// `spacing` ångström or coarser on a large structure.
pub fn gaussian(points: &[[f32; 3]], radii: &[f32], spacing: f32) -> Vec<Facet> {
    if points.is_empty() {
        return Vec::new();
    }
    let radii: Vec<f32> = radii.iter().map(|r| r.max(0.3) + SWELL).collect();
    let rmax = radii.iter().cloned().fold(0.0, f32::max);
    let pad = rmax * REACH + 1.0;
    let mut lo = [f32::MAX; 3];
    let mut hi = [f32::MIN; 3];
    for p in points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k] - pad);
            hi[k] = hi[k].max(p[k] + pad);
        }
    }
    let volume = (0..3).map(|k| hi[k] - lo[k]).product::<f32>();
    let h = spacing.max((volume / MAX_POINTS).cbrt());
    let dims: [usize; 3] = std::array::from_fn(|k| ((hi[k] - lo[k]) / h).ceil() as usize + 1);
    let [nx, ny, nz] = dims;
    let slab = nx * ny;

    // Each slice of the grid, filled from the atoms that reach it.
    let mut by_slice: Vec<Vec<u32>> = vec![Vec::new(); nz];
    for (i, p) in points.iter().enumerate() {
        let reach = radii[i] * REACH;
        let z0 = (((p[2] - reach - lo[2]) / h).floor().max(0.0)) as usize;
        let z1 = ((((p[2] + reach - lo[2]) / h).ceil()) as usize).min(nz - 1);
        for list in &mut by_slice[z0..=z1] {
            list.push(i as u32);
        }
    }
    let mut density = vec![0.0f32; slab * nz];
    density
        .par_chunks_mut(slab)
        .enumerate()
        .for_each(|(iz, dens)| {
            let z = lo[2] + iz as f32 * h;
            for &i in &by_slice[iz] {
                let (p, r) = (points[i as usize], radii[i as usize]);
                let reach = r * REACH;
                let dz = z - p[2];
                if dz.abs() > reach {
                    continue;
                }
                let span = (reach * reach - dz * dz).sqrt();
                let x0 = (((p[0] - span - lo[0]) / h).floor().max(0.0)) as usize;
                let x1 = ((((p[0] + span - lo[0]) / h).ceil()) as usize).min(nx - 1);
                let y0 = (((p[1] - span - lo[1]) / h).floor().max(0.0)) as usize;
                let y1 = ((((p[1] + span - lo[1]) / h).ceil()) as usize).min(ny - 1);
                let inv = 1.0 / (r * r);
                for iy in y0..=y1 {
                    let dy = lo[1] + iy as f32 * h - p[1];
                    for ix in x0..=x1 {
                        let dx = lo[0] + ix as f32 * h - p[0];
                        let d2 = (dx * dx + dy * dy + dz * dz) * inv;
                        if d2 > REACH * REACH {
                            continue;
                        }
                        dens[iy * nx + ix] += (-SHARPNESS * (d2 - 1.0)).exp();
                    }
                }
            }
        });

    let at = |x: usize, y: usize, z: usize| density[z * slab + y * nx + x];
    // The outward normal is down the density's slope.
    let gradient = |x: usize, y: usize, z: usize| {
        let d = |k: usize, i: usize, n: usize| -> f32 {
            let (a, b) = (i.saturating_sub(1), (i + 1).min(n - 1));
            let (pa, pb) = match k {
                0 => (at(a, y, z), at(b, y, z)),
                1 => (at(x, a, z), at(x, b, z)),
                _ => (at(x, y, a), at(x, y, b)),
            };
            pa - pb
        };
        [d(0, x, nx), d(1, y, ny), d(2, z, nz)]
    };

    let world = |g: [usize; 3]| [0, 1, 2].map(|k| lo[k] + g[k] as f32 * h);
    let mut facets: Vec<Facet> = (0..nz - 1)
        .into_par_iter()
        .flat_map_iter(|z| {
            let mut out = Vec::new();
            for y in 0..ny - 1 {
                for x in 0..nx - 1 {
                    let cell: [[usize; 3]; 8] = CORNER.map(|c| [x + c[0], y + c[1], z + c[2]]);
                    let value: [f32; 8] = cell.map(|c| at(c[0], c[1], c[2]) - LEVEL);
                    if value.iter().all(|v| *v <= 0.0) || value.iter().all(|v| *v > 0.0) {
                        continue;
                    }
                    for tet in TETS {
                        march(&tet, &cell, &value, &world, &gradient, &mut out);
                    }
                }
            }
            out
        })
        .collect();

    // Each triangle belongs to the atom whose Gaussian is strongest at its
    // middle, so colours meet along the triangles rather than the coarser
    // grid. Atoms are found through buckets the size of the widest reach.
    let cell = rmax * REACH;
    let key = |p: [f32; 3]| [0, 1, 2].map(|k| ((p[k] - lo[k]) / cell).floor() as i32);
    let mut buckets: std::collections::HashMap<[i32; 3], Vec<u32>> =
        std::collections::HashMap::new();
    for (i, p) in points.iter().enumerate() {
        buckets.entry(key(*p)).or_default().push(i as u32);
    }
    facets.par_iter_mut().for_each(|f| {
        let m = [0, 1, 2].map(|k| (f.p[0][k] + f.p[1][k] + f.p[2][k]) / 3.0);
        let c = key(m);
        let mut best = (f32::MAX, 0u32);
        for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let Some(list) = buckets.get(&[c[0] + dx, c[1] + dy, c[2] + dz]) else {
                        continue;
                    };
                    for &i in list {
                        let p = points[i as usize];
                        let r = radii[i as usize];
                        let d2 =
                            ((m[0] - p[0]).powi(2) + (m[1] - p[1]).powi(2) + (m[2] - p[2]).powi(2))
                                / (r * r);
                        if d2 < best.0 {
                            best = (d2, i);
                        }
                    }
                }
            }
        }
        f.atom = best.1;
    });
    facets
}

/// The cube's corners, and its six tetrahedra about the 0–6 diagonal,
/// which share faces with their neighbours' so the skin has no cracks.
pub(crate) const CORNER: [[usize; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
    [0, 1, 1],
];
pub(crate) const TETS: [[usize; 4]; 6] = [
    [0, 5, 1, 6],
    [0, 1, 2, 6],
    [0, 2, 3, 6],
    [0, 3, 7, 6],
    [0, 7, 4, 6],
    [0, 4, 5, 6],
];

/// The surface's crossing of one tetrahedron: none, one triangle where one
/// corner is on its own side, or two where the corners split two and two.
/// `world` places a grid point, and `gradient` points out of the surface.
pub(crate) fn march(
    tet: &[usize; 4],
    cell: &[[usize; 3]; 8],
    value: &[f32; 8],
    world: &impl Fn([usize; 3]) -> [f32; 3],
    gradient: &impl Fn(usize, usize, usize) -> [f32; 3],
    out: &mut Vec<Facet>,
) {
    let inside: Vec<usize> = tet.iter().copied().filter(|&c| value[c] > 0.0).collect();
    let outside: Vec<usize> = tet.iter().copied().filter(|&c| value[c] <= 0.0).collect();
    if inside.is_empty() || outside.is_empty() {
        return;
    }
    // Where the surface crosses the edge from an inside corner to an outside one.
    let cross = |a: usize, b: usize| {
        let t = value[a] / (value[a] - value[b]);
        let (pa, pb) = (world(cell[a]), world(cell[b]));
        let (ga, gb) = (
            gradient(cell[a][0], cell[a][1], cell[a][2]),
            gradient(cell[b][0], cell[b][1], cell[b][2]),
        );
        let p = [0, 1, 2].map(|k| pa[k] + (pb[k] - pa[k]) * t);
        let n = unit([0, 1, 2].map(|k| ga[k] + (gb[k] - ga[k]) * t));
        (p, n)
    };
    // Its atom is settled once every triangle is made.
    let mut push = |v: [([f32; 3], [f32; 3]); 3]| {
        out.push(Facet {
            p: v.map(|(p, _)| p),
            n: v.map(|(_, n)| n),
            atom: 0,
        });
    };
    match (inside.len(), outside.len()) {
        (1, 3) => push(
            outside
                .iter()
                .map(|&o| cross(inside[0], o))
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),
        ),
        (3, 1) => push(
            inside
                .iter()
                .map(|&i| cross(i, outside[0]))
                .collect::<Vec<_>>()
                .try_into()
                .unwrap(),
        ),
        _ => {
            let (a, b, c, d) = (
                cross(inside[0], outside[0]),
                cross(inside[0], outside[1]),
                cross(inside[1], outside[1]),
                cross(inside[1], outside[0]),
            );
            push([a, b, c]);
            push([a, c, d]);
        }
    }
}

fn unit(v: [f32; 3]) -> [f32; 3] {
    let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if l > 1e-12 {
        [v[0] / l, v[1] / l, v[2] / l]
    } else {
        [0.0, 0.0, 1.0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_lone_atom_is_a_ball_at_its_radius() {
        let facets = gaussian(&[[1.0, 2.0, 3.0]], &[1.5], 0.25);
        assert!(facets.len() > 100);
        let r = 1.5 + SWELL;
        for f in &facets {
            for (p, n) in f.p.iter().zip(&f.n) {
                let d = [p[0] - 1.0, p[1] - 2.0, p[2] - 3.0];
                let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                assert!((len - r).abs() < 0.05, "a corner {len} from the centre");
                // The normal points out.
                assert!(d[0] * n[0] + d[1] * n[1] + d[2] * n[2] > 0.9 * len);
            }
            assert_eq!(f.atom, 0);
        }
    }

    #[test]
    fn two_atoms_close_together_are_one_skin() {
        let far = gaussian(&[[0.0; 3], [12.0, 0.0, 0.0]], &[1.5, 1.5], 0.4);
        let near = gaussian(&[[0.0; 3], [2.5, 0.0, 0.0]], &[1.5, 1.5], 0.4);
        // Between two near atoms the skin is bridged: nothing crosses the
        // line joining them, where two separate balls would have each a cap.
        let on_axis = |fs: &[Facet]| {
            fs.iter().any(|f| {
                f.p.iter()
                    .any(|p| p[1].abs() < 0.2 && p[2].abs() < 0.2 && p[0] > 0.5 && p[0] < 2.0)
            })
        };
        assert!(!on_axis(&near));
        assert!(far.iter().any(|f| f.atom == 0) && far.iter().any(|f| f.atom == 1));
        // Each ball's triangles are its own atom's.
        assert!(far.iter().all(|f| (f.p[0][0] < 6.0) == (f.atom == 0)));
    }
}
