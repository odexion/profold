//! Volumetric data: a value on every point of a grid, as a Gaussian cube's
//! orbital or density, a cryo-EM or crystallographic map (CCP4, MRC), or an
//! electrostatic potential (OpenDX). Drawn as the surface where it crosses a
//! level, and for data with a sign, as a second surface where it crosses
//! the level's negative, the two lobes of an orbital.

use anyhow::{anyhow, bail, Context, Result};
use rayon::prelude::*;

use crate::surface::{march, Facet, CORNER, TETS};

type V = [f32; 3];

/// A grid of values. Point `(i, j, k)` is at `origin + i·axes[0] +
/// j·axes[1] + k·axes[2]`, and its value at `(k·dims[1] + j)·dims[0] + i`:
/// the first axis runs fastest.
#[derive(Clone, Debug)]
pub struct Volume {
    pub origin: V,
    pub axes: [V; 3],
    pub dims: [usize; 3],
    pub data: Vec<f32>,
}

impl Volume {
    pub fn at(&self, i: usize, j: usize, k: usize) -> f32 {
        self.data[(k * self.dims[1] + j) * self.dims[0] + i]
    }

    pub fn point(&self, g: [usize; 3]) -> V {
        let a = &self.axes;
        std::array::from_fn(|c| {
            self.origin[c] + g[0] as f32 * a[0][c] + g[1] as f32 * a[1][c] + g[2] as f32 * a[2][c]
        })
    }

    /// The box the grid fills, as a cell.
    pub fn cell(&self) -> crate::protein::Cell {
        crate::protein::Cell {
            origin: self.origin,
            vectors: std::array::from_fn(|i| self.axes[i].map(|x| x * (self.dims[i].max(2) - 1) as f32)),
        }
    }

    /// Whether the values run both ways from zero, as an orbital's do,
    /// rather than all one way, as a density's.
    pub fn signed(&self) -> bool {
        let (lo, hi) = self.data.par_iter().fold(|| (0f32, 0f32), |(a, b), &x| (a.min(x), b.max(x))).reduce(
            || (0.0, 0.0),
            |(a, b), (c, d)| (a.min(c), b.max(d)),
        );
        lo < 0.0 && hi > 0.0 && lo.abs().min(hi) > 0.1 * lo.abs().max(hi)
    }

    /// A level to draw at first: the mean plus two and a half standard
    /// deviations, as a map is usually contoured; for signed data, of the
    /// values' sizes.
    pub fn default_level(&self) -> f32 {
        let signed = self.signed();
        let n = self.data.len().max(1) as f64;
        let value = |x: f32| if signed { x.abs() as f64 } else { x as f64 };
        let mean = self.data.par_iter().map(|&x| value(x)).sum::<f64>() / n;
        let var = self.data.par_iter().map(|&x| (value(x) - mean).powi(2)).sum::<f64>() / n;
        let level = (mean + 2.5 * var.sqrt()) as f32;
        if level.is_finite() && level != 0.0 { level } else { 1.0 }
    }

    /// The surface where the values cross `level`, facing out of where they
    /// are above it; with `below`, where they fall under it instead.
    pub fn isosurface(&self, level: f32, below: bool) -> Vec<Facet> {
        let [nx, ny, nz] = self.dims;
        if nx < 2 || ny < 2 || nz < 2 {
            return Vec::new();
        }
        let sign = if below { -1.0 } else { 1.0 };
        let value = |i: usize, j: usize, k: usize| sign * (self.at(i, j, k) - level);
        // Normals are found on the grid's own axes, then turned into space by
        // the inverse transpose of the axes, so a skewed grid shades true.
        let m: [[f32; 3]; 3] = std::array::from_fn(|r| std::array::from_fn(|c| self.axes[c][r]));
        let inv = invert(m).unwrap_or([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        let gradient = |x: usize, y: usize, z: usize| {
            let d = |axis: usize| -> f32 {
                let n = self.dims[axis];
                let i = [x, y, z][axis];
                let (a, b) = (i.saturating_sub(1), (i + 1).min(n - 1));
                let mut p = [x, y, z];
                let mut q = [x, y, z];
                p[axis] = a;
                q[axis] = b;
                value(p[0], p[1], p[2]) - value(q[0], q[1], q[2])
            };
            let g = [d(0), d(1), d(2)];
            // inv^T · g
            std::array::from_fn::<f32, 3, _>(|r| inv[0][r] * g[0] + inv[1][r] * g[1] + inv[2][r] * g[2])
        };
        let world = |g: [usize; 3]| self.point(g);
        (0..nz - 1)
            .into_par_iter()
            .flat_map_iter(|z| {
                let mut out = Vec::new();
                for y in 0..ny - 1 {
                    for x in 0..nx - 1 {
                        let cell: [[usize; 3]; 8] = CORNER.map(|c| [x + c[0], y + c[1], z + c[2]]);
                        let v: [f32; 8] = cell.map(|c| value(c[0], c[1], c[2]));
                        if v.iter().all(|v| *v <= 0.0) || v.iter().all(|v| *v > 0.0) {
                            continue;
                        }
                        for tet in TETS {
                            march(&tet, &cell, &v, &world, &gradient, &mut out);
                        }
                    }
                }
                out
            })
            .collect()
    }
}

fn invert(m: [[f32; 3]; 3]) -> Option<[[f32; 3]; 3]> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-20 {
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

/// Whether a file's name says it is a map with no atoms of its own.
pub fn is_map(name: &str) -> bool {
    let lower = name.to_ascii_lowercase();
    let lower = lower.strip_suffix(".gz").unwrap_or(&lower);
    [".mrc", ".map", ".ccp4", ".dx", ".mrcs"].iter().any(|e| lower.ends_with(e))
}

/// A map from its file's bytes, gzipped or not.
pub fn read_map(name: &str, bytes: Vec<u8>) -> Result<Volume> {
    let bytes = crate::formats::gunzip(bytes)?;
    let lower = name.to_ascii_lowercase();
    if lower.trim_end_matches(".gz").ends_with(".dx") {
        let text = String::from_utf8_lossy(&bytes);
        dx(&text).with_context(|| format!("reading {name} as OpenDX"))
    } else {
        mrc(&bytes).with_context(|| format!("reading {name} as a CCP4/MRC map"))
    }
}

/// A CCP4 or MRC map: a 1024-byte header, any extended header, then the
/// values section by section, row by row.
pub fn mrc(b: &[u8]) -> Result<Volume> {
    if b.len() < 1024 {
        bail!("too short for a map header");
    }
    // The machine stamp at 212 says the byte order; little-endian unless it
    // says otherwise, as nearly every map is.
    let big = b[212] == 0x11;
    let int = |o: usize| {
        let w = [b[o], b[o + 1], b[o + 2], b[o + 3]];
        if big { i32::from_be_bytes(w) } else { i32::from_le_bytes(w) }
    };
    let float = |o: usize| f32::from_bits(int(o) as u32);
    let (nc, nr, ns) = (int(0), int(4), int(8));
    let mode = int(12);
    if nc <= 0 || nr <= 0 || ns <= 0 {
        bail!("a grid of {nc} × {nr} × {ns}");
    }
    let dims = [nc as usize, nr as usize, ns as usize];
    let start = [int(16), int(20), int(24)];
    let m = [int(28), int(32), int(36)].map(|x| x.max(1) as f32);
    let cell = [float(40), float(44), float(48)];
    let angles = [float(52), float(56), float(60)].map(|a| if a > 0.0 { a } else { 90.0 });
    let order = [int(64), int(68), int(72)];
    let order = if order.iter().all(|o| (1..=3).contains(o)) { order } else { [1, 2, 3] };
    let extended = int(92).max(0) as usize;
    let origin = [float(196), float(200), float(204)];

    // The cell's edges from its lengths and angles, a along x, b in the xy
    // plane; one grid step along each is the edge over its sampling.
    // A right angle exactly, as nearly every map's are, so its axes are square.
    let trig = |deg: f32| if (deg - 90.0).abs() < 1e-4 { (0.0, 1.0) } else { (deg.to_radians().cos(), deg.to_radians().sin()) };
    let [(cal, _), (cbe, _), (cga, sga)] = angles.map(trig);
    let a = [cell[0], 0.0, 0.0];
    let bv = [cell[1] * cga, cell[1] * sga, 0.0];
    let cx = cell[2] * cbe;
    let cy = cell[2] * (cal - cbe * cga) / sga.max(1e-6);
    let cz = (cell[2] * cell[2] - cx * cx - cy * cy).max(0.0).sqrt();
    let edges = [a, bv, [cx, cy, cz]];
    let step: [V; 3] = std::array::from_fn(|i| edges[i].map(|x| x / m[i]));
    // Columns, rows and sections each run along the cell axis the header names.
    let axes: [V; 3] = std::array::from_fn(|i| step[(order[i] - 1) as usize]);
    let corner = if origin.iter().any(|x| *x != 0.0 && x.is_finite()) {
        origin
    } else {
        let mut o = [0.0f32; 3];
        for i in 0..3 {
            for c in 0..3 {
                o[c] += start[i] as f32 * axes[i][c];
            }
        }
        o
    };

    let count = dims[0] * dims[1] * dims[2];
    let body = b.get(1024 + extended..).ok_or_else(|| anyhow!("the extended header runs past the end"))?;
    let need = |size: usize| -> Result<&[u8]> {
        body.get(..count * size).ok_or_else(|| anyhow!("the map ends before its {count} values"))
    };
    let data: Vec<f32> = match mode {
        0 => need(1)?.iter().map(|&x| x as i8 as f32).collect(),
        1 => need(2)?
            .chunks_exact(2)
            .map(|w| if big { i16::from_be_bytes([w[0], w[1]]) } else { i16::from_le_bytes([w[0], w[1]]) } as f32)
            .collect(),
        2 => need(4)?
            .chunks_exact(4)
            .map(|w| {
                let w = [w[0], w[1], w[2], w[3]];
                if big { f32::from_be_bytes(w) } else { f32::from_le_bytes(w) }
            })
            .collect(),
        6 => need(2)?
            .chunks_exact(2)
            .map(|w| if big { u16::from_be_bytes([w[0], w[1]]) } else { u16::from_le_bytes([w[0], w[1]]) } as f32)
            .collect(),
        other => bail!("values of mode {other}, which is not read"),
    };
    Ok(Volume { origin: corner, axes, dims, data })
}

/// An OpenDX grid, as APBS writes a potential: the counts, origin and
/// steps, then the values with the last axis running fastest.
pub fn dx(text: &str) -> Result<Volume> {
    let mut counts = None;
    let mut origin = None;
    let mut deltas: Vec<V> = Vec::new();
    let mut data: Vec<f32> = Vec::new();
    let mut reading = false;
    let numbers = |s: &str| s.split_whitespace().filter_map(|x| x.parse::<f32>().ok()).collect::<Vec<f32>>();
    for line in text.lines() {
        let line = line.trim();
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        if reading {
            if line.starts_with("attribute") || line.starts_with("object") || line.starts_with("component") {
                reading = false;
                continue;
            }
            data.extend(numbers(line));
            continue;
        }
        if line.starts_with("object") && line.contains("gridpositions") {
            let v: Vec<usize> = line
                .split("counts")
                .nth(1)
                .map(|r| r.split_whitespace().filter_map(|x| x.parse().ok()).collect())
                .unwrap_or_default();
            if v.len() == 3 {
                counts = Some([v[0], v[1], v[2]]);
            }
        } else if let Some(rest) = line.strip_prefix("origin") {
            let v = numbers(rest);
            if v.len() == 3 {
                origin = Some([v[0], v[1], v[2]]);
            }
        } else if let Some(rest) = line.strip_prefix("delta") {
            let v = numbers(rest);
            if v.len() == 3 {
                deltas.push([v[0], v[1], v[2]]);
            }
        } else if line.starts_with("object") && line.contains("data follows") {
            reading = true;
        }
    }
    let counts = counts.ok_or_else(|| anyhow!("no grid counts"))?;
    if deltas.len() != 3 {
        bail!("{} grid steps where three were expected", deltas.len());
    }
    let count = counts.iter().product::<usize>();
    if data.len() < count {
        bail!("{} values for a grid of {count}", data.len());
    }
    data.truncate(count);
    // The file's last axis runs fastest, so it is the volume's first.
    Ok(Volume {
        origin: origin.unwrap_or([0.0; 3]),
        axes: [deltas[2], deltas[1], deltas[0]],
        dims: [counts[2], counts[1], counts[0]],
        data,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A ball of radius `r` about the middle of an n-cube grid of spacing `h`,
    /// as values falling from 1 at the centre.
    fn ball(n: usize, h: f32, r: f32) -> Volume {
        let c = (n - 1) as f32 * h / 2.0;
        let mut data = Vec::with_capacity(n * n * n);
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    let d = ((i as f32 * h - c).powi(2) + (j as f32 * h - c).powi(2) + (k as f32 * h - c).powi(2)).sqrt();
                    data.push(1.0 - d / r);
                }
            }
        }
        Volume { origin: [0.0; 3], axes: [[h, 0.0, 0.0], [0.0, h, 0.0], [0.0, 0.0, h]], dims: [n; 3], data }
    }

    #[test]
    fn a_ball_is_drawn_at_its_radius_with_normals_out() {
        let v = ball(24, 0.5, 4.0);
        let c = 23.0 * 0.5 / 2.0;
        let facets = v.isosurface(0.0, false);
        assert!(facets.len() > 200);
        for f in &facets {
            for (p, n) in f.p.iter().zip(&f.n) {
                let d = [p[0] - c, p[1] - c, p[2] - c];
                let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                assert!((len - 4.0).abs() < 0.1, "{len}");
                assert!(d[0] * n[0] + d[1] * n[1] + d[2] * n[2] > 0.8 * len);
            }
        }
        // Drawn from below, the same place faces in.
        let inner = v.isosurface(0.0, true);
        assert_eq!(inner.len(), facets.len());
        assert!(!v.signed() || v.data.iter().any(|x| *x < 0.0));
    }

    #[test]
    fn an_mrc_map_is_read_with_its_origin_and_spacing() {
        let (nx, ny, nz) = (3i32, 2i32, 2i32);
        let mut b = vec![0u8; 1024];
        let put = |b: &mut Vec<u8>, o: usize, v: [u8; 4]| b[o..o + 4].copy_from_slice(&v);
        for (o, v) in [(0, nx), (4, ny), (8, nz), (12, 2), (28, nx), (32, ny), (36, nz), (64, 1), (68, 2), (72, 3)] {
            put(&mut b, o, v.to_le_bytes());
        }
        for (o, v) in [(40, 6.0f32), (44, 4.0), (48, 2.0), (52, 90.0), (56, 90.0), (60, 90.0), (196, 10.0), (200, 20.0), (204, 30.0)] {
            put(&mut b, o, v.to_le_bytes());
        }
        b[212] = 0x44;
        for i in 0..(nx * ny * nz) {
            b.extend((i as f32).to_le_bytes());
        }
        let v = mrc(&b).unwrap();
        assert_eq!(v.dims, [3, 2, 2]);
        assert_eq!(v.axes, [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(v.origin, [10.0, 20.0, 30.0]);
        assert_eq!(v.at(2, 1, 1), 11.0);
        assert_eq!(v.point([2, 1, 1]), [14.0, 22.0, 31.0]);
    }

    #[test]
    fn an_opendx_grid_runs_its_last_axis_fastest() {
        let text = "# APBS\nobject 1 class gridpositions counts 2 2 3\norigin 1 2 3\ndelta 0.5 0 0\ndelta 0 0.5 0\ndelta 0 0 0.5\n\
                    object 2 class gridconnections counts 2 2 3\nobject 3 class array type double rank 0 items 12 data follows\n\
                    0 1 2\n3 4 5\n6 7 8\n9 10 11\nattribute \"dep\" string \"positions\"\n";
        let v = dx(text).unwrap();
        assert_eq!(v.dims, [3, 2, 2]);
        // The file's (x=1, y=0, z=2) is the 8th value along: 1·6 + 0·3 + 2.
        assert_eq!(v.at(2, 0, 1), 8.0);
        assert_eq!(v.point([2, 0, 1]), [1.5, 2.0, 4.0]);
    }
}
