//! Structural superposition.
//!
//! Kabsch by way of Horn's quaternion form: build the 4x4 correlation matrix
//! from the matched alpha carbons and take the eigenvector of its largest
//! eigenvalue. That avoids needing a general SVD, and unlike the determinant
//! fix-up on a rotation matrix it can never fall into a reflection.

use anyhow::{bail, Result};

use crate::protein::Structure;

/// Reference colour, and the prediction laid over it.
pub const REFERENCE_TINT: [f32; 3] = [0.36, 0.46, 0.64];
pub const PREDICTION_TINT: [f32; 3] = [0.88, 0.49, 0.22];

pub struct Fit {
    /// Root-mean-square deviation over the matched alpha carbons, in angstrom.
    pub rmsd: f32,
    pub matched: usize,
}

/// Move `prediction` onto `reference` and record how far apart they end up.
///
/// The prediction was folded from the reference's own sequence, so residue *i*
/// of one is residue *i* of the other and no sequence alignment is needed.
pub fn superimpose(reference: &mut Structure, prediction: &mut Structure) -> Result<Fit> {
    let (ref_res, _) = reference.primary_chain();
    let pred_res: Vec<u32> = prediction
        .residues
        .iter()
        .enumerate()
        .filter(|(_, r)| r.ca.is_some())
        .map(|(i, _)| i as u32)
        .collect();

    let n = ref_res.len().min(pred_res.len());
    if n < 3 {
        bail!("not enough matched residues to superimpose");
    }

    let ca = |s: &Structure, ri: u32| s.atoms[s.residues[ri as usize].ca.unwrap() as usize].pos;
    let q: Vec<[f64; 3]> = ref_res[..n]
        .iter()
        .map(|&ri| ca(reference, ri).map(f64::from))
        .collect();
    let p: Vec<[f64; 3]> = pred_res[..n]
        .iter()
        .map(|&ri| ca(prediction, ri).map(f64::from))
        .collect();

    let cp = centroid(&p);
    let cq = centroid(&q);

    // Correlation of the centred point sets.
    let mut s = [[0.0f64; 3]; 3];
    for (a, b) in p.iter().zip(&q) {
        for i in 0..3 {
            for j in 0..3 {
                s[i][j] += (a[i] - cp[i]) * (b[j] - cq[j]);
            }
        }
    }

    let f = [
        [
            s[0][0] + s[1][1] + s[2][2],
            s[1][2] - s[2][1],
            s[2][0] - s[0][2],
            s[0][1] - s[1][0],
        ],
        [
            s[1][2] - s[2][1],
            s[0][0] - s[1][1] - s[2][2],
            s[0][1] + s[1][0],
            s[0][2] + s[2][0],
        ],
        [
            s[2][0] - s[0][2],
            s[0][1] + s[1][0],
            -s[0][0] + s[1][1] - s[2][2],
            s[1][2] + s[2][1],
        ],
        [
            s[0][1] - s[1][0],
            s[0][2] + s[2][0],
            s[1][2] + s[2][1],
            -s[0][0] - s[1][1] + s[2][2],
        ],
    ];

    let (values, vectors) = jacobi(f);
    let best = (0..4).max_by(|a, b| values[*a].total_cmp(&values[*b])).unwrap();
    let quat = [
        vectors[0][best],
        vectors[1][best],
        vectors[2][best],
        vectors[3][best],
    ];
    let r = rotation(quat);

    // Apply to every atom, not just the matched ones.
    for atom in &mut prediction.atoms {
        let v = [
            f64::from(atom.pos[0]) - cp[0],
            f64::from(atom.pos[1]) - cp[1],
            f64::from(atom.pos[2]) - cp[2],
        ];
        for i in 0..3 {
            atom.pos[i] = (r[i][0] * v[0] + r[i][1] * v[1] + r[i][2] * v[2] + cq[i]) as f32;
        }
    }

    // Per-residue deviation, recorded on both so either can be coloured by it.
    reference.deviation = vec![-1.0; reference.residues.len()];
    prediction.deviation = vec![-1.0; prediction.residues.len()];
    let mut total = 0.0f64;
    for k in 0..n {
        let a = ca(prediction, pred_res[k]);
        let b = ca(reference, ref_res[k]);
        let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
        total += f64::from(d) * f64::from(d);
        reference.deviation[ref_res[k] as usize] = d;
        prediction.deviation[pred_res[k] as usize] = d;
    }

    Ok(Fit {
        rmsd: (total / n as f64).sqrt() as f32,
        matched: n,
    })
}

fn centroid(v: &[[f64; 3]]) -> [f64; 3] {
    let mut c = [0.0; 3];
    for p in v {
        for i in 0..3 {
            c[i] += p[i];
        }
    }
    for x in &mut c {
        *x /= v.len() as f64;
    }
    c
}

fn rotation(q: [f64; 4]) -> [[f64; 3]; 3] {
    let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    let (w, x, y, z) = (q[0] / n, q[1] / n, q[2] / n, q[3] / n);
    [
        [
            w * w + x * x - y * y - z * z,
            2.0 * (x * y - w * z),
            2.0 * (x * z + w * y),
        ],
        [
            2.0 * (x * y + w * z),
            w * w - x * x + y * y - z * z,
            2.0 * (y * z - w * x),
        ],
        [
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            w * w - x * x - y * y + z * z,
        ],
    ]
}

/// Cyclic Jacobi for a symmetric 4x4. Returns the eigenvalues and, as columns,
/// their eigenvectors.
fn jacobi(mut a: [[f64; 4]; 4]) -> ([f64; 4], [[f64; 4]; 4]) {
    let mut v = [[0.0f64; 4]; 4];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }

    for _ in 0..64 {
        let off = (0..4)
            .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum::<f64>();
        if off < 1e-24 {
            break;
        }
        for p in 0..3 {
            for q in p + 1..4 {
                if a[p][q].abs() < 1e-30 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;

                for k in 0..4 {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..4 {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for k in 0..4 {
                    let (vkp, vkq) = (v[k][p], v[k][q]);
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }

    ([a[0][0], a[1][1], a[2][2], a[3][3]], v)
}
