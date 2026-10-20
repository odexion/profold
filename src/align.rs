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

// The fold plugin's own comparison, residue i against residue i.
#[cfg_attr(not(feature = "simplefold"), allow(dead_code))]
pub struct Fit {
    /// Root-mean-square deviation over the matched alpha carbons, in angstrom.
    pub rmsd: f32,
    pub matched: usize,
}

/// Move `prediction` onto `reference` and record how far apart they end up.
///
/// The prediction was folded from the reference's own sequence, so residue *i*
/// of one is residue *i* of the other and no sequence alignment is needed.
#[cfg_attr(not(feature = "simplefold"), allow(dead_code))]
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

    let (r, cp, cq) = kabsch(&p, &q);

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

/// The rotation that best lays points `p` over points `q`, after moving
/// each to its centroid; and the two centroids.
pub(crate) fn kabsch(p: &[[f64; 3]], q: &[[f64; 3]]) -> ([[f64; 3]; 3], [f64; 3], [f64; 3]) {
    let cp = centroid(p);
    let cq = centroid(q);

    // Correlation of the centred point sets.
    let mut s = [[0.0f64; 3]; 3];
    for (a, b) in p.iter().zip(q) {
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
    (r, cp, cq)
}

/// How well two unrelated structures were laid over each other.
pub struct Overlay {
    pub rmsd: f32,
    /// The pairs the fit kept, once outliers were let go.
    pub matched: usize,
    /// The pairs the sequence alignment found.
    pub aligned: usize,
    /// Of those, how many are the same residue.
    pub identical: usize,
}

/// Move `mobile` onto `reference` where their sequences align, as PyMOL's
/// `align` does: the residues are paired by a global alignment of their
/// sequences, the alpha carbons of the pairs fitted, and the pairs that sit
/// furthest apart let go before fitting again, five times or until none
/// is. Every residue paired is given its deviation, kept or not.
pub fn superpose(reference: &mut Structure, mobile: &mut Structure) -> Result<Overlay> {
    let traced = |s: &Structure| -> Vec<usize> {
        (0..s.residues.len()).filter(|&r| s.residues[r].ca.is_some()).collect()
    };
    let (ra, rb) = (traced(reference), traced(mobile));
    let seq = |s: &Structure, r: &[usize]| -> Vec<char> { r.iter().map(|&i| s.residues[i].code).collect() };
    let pairs: Vec<(usize, usize)> = align_sequences(&seq(reference, &ra), &seq(mobile, &rb))
        .into_iter()
        .map(|(i, j)| (ra[i], rb[j]))
        .collect();
    if pairs.len() < 3 {
        bail!("sequence identity is too low for superposition");
    }
    let identical = pairs.iter().filter(|(a, b)| reference.residues[*a].code == mobile.residues[*b].code).count();
    let ca = |s: &Structure, r: usize| s.atoms[s.residues[r].ca.unwrap() as usize].pos.map(f64::from);

    let mut kept = pairs.clone();
    let mut fit = ([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [0.0; 3], [0.0; 3]);
    let moved = |fit: &([[f64; 3]; 3], [f64; 3], [f64; 3]), p: [f64; 3]| -> [f64; 3] {
        let (r, cp, cq) = fit;
        let v = [p[0] - cp[0], p[1] - cp[1], p[2] - cp[2]];
        std::array::from_fn(|i| r[i][0] * v[0] + r[i][1] * v[1] + r[i][2] * v[2] + cq[i])
    };
    let distance = |a: [f64; 3], b: [f64; 3]| ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
    for _ in 0..5 {
        let p: Vec<[f64; 3]> = kept.iter().map(|&(_, b)| ca(mobile, b)).collect();
        let q: Vec<[f64; 3]> = kept.iter().map(|&(a, _)| ca(reference, a)).collect();
        fit = kabsch(&p, &q);
        let d: Vec<f64> = p.iter().zip(&q).map(|(p, q)| distance(moved(&fit, *p), *q)).collect();
        let rmsd = (d.iter().map(|d| d * d).sum::<f64>() / d.len() as f64).sqrt();
        // Twice the deviation, or 2 Å for a fit already close.
        let cut = (2.0 * rmsd).max(2.0);
        let next: Vec<(usize, usize)> = kept.iter().zip(&d).filter(|(_, d)| **d <= cut).map(|(k, _)| *k).collect();
        if next.len() == kept.len() || next.len() < 3 {
            break;
        }
        kept = next;
    }
    for atom in &mut mobile.atoms {
        atom.pos = moved(&fit, atom.pos.map(f64::from)).map(|v| v as f32);
    }
    reference.deviation = vec![-1.0; reference.residues.len()];
    mobile.deviation = vec![-1.0; mobile.residues.len()];
    for &(a, b) in &pairs {
        let d = distance(ca(reference, a), ca(mobile, b)) as f32;
        reference.deviation[a] = d;
        mobile.deviation[b] = d;
    }
    let total: f64 = kept.iter().map(|&(a, b)| distance(ca(reference, a), ca(mobile, b)).powi(2)).sum();
    Ok(Overlay {
        rmsd: (total / kept.len() as f64).sqrt() as f32,
        matched: kept.len(),
        aligned: pairs.len(),
        identical,
    })
}

/// How alike two residues are, by one-letter code: a match scores best,
/// residues of a kind (hydrophobic, aromatic, acidic, basic, small polar)
/// a little, anything else counts against.
fn likeness(a: char, b: char) -> i32 {
    if a == b {
        return if a == 'X' { 0 } else { 5 };
    }
    const KINDS: [&str; 6] = ["AVLIMC", "FWYH", "DE", "KR", "STNQ", "GP"];
    if KINDS.iter().any(|k| k.contains(a) && k.contains(b)) {
        1
    } else {
        -2
    }
}

/// A global alignment of two sequences, Gotoh's: Needleman and Wunsch with
/// a gap costing more to open than to grow, and nothing for overhanging
/// ends, so a domain finds its place in a longer chain. Returns the pairs.
pub fn align_sequences(a: &[char], b: &[char]) -> Vec<(usize, usize)> {
    const OPEN: i32 = -10;
    const EXTEND: i32 = -1;
    const NONE: i32 = i32::MIN / 4;
    let (n, m) = (a.len(), b.len());
    if n == 0 || m == 0 {
        return Vec::new();
    }
    // Best scores ending in a match, a gap in b, a gap in a.
    let w = m + 1;
    let mut h = vec![NONE; (n + 1) * w];
    let mut e = vec![NONE; (n + 1) * w];
    let mut f = vec![NONE; (n + 1) * w];
    // Free end gaps: any start along either edge costs nothing.
    for i in 0..=n {
        h[i * w] = 0;
    }
    for j in 0..=m {
        h[j] = 0;
    }
    for i in 1..=n {
        for j in 1..=m {
            let k = i * w + j;
            e[k] = (h[k - w] + OPEN).max(e[k - w] + EXTEND);
            f[k] = (h[k - 1] + OPEN).max(f[k - 1] + EXTEND);
            let diag = h[k - w - 1] + likeness(a[i - 1], b[j - 1]);
            h[k] = diag.max(e[k]).max(f[k]);
        }
    }
    // Start the trace from the best score on the last row or column.
    let mut best = (h[n * w + m], n, m);
    for i in 0..=n {
        if h[i * w + m] > best.0 {
            best = (h[i * w + m], i, m);
        }
    }
    for j in 0..=m {
        if h[n * w + j] > best.0 {
            best = (h[n * w + j], n, j);
        }
    }
    let (_, mut i, mut j) = best;
    let mut pairs = Vec::new();
    // 0 in a match, 1 in a gap in b, 2 in a gap in a.
    let mut state = 0;
    while i > 0 && j > 0 {
        let k = i * w + j;
        match state {
            0 => {
                if h[k] == e[k] {
                    state = 1;
                } else if h[k] == f[k] {
                    state = 2;
                } else {
                    pairs.push((i - 1, j - 1));
                    i -= 1;
                    j -= 1;
                }
            }
            1 => {
                if e[k] != e[k - w] + EXTEND {
                    state = 0;
                }
                i -= 1;
            }
            _ => {
                if f[k] != f[k - 1] + EXTEND {
                    state = 0;
                }
                j -= 1;
            }
        }
    }
    pairs.reverse();
    pairs
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

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(name: &str) -> Structure {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("samples").join(name);
        crate::formats::load_bytes(name, std::fs::read(path).unwrap()).unwrap()
    }

    #[test]
    fn sequences_align_round_an_insertion() {
        let a: Vec<char> = "MKTAYIAKQRQISFVKSHFSRQ".chars().collect();
        let b: Vec<char> = "MKTAYIAKQGGGRQISFVKSHFSRQ".chars().collect();
        let pairs = align_sequences(&a, &b);
        assert_eq!(pairs.len(), a.len());
        assert!(pairs.iter().all(|&(i, j)| a[i] == b[j]));
        // A fragment finds its place in the longer chain.
        let frag: Vec<char> = "ISFVKSH".chars().collect();
        let pairs = align_sequences(&frag, &a);
        assert_eq!(pairs.first(), Some(&(0, 11)));
        assert_eq!(pairs.len(), frag.len());
    }

    #[test]
    fn a_structure_moved_away_comes_back_onto_itself() {
        let reference = sample("4hhb.pdb.gz");
        let mut moved = reference.clone();
        // Turn it a quarter about z and shift it.
        for a in &mut moved.atoms {
            let [x, y, z] = a.pos;
            a.pos = [-y + 20.0, x - 5.0, z + 3.0];
        }
        let mut r = reference.clone();
        let fit = superpose(&mut r, &mut moved).unwrap();
        assert!(fit.rmsd < 0.01, "{}", fit.rmsd);
        assert_eq!(fit.identical, fit.aligned);
        let a = reference.atoms[100].pos;
        let b = moved.atoms[100].pos;
        assert!((0..3).all(|k| (a[k] - b[k]).abs() < 0.01));
    }

    #[test]
    fn haemoglobin_alpha_lies_on_beta() {
        // Alpha and beta chains are 43% identical and fold the same way.
        let s = sample("4hhb.pdb.gz");
        let chain = |id: &str| {
            let mut c = s.clone();
            let keep = c.chains.iter().position(|ch| ch.id == id).unwrap() as u16;
            // The other chains lose their trace, so only this one is paired.
            for r in &mut c.residues {
                if r.chain != keep {
                    r.ca = None;
                }
            }
            c
        };
        let (mut alpha, mut beta) = (chain("A"), chain("B"));
        let fit = superpose(&mut alpha, &mut beta).unwrap();
        assert!(fit.aligned > 120, "{}", fit.aligned);
        assert!(fit.rmsd < 2.0, "{}", fit.rmsd);
        assert!(fit.identical * 100 / fit.aligned > 35);
    }
}
