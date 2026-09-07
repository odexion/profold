//! A small software rasterizer for molecules.
//!
//! GPUI has no 3D pipeline, so we build the frame ourselves: orthographic
//! camera, sphere impostors with a z-buffer, then two screen-space passes —
//! depth cueing and a depth-discontinuity outline — that give the result the
//! flat, printed look of a molecular illustration.

use crate::protein::{Ss, Structure};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Repr {
    /// Ribbons for helices, arrows for strands, cord for the rest.
    Cartoon,
    /// Smooth tube through the alpha carbons.
    Tube,
    /// Every atom at its van der Waals radius.
    Spacefill,
}

impl Repr {
    pub const ALL: [Repr; 3] = [Repr::Cartoon, Repr::Tube, Repr::Spacefill];
    pub fn label(self) -> &'static str {
        match self {
            Repr::Cartoon => "Cartoon",
            Repr::Tube => "Tube",
            Repr::Spacefill => "Spacefill",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Coloring {
    Chain,
    /// Blue at the N terminus through to red at the C terminus.
    Rainbow,
    /// pLDDT, for predicted models.
    Confidence,
    Element,
    /// One colour per structure, when two are superimposed.
    Pair,
    /// Distance between the superimposed counterparts, per residue.
    Deviation,
}

impl Coloring {
    pub const ALL: [Coloring; 4] = [
        Coloring::Chain,
        Coloring::Rainbow,
        Coloring::Confidence,
        Coloring::Element,
    ];
    /// The options that mean something while two structures are overlaid.
    pub const COMPARE: [Coloring; 4] = [
        Coloring::Pair,
        Coloring::Deviation,
        Coloring::Rainbow,
        Coloring::Element,
    ];
    pub fn label(self) -> &'static str {
        match self {
            Coloring::Chain => "Chain",
            Coloring::Rainbow => "Rainbow",
            Coloring::Confidence => "pLDDT",
            Coloring::Element => "Element",
            Coloring::Pair => "Pair",
            Coloring::Deviation => "Deviation",
        }
    }
}

/// The colour ramps a structure can be painted with.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Palette {
    /// The classic molecular rainbow.
    Spectrum,
    /// Perceptually even and safe for colour vision deficiency.
    Viridis,
    Ember,
    Ice,
}

impl Palette {
    pub const ALL: [Palette; 4] = [
        Palette::Spectrum,
        Palette::Viridis,
        Palette::Ember,
        Palette::Ice,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Palette::Spectrum => "Spectrum",
            Palette::Viridis => "Viridis",
            Palette::Ember => "Ember",
            Palette::Ice => "Ice",
        }
    }

    /// Position `t` in 0..=1 along the ramp.
    pub fn ramp(self, t: f32) -> [f32; 3] {
        let t = t.clamp(0.0, 1.0);
        match self {
            // Hue sweep rather than a table, so it stays smooth end to end.
            Palette::Spectrum => hsl_to_rgb((1.0 - t) * 250.0, 0.62, 0.52),
            Palette::Viridis => sample(
                &[
                    [0.267, 0.005, 0.329],
                    [0.229, 0.322, 0.545],
                    [0.128, 0.567, 0.551],
                    [0.369, 0.789, 0.383],
                    [0.993, 0.906, 0.144],
                ],
                t,
            ),
            Palette::Ember => sample(
                &[
                    [0.20, 0.07, 0.32],
                    [0.55, 0.15, 0.40],
                    [0.83, 0.33, 0.32],
                    [0.96, 0.60, 0.25],
                    [0.99, 0.87, 0.44],
                ],
                t,
            ),
            Palette::Ice => sample(
                &[
                    [0.07, 0.12, 0.34],
                    [0.11, 0.37, 0.62],
                    [0.24, 0.64, 0.75],
                    [0.72, 0.92, 0.92],
                ],
                t,
            ),
        }
    }

    /// A distinct colour per chain. Spectrum keeps its hand-picked hues and
    /// Viridis borrows the Okabe–Ito set, which is built to stay separable
    /// under colour vision deficiency; the rest step along their own ramp,
    /// interleaved so neighbouring chains do not land next to each other.
    pub fn chain(self, i: usize) -> [f32; 3] {
        const SPECTRUM: [[f32; 3]; 8] = [
            [0.29, 0.47, 0.78],
            [0.85, 0.44, 0.31],
            [0.36, 0.66, 0.48],
            [0.72, 0.52, 0.80],
            [0.88, 0.70, 0.26],
            [0.36, 0.70, 0.76],
            [0.83, 0.47, 0.62],
            [0.52, 0.55, 0.62],
        ];
        const OKABE_ITO: [[f32; 3]; 8] = [
            [0.000, 0.447, 0.698],
            [0.902, 0.624, 0.000],
            [0.000, 0.620, 0.451],
            [0.800, 0.475, 0.655],
            [0.337, 0.706, 0.914],
            [0.835, 0.369, 0.000],
            [0.941, 0.894, 0.259],
            [0.350, 0.350, 0.350],
        ];
        const ORDER: [usize; 8] = [0, 4, 2, 6, 1, 5, 3, 7];
        match self {
            Palette::Spectrum => SPECTRUM[i % 8],
            Palette::Viridis => OKABE_ITO[i % 8],
            _ => self.ramp(ORDER[i % 8] as f32 / 7.0),
        }
    }
}

/// Even stops along a ramp table.
fn sample(stops: &[[f32; 3]], t: f32) -> [f32; 3] {
    let last = stops.len() - 1;
    let x = t.clamp(0.0, 1.0) * last as f32;
    let i = (x.floor() as usize).min(last - 1);
    mix(stops[i], stops[i + 1], x - i as f32)
}

/// How the frame is shaded.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Style {
    /// Smooth shading on paper, with distance washing out.
    Soft,
    /// Ink outlines, vivid colour, banded light: the look of a drawn figure.
    Illustrative,
}

impl Style {
    pub const ALL: [Style; 2] = [Style::Soft, Style::Illustrative];
    pub fn label(self) -> &'static str {
        match self {
            Style::Soft => "Soft",
            Style::Illustrative => "Illustrative",
        }
    }
}

#[derive(Clone, Copy)]
pub struct Camera {
    pub yaw: f32,
    pub pitch: f32,
    /// Multiplier on the auto-fit scale.
    pub zoom: f32,
    pub pan: [f32; 2],
    pub center: [f32; 3],
    pub radius: f32,
}

impl Camera {
    pub fn frame(structure: &Structure) -> Self {
        Self::frame_all(&[structure])
    }

    /// Frame everything on screen, so a superimposed pair stays in view even
    /// where the prediction wanders away from the reference.
    pub fn frame_all(structures: &[&Structure]) -> Self {
        let mut lo = [f32::MAX; 3];
        let mut hi = [f32::MIN; 3];
        for s in structures {
            for a in &s.atoms {
                for i in 0..3 {
                    lo[i] = lo[i].min(a.pos[i]);
                    hi[i] = hi[i].max(a.pos[i]);
                }
            }
        }
        if lo[0] > hi[0] {
            return Camera {
                yaw: 0.55,
                pitch: 0.25,
                zoom: 1.0,
                pan: [0.0, 0.0],
                center: [0.0; 3],
                radius: 1.0,
            };
        }
        let center = [
            (lo[0] + hi[0]) * 0.5,
            (lo[1] + hi[1]) * 0.5,
            (lo[2] + hi[2]) * 0.5,
        ];
        let mut r2: f32 = 0.0;
        for s in structures {
            for a in &s.atoms {
                let d = [
                    a.pos[0] - center[0],
                    a.pos[1] - center[1],
                    a.pos[2] - center[2],
                ];
                r2 = r2.max(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            }
        }
        let radius = r2.sqrt().max(1.0);
        Camera {
            yaw: 0.55,
            pitch: 0.25,
            zoom: 1.0,
            pan: [0.0, 0.0],
            center,
            radius,
        }
    }

    fn basis(&self) -> [[f32; 3]; 3] {
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();
        // Rows of Rx(pitch) * Ry(yaw): world -> view.
        [
            [cy, 0.0, -sy],
            [sp * sy, cp, sp * cy],
            [cp * sy, -sp, cp * cy],
        ]
    }
}

#[derive(Clone, Copy)]
pub struct Sphere {
    pub c: [f32; 3],
    pub r: f32,
    pub color: [f32; 3],
}

/// A sphere-swept segment. The tube is built from these rather than from a
/// chain of spheres: consecutive capsules share an end cap exactly, so the
/// surface has no seams to shade.
#[derive(Clone, Copy)]
pub struct Capsule {
    pub a: [f32; 3],
    pub b: [f32; 3],
    pub r: f32,
    pub ca: [f32; 3],
    pub cb: [f32; 3],
}

/// A shaded triangle with per-vertex normals, used for cartoon ribbons.
#[derive(Clone, Copy)]
pub struct Tri {
    pub p: [[f32; 3]; 3],
    pub n: [[f32; 3]; 3],
    pub c: [f32; 3],
}

#[derive(Default)]
pub struct Scene {
    pub spheres: Vec<Sphere>,
    pub capsules: Vec<Capsule>,
    pub tris: Vec<Tri>,
}

impl Scene {
    pub fn len(&self) -> usize {
        self.spheres.len() + self.capsules.len() + self.tris.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn merge(&mut self, mut other: Scene) {
        self.spheres.append(&mut other.spheres);
        self.capsules.append(&mut other.capsules);
        self.tris.append(&mut other.tris);
    }
}

// ---------------------------------------------------------------------------
// Scene construction
// ---------------------------------------------------------------------------

/// The AlphaFold pLDDT bands, which people read at a glance.
fn plddt_color(v: f32) -> [f32; 3] {
    if v >= 90.0 {
        [0.05, 0.34, 0.77]
    } else if v >= 70.0 {
        [0.40, 0.79, 0.94]
    } else if v >= 50.0 {
        [1.00, 0.86, 0.36]
    } else {
        [0.98, 0.49, 0.24]
    }
}

fn rgb_to_hsl(c: [f32; 3]) -> (f32, f32, f32) {
    let (r, g, b) = (c[0], c[1], c[2]);
    let max = r.max(g).max(b);
    let min = r.min(g).min(b);
    let l = (max + min) * 0.5;
    if (max - min).abs() < 1e-6 {
        return (0.0, 0.0, l);
    }
    let d = max - min;
    let s = if l > 0.5 {
        d / (2.0 - max - min)
    } else {
        d / (max + min)
    };
    let h = if max == r {
        ((g - b) / d + if g < b { 6.0 } else { 0.0 }) * 60.0
    } else if max == g {
        ((b - r) / d + 2.0) * 60.0
    } else {
        ((r - g) / d + 4.0) * 60.0
    };
    (h, s, l)
}

/// Push a colour toward poster ink: saturated, and away from the extremes of
/// lightness so it still holds an outline.
fn vivid(c: [f32; 3]) -> [f32; 3] {
    let (h, s, l) = rgb_to_hsl(c);
    if s < 0.02 {
        return [l, l, l];
    }
    hsl_to_rgb(h, (s * 1.5 + 0.18).min(1.0), 0.5 + (l - 0.5) * 0.62)
}

fn hsl_to_rgb(h: f32, s: f32, l: f32) -> [f32; 3] {
    let c = (1.0 - (2.0 * l - 1.0).abs()) * s;
    let hp = h / 60.0;
    let x = c * (1.0 - (hp % 2.0 - 1.0).abs());
    let (r, g, b) = match hp as i32 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };
    let m = l - c * 0.5;
    [r + m, g + m, b + m]
}

/// Deviation in angstrom, along the chosen ramp. Compressed at the low end,
/// where a tenth of an angstrom still means something, and saturating by 8 A,
/// past which the two have simply gone their separate ways.
fn deviation_color(d: f32, palette: Palette) -> [f32; 3] {
    let t = (d.max(0.0) / 8.0).powf(0.6).clamp(0.0, 1.0);
    palette.ramp(t)
}

fn residue_color(
    s: &Structure,
    residue: usize,
    coloring: Coloring,
    palette: Palette,
) -> [f32; 3] {
    let r = &s.residues[residue];
    match coloring {
        Coloring::Pair => s.tint.unwrap_or(palette.chain(r.chain as usize)),
        Coloring::Deviation => match s.deviation.get(residue) {
            Some(d) if *d >= 0.0 => deviation_color(*d, palette),
            // No counterpart to compare against.
            _ => [0.74, 0.73, 0.71],
        },
        Coloring::Chain => palette.chain(r.chain as usize),
        Coloring::Rainbow => palette.ramp(residue as f32 / (s.residues.len().max(2) - 1) as f32),
        // pLDDT keeps AlphaFold's own bands: they are a published convention
        // and a reader knows them by sight.
        Coloring::Confidence => {
            if s.has_confidence {
                plddt_color(r.b)
            } else {
                palette.chain(r.chain as usize)
            }
        }
        Coloring::Element => [0.42, 0.44, 0.47],
    }
}

fn catmull_rom(p0: [f32; 3], p1: [f32; 3], p2: [f32; 3], p3: [f32; 3], t: f32) -> [f32; 3] {
    let t2 = t * t;
    let t3 = t2 * t;
    let mut out = [0.0; 3];
    for i in 0..3 {
        out[i] = 0.5
            * ((2.0 * p1[i])
                + (-p0[i] + p2[i]) * t
                + (2.0 * p0[i] - 5.0 * p1[i] + 4.0 * p2[i] - p3[i]) * t2
                + (-p0[i] + 3.0 * p1[i] - 3.0 * p2[i] + p3[i]) * t3);
    }
    out
}

/// Build the drawable primitives for a structure. `quality` in 0..=1 trades
/// tube smoothness for speed on very long chains.
pub fn build_scene(
    s: &Structure,
    repr: Repr,
    coloring: Coloring,
    quality: f32,
    style: Style,
    palette: Palette,
) -> Scene {
    let mut scene = build_primitives(s, repr, coloring, quality, palette);
    if style == Style::Illustrative {
        for c in scene
            .spheres
            .iter_mut()
            .map(|p| &mut p.color)
            .chain(scene.capsules.iter_mut().flat_map(|c| [&mut c.ca, &mut c.cb]))
            .chain(scene.tris.iter_mut().map(|t| &mut t.c))
        {
            *c = vivid(*c);
        }
    }
    scene
}

fn build_primitives(
    s: &Structure,
    repr: Repr,
    coloring: Coloring,
    quality: f32,
    palette: Palette,
) -> Scene {
    match repr {
        Repr::Spacefill => Scene {
            spheres: s
                .atoms
                .iter()
                .map(|a| Sphere {
                    c: a.pos,
                    r: a.element.vdw(),
                    color: if coloring == Coloring::Element {
                        a.element.color()
                    } else {
                        residue_color(s, a.residue as usize, coloring, palette)
                    },
                })
                .collect(),
            capsules: Vec::new(),
            tris: Vec::new(),
        },
        Repr::Tube => Scene {
            spheres: Vec::new(),
            capsules: build_tube(s, coloring, quality, palette),
            tris: Vec::new(),
        },
        Repr::Cartoon => Scene {
            spheres: Vec::new(),
            capsules: Vec::new(),
            tris: build_cartoon(s, coloring, quality, palette),
        },
    }
}

// -- small vector helpers ----------------------------------------------------

#[inline]
fn sub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
fn unit(v: [f32; 3]) -> [f32; 3] {
    let n = dot(v, v).sqrt();
    if n < 1e-6 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}
#[inline]
fn mix(a: [f32; 3], b: [f32; 3], t: f32) -> [f32; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

/// Alpha-carbon runs for one chain, split wherever the trace jumps.
fn ca_runs(s: &Structure, chain: &crate::protein::Chain) -> Vec<Vec<usize>> {
    let mut runs: Vec<Vec<usize>> = Vec::new();
    let mut run: Vec<usize> = Vec::new();
    let cut = |run: &mut Vec<usize>, runs: &mut Vec<Vec<usize>>| {
        if run.len() > 1 {
            runs.push(std::mem::take(run));
        } else {
            run.clear();
        }
    };
    for &ri in &chain.residues {
        let ri = ri as usize;
        let Some(ca) = s.residues[ri].ca else {
            cut(&mut run, &mut runs);
            continue;
        };
        let p = s.atoms[ca as usize].pos;
        if let Some(&prev) = run.last() {
            let q = s.atoms[s.residues[prev].ca.unwrap() as usize].pos;
            // Consecutive alpha carbons sit 3.8 A apart; 5.5 A means a gap.
            if dot(sub(p, q), sub(p, q)) > 5.5 * 5.5 {
                cut(&mut run, &mut runs);
            }
        }
        run.push(ri);
    }
    cut(&mut run, &mut runs);
    runs
}

const TUBE_RADIUS: f32 = 0.75;

fn build_tube(
    s: &Structure,
    coloring: Coloring,
    quality: f32,
    palette: Palette,
) -> Vec<Capsule> {
    let mut out = Vec::new();

    for chain in &s.chains {
        let runs = ca_runs(s, chain);
        for run in runs {
            let run: Vec<(usize, [f32; 3])> = run
                .iter()
                .map(|&ri| (ri, s.atoms[s.residues[ri].ca.unwrap() as usize].pos))
                .collect();
            let steps = (14.0 * quality).round().max(2.0) as usize;
            let at = |i: isize| -> [f32; 3] {
                run[i.clamp(0, run.len() as isize - 1) as usize].1
            };

            // Sample the Catmull-Rom spline through the trace, then join the
            // samples with capsules.
            let mut points: Vec<([f32; 3], [f32; 3])> = Vec::with_capacity(run.len() * steps + 1);
            for i in 0..run.len() - 1 {
                let (p0, p1, p2, p3) = (
                    at(i as isize - 1),
                    at(i as isize),
                    at(i as isize + 1),
                    at(i as isize + 2),
                );
                let c0 = residue_color(s, run[i].0, coloring, palette);
                let c1 = residue_color(s, run[i + 1].0, coloring, palette);
                for k in 0..steps {
                    let t = k as f32 / steps as f32;
                    points.push((
                        catmull_rom(p0, p1, p2, p3, t),
                        [
                            c0[0] + (c1[0] - c0[0]) * t,
                            c0[1] + (c1[1] - c0[1]) * t,
                            c0[2] + (c1[2] - c0[2]) * t,
                        ],
                    ));
                }
            }
            let last = run[run.len() - 1];
            points.push((last.1, residue_color(s, last.0, coloring, palette)));

            for pair in points.windows(2) {
                out.push(Capsule {
                    a: pair[0].0,
                    b: pair[1].0,
                    r: TUBE_RADIUS,
                    ca: pair[0].1,
                    cb: pair[1].1,
                });
            }
        }
    }

    out
}


// ---------------------------------------------------------------------------
// Cartoon
// ---------------------------------------------------------------------------

/// Half-width and half-thickness of the swept cross-section, per structure.
const COIL_R: f32 = 0.30;
const HELIX: (f32, f32) = (1.15, 0.22);
const SHEET: (f32, f32) = (1.20, 0.22);
/// The arrowhead flares to this before tapering to a point.
const ARROW_HW: f32 = 1.95;
/// Sides of the swept cross-section. Enough that a coil reads as round.
const SIDES: usize = 10;

/// A ribbon that follows the raw alpha carbons of a strand shows the pleat, so
/// sheet residues get a light smoothing pass first. Helices keep their spiral.
fn smooth_sheets(s: &Structure, run: &[usize], p: &mut [[f32; 3]]) {
    for _ in 0..2 {
        let src = p.to_vec();
        for k in 1..run.len() - 1 {
            if s.residues[run[k]].ss != Ss::Sheet {
                continue;
            }
            p[k] = [
                (src[k - 1][0] + 2.0 * src[k][0] + src[k + 1][0]) * 0.25,
                (src[k - 1][1] + 2.0 * src[k][1] + src[k + 1][1]) * 0.25,
                (src[k - 1][2] + 2.0 * src[k][2] + src[k + 1][2]) * 0.25,
            ];
        }
    }
}

fn build_cartoon(
    s: &Structure,
    coloring: Coloring,
    quality: f32,
    palette: Palette,
) -> Vec<Tri> {
    let mut tris = Vec::new();
    let steps = (8.0 * quality).round().max(2.0) as usize;

    for chain in &s.chains {
        for run in ca_runs(s, chain) {
            if run.len() < 3 {
                continue;
            }
            let n = run.len();
            let mut p: Vec<[f32; 3]> = run
                .iter()
                .map(|&ri| s.atoms[s.residues[ri].ca.unwrap() as usize].pos)
                .collect();
            smooth_sheets(s, &run, &mut p);

            // Per-residue frame. The ribbon's flat face is set by the carbonyl
            // direction, which is what makes a helix read as a coiled ribbon
            // rather than a twisted one.
            let mut side: Vec<[f32; 3]> = Vec::with_capacity(n);
            for k in 0..n {
                let t = if k == 0 {
                    unit(sub(p[1], p[0]))
                } else if k == n - 1 {
                    unit(sub(p[n - 1], p[n - 2]))
                } else {
                    unit(sub(p[k + 1], p[k - 1]))
                };
                let r = &s.residues[run[k]];
                let raw = match (r.c, r.o) {
                    (Some(c), Some(o)) => sub(s.atoms[o as usize].pos, s.atoms[c as usize].pos),
                    // No carbonyl: fall back to the local bend.
                    _ if k > 0 && k + 1 < n => cross(sub(p[k + 1], p[k]), sub(p[k - 1], p[k])),
                    _ => [0.0, 0.0, 1.0],
                };
                let mut w = unit(sub(raw, [t[0] * dot(raw, t), t[1] * dot(raw, t), t[2] * dot(raw, t)]));
                if dot(w, w) < 0.5 {
                    // Degenerate: any vector perpendicular to the tangent.
                    let a = if t[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
                    w = unit(cross(t, a));
                }
                // A strand alternates its carbonyls; keep the face continuous.
                if let Some(&prev) = side.last() {
                    if dot(w, prev) < 0.0 {
                        w = [-w[0], -w[1], -w[2]];
                    }
                }
                side.push(w);
            }

            let ss = |k: usize| s.residues[run[k]].ss;
            let profile = |k: usize| match ss(k) {
                Ss::Helix => HELIX,
                Ss::Sheet => SHEET,
                Ss::Coil => (COIL_R, COIL_R),
            };
            let arrow = |k: usize| ss(k) == Ss::Sheet && (k + 1 >= n || ss(k + 1) != Ss::Sheet);

            // Sample the spline, carrying the frame and cross-section with it.
            struct Sample {
                pos: [f32; 3],
                side: [f32; 3],
                hw: f32,
                ht: f32,
                color: [f32; 3],
            }
            let at = |i: isize| p[i.clamp(0, n as isize - 1) as usize];
            let mut samples: Vec<Sample> = Vec::with_capacity(n * steps + 1);
            for k in 0..n - 1 {
                let c0 = residue_color(s, run[k], coloring, palette);
                let c1 = residue_color(s, run[k + 1], coloring, palette);
                let (hw0, ht0) = profile(k);
                let (hw1, ht1) = profile(k + 1);
                for step in 0..steps {
                    let f = step as f32 / steps as f32;
                    let (hw, ht) = if arrow(k) {
                        // Flare at the head, then taper to a point.
                        (ARROW_HW * (1.0 - f) + 0.12 * f, ht0)
                    } else if arrow(k + 1) {
                        // Hold the strand width; the flare is a deliberate step.
                        (hw0, ht0)
                    } else {
                        (hw0 + (hw1 - hw0) * f, ht0 + (ht1 - ht0) * f)
                    };
                    samples.push(Sample {
                        pos: catmull_rom(
                            at(k as isize - 1),
                            at(k as isize),
                            at(k as isize + 1),
                            at(k as isize + 2),
                            f,
                        ),
                        side: mix(side[k], side[k + 1], f),
                        hw,
                        ht,
                        color: mix(c0, c1, f),
                    });
                }
            }
            let (hw, ht) = profile(n - 1);
            samples.push(Sample {
                pos: p[n - 1],
                side: side[n - 1],
                hw: if arrow(n - 1) { 0.12 } else { hw },
                ht,
                color: residue_color(s, run[n - 1], coloring, palette),
            });
            if samples.len() < 2 {
                continue;
            }

            // Turn each sample into a ring of points on an elliptical section.
            let count = samples.len();
            let mut rings: Vec<[([f32; 3], [f32; 3]); SIDES]> = Vec::with_capacity(count);
            for (i, sm) in samples.iter().enumerate() {
                let tan = if i == 0 {
                    unit(sub(samples[1].pos, sm.pos))
                } else if i == count - 1 {
                    unit(sub(sm.pos, samples[count - 2].pos))
                } else {
                    unit(sub(samples[i + 1].pos, samples[i - 1].pos))
                };
                let w = unit(sub(
                    sm.side,
                    [
                        tan[0] * dot(sm.side, tan),
                        tan[1] * dot(sm.side, tan),
                        tan[2] * dot(sm.side, tan),
                    ],
                ));
                let up = cross(tan, w);
                rings.push(std::array::from_fn(|j| {
                    let a = std::f32::consts::TAU * j as f32 / SIDES as f32;
                    let (sa, ca) = a.sin_cos();
                    let pos = [
                        sm.pos[0] + w[0] * sm.hw * ca + up[0] * sm.ht * sa,
                        sm.pos[1] + w[1] * sm.hw * ca + up[1] * sm.ht * sa,
                        sm.pos[2] + w[2] * sm.hw * ca + up[2] * sm.ht * sa,
                    ];
                    // Outward normal of an ellipse swaps the two half-axes.
                    let nrm = unit([
                        w[0] * sm.ht * ca + up[0] * sm.hw * sa,
                        w[1] * sm.ht * ca + up[1] * sm.hw * sa,
                        w[2] * sm.ht * ca + up[2] * sm.hw * sa,
                    ]);
                    (pos, nrm)
                }));
            }

            for i in 0..count - 1 {
                let (a, b) = (&rings[i], &rings[i + 1]);
                let color = samples[i].color;
                for j in 0..SIDES {
                    let k = (j + 1) % SIDES;
                    tris.push(Tri {
                        p: [a[j].0, a[k].0, b[k].0],
                        n: [a[j].1, a[k].1, b[k].1],
                        c: color,
                    });
                    tris.push(Tri {
                        p: [a[j].0, b[k].0, b[j].0],
                        n: [a[j].1, b[k].1, b[j].1],
                        c: color,
                    });
                }
            }

            // Close both ends so the ribbon is not hollow.
            for (ring, color) in [
                (&rings[0], samples[0].color),
                (&rings[count - 1], samples[count - 1].color),
            ] {
                let nrm = unit(cross(
                    sub(ring[1].0, ring[0].0),
                    sub(ring[2].0, ring[0].0),
                ));
                for j in 1..SIDES - 1 {
                    tris.push(Tri {
                        p: [ring[0].0, ring[j].0, ring[j + 1].0],
                        n: [nrm; 3],
                        c: color,
                    });
                }
            }
        }
    }

    tris
}

// ---------------------------------------------------------------------------
// Rasterization
// ---------------------------------------------------------------------------

pub struct Frame {
    pub width: usize,
    pub height: usize,
    /// BGRA, which is what GPUI uploads directly.
    pub bgra: Vec<u8>,
}

/// Vertical paper gradient behind the molecule.
#[inline]
fn background(y: usize, height: usize, top: [f32; 3], bottom: [f32; 3]) -> [f32; 3] {
    let t = y as f32 / height.max(1) as f32;
    [
        top[0] + (bottom[0] - top[0]) * t,
        top[1] + (bottom[1] - top[1]) * t,
        top[2] + (bottom[2] - top[2]) * t,
    ]
}

/// Light from the upper left, slightly in front of the viewer.
const LIGHT: [f32; 3] = [-0.4082, 0.5585, 0.7226];
/// Halfway between the light and the view direction, for the highlight.
const HALF: [f32; 3] = [-0.2225, 0.3044, 0.9262];

#[inline]
fn shade(n: [f32; 3], base: [f32; 3], style: Style) -> [f32; 3] {
    let diff = (n[0] * LIGHT[0] + n[1] * LIGHT[1] + n[2] * LIGHT[2]).max(0.0);
    let nh = (n[0] * HALF[0] + n[1] * HALF[1] + n[2] * HALF[2]).max(0.0);
    let n2 = nh * nh;
    let n8 = n2 * n2 * n2 * n2;

    let (k, spec) = match style {
        Style::Soft => (0.34 + 0.70 * diff, n8 * n8 * n8 * 0.28),
        Style::Illustrative => {
            // Quantise the light into a few steps, then ease most of the way
            // back: enough banding to read as drawn, not enough to look harsh.
            const STEPS: f32 = 4.0;
            let banded = (diff * STEPS).floor() / (STEPS - 1.0);
            let mixed = diff * 0.35 + banded.min(1.0) * 0.65;
            (0.62 + 0.42 * mixed, n8 * n8 * 0.10)
        }
    };
    [base[0] * k + spec, base[1] * k + spec, base[2] * k + spec]
}

struct Target {
    color: Vec<f32>,
    depth: Vec<f32>,
    style: Style,
}

impl Target {
    #[inline]
    fn plot(&mut self, idx: usize, z: f32, n: [f32; 3], base: [f32; 3]) {
        if z <= self.depth[idx] {
            return;
        }
        self.depth[idx] = z;
        let c = shade(n, base, self.style);
        let o = idx * 3;
        self.color[o] = c[0];
        self.color[o + 1] = c[1];
        self.color[o + 2] = c[2];
    }
}

/// Ray-sphere in orthographic view space: returns the near-surface depth and
/// normal for the ray through `(x, y)`.
#[inline]
fn sphere_hit(c: [f32; 3], r: f32, x: f32, y: f32) -> Option<(f32, [f32; 3])> {
    let dx = x - c[0];
    let dy = y - c[1];
    let d2 = dx * dx + dy * dy;
    let rr = r * r;
    if d2 > rr {
        return None;
    }
    let dz = (rr - d2).sqrt();
    Some((c[2] + dz, [dx / r, dy / r, dz / r]))
}

pub fn render(
    scene: &Scene,
    cam: &Camera,
    width: usize,
    height: usize,
    bg_top: [f32; 3],
    bg_bottom: [f32; 3],
    style: Style,
) -> Frame {
    let n = width * height;
    let mut t = Target {
        color: vec![0.0; n * 3],
        depth: vec![f32::NEG_INFINITY; n],
        style,
    };

    for y in 0..height {
        let bg = background(y, height, bg_top, bg_bottom);
        for x in 0..width {
            let o = (y * width + x) * 3;
            t.color[o] = bg[0];
            t.color[o + 1] = bg[1];
            t.color[o + 2] = bg[2];
        }
    }

    if scene.is_empty() || width == 0 || height == 0 {
        return to_bgra(t.color, width, height);
    }

    let m = cam.basis();
    let scale = (width.min(height) as f32) / (2.2 * cam.radius) * cam.zoom;
    let cx = width as f32 * 0.5 + cam.pan[0];
    let cy = height as f32 * 0.5 + cam.pan[1];
    let inv = 1.0 / scale;

    let to_view = |p: [f32; 3]| {
        let q = [p[0] - cam.center[0], p[1] - cam.center[1], p[2] - cam.center[2]];
        [
            m[0][0] * q[0] + m[0][1] * q[1] + m[0][2] * q[2],
            m[1][0] * q[0] + m[1][1] * q[1] + m[1][2] * q[2],
            m[2][0] * q[0] + m[2][1] * q[1] + m[2][2] * q[2],
        ]
    };
    let rot = |v: [f32; 3]| {
        [
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
        ]
    };
    // View-space coordinates of a pixel centre.
    let view_x = |px: usize| (px as f32 + 0.5 - cx) * inv;
    let view_y = |py: usize| -(py as f32 + 0.5 - cy) * inv;

    let mut zmin = f32::MAX;
    let mut zmax = f32::MIN;
    let mut track = |lo: f32, hi: f32| {
        if lo < zmin {
            zmin = lo;
        }
        if hi > zmax {
            zmax = hi;
        }
    };

    // -- spheres ----------------------------------------------------------
    for s in &scene.spheres {
        let v = to_view(s.c);
        track(v[2] - s.r, v[2] + s.r);

        let sx = cx + v[0] * scale;
        let sy = cy - v[1] * scale;
        let rp = s.r * scale;
        let Some((x0, x1, y0, y1)) = bbox(sx - rp, sx + rp, sy - rp, sy + rp, width, height) else {
            continue;
        };
        for py in y0..y1 {
            let y = view_y(py);
            let row = py * width;
            for px in x0..x1 {
                if let Some((z, n)) = sphere_hit(v, s.r, view_x(px), y) {
                    t.plot(row + px, z, n, s.color);
                }
            }
        }
    }

    // -- capsules ---------------------------------------------------------
    for c in &scene.capsules {
        let a = to_view(c.a);
        let b = to_view(c.b);
        track(a[2].min(b[2]) - c.r, a[2].max(b[2]) + c.r);

        let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let u = if len > 1e-6 {
            [d[0] / len, d[1] / len, d[2] / len]
        } else {
            [0.0, 0.0, 1.0]
        };

        let sax = cx + a[0] * scale;
        let say = cy - a[1] * scale;
        let sbx = cx + b[0] * scale;
        let sby = cy - b[1] * scale;
        let rp = c.r * scale;
        let Some((x0, x1, y0, y1)) = bbox(
            sax.min(sbx) - rp,
            sax.max(sbx) + rp,
            say.min(sby) - rp,
            say.max(sby) + rp,
            width,
            height,
        ) else {
            continue;
        };

        // Coefficients of the ray-cylinder quadratic that do not vary per
        // pixel. `q` collapses to zero when the axis points at the viewer, in
        // which case the end caps alone describe the shape.
        let q = 1.0 - u[2] * u[2];
        let rr = c.r * c.r;

        for py in y0..y1 {
            let y = view_y(py);
            let row = py * width;
            for px in x0..x1 {
                let x = view_x(px);
                let wx = x - a[0];
                let wy = y - a[1];

                if q > 1e-6 && len > 1e-6 {
                    let k = wx * u[0] + wy * u[1];
                    let bq = -2.0 * k * u[2];
                    let cq = wx * wx + wy * wy - k * k - rr;
                    let disc = bq * bq - 4.0 * q * cq;
                    if disc >= 0.0 {
                        let cz = (-bq + disc.sqrt()) / (2.0 * q);
                        let along = k + cz * u[2];
                        let s = along / len;
                        if (0.0..=1.0).contains(&s) {
                            let nrm = [
                                (wx - u[0] * along) / c.r,
                                (wy - u[1] * along) / c.r,
                                (cz - u[2] * along) / c.r,
                            ];
                            let base = [
                                c.ca[0] + (c.cb[0] - c.ca[0]) * s,
                                c.ca[1] + (c.cb[1] - c.ca[1]) * s,
                                c.ca[2] + (c.cb[2] - c.ca[2]) * s,
                            ];
                            t.plot(row + px, a[2] + cz, nrm, base);
                            continue;
                        }
                    }
                }

                // Outside the cylinder's span, the caps take over.
                let mut best: Option<(f32, [f32; 3], [f32; 3])> = None;
                for (centre, color) in [(a, c.ca), (b, c.cb)] {
                    if let Some((z, nrm)) = sphere_hit(centre, c.r, x, y) {
                        if best.as_ref().is_none_or(|(bz, _, _)| z > *bz) {
                            best = Some((z, nrm, color));
                        }
                    }
                }
                if let Some((z, nrm, color)) = best {
                    t.plot(row + px, z, nrm, color);
                }
            }
        }
    }

    // -- ribbon meshes ----------------------------------------------------
    for tri in &scene.tris {
        let v = [to_view(tri.p[0]), to_view(tri.p[1]), to_view(tri.p[2])];
        let nv = [rot(tri.n[0]), rot(tri.n[1]), rot(tri.n[2])];
        track(
            v[0][2].min(v[1][2]).min(v[2][2]),
            v[0][2].max(v[1][2]).max(v[2][2]),
        );

        let sx: [f32; 3] = std::array::from_fn(|i| cx + v[i][0] * scale);
        let sy: [f32; 3] = std::array::from_fn(|i| cy - v[i][1] * scale);

        let d = (sy[1] - sy[2]) * (sx[0] - sx[2]) + (sx[2] - sx[1]) * (sy[0] - sy[2]);
        if d.abs() < 1e-9 {
            continue;
        }
        let inv = 1.0 / d;

        let Some((x0, x1, y0, y1)) = bbox(
            sx[0].min(sx[1]).min(sx[2]),
            sx[0].max(sx[1]).max(sx[2]),
            sy[0].min(sy[1]).min(sy[2]),
            sy[0].max(sy[1]).max(sy[2]),
            width,
            height,
        ) else {
            continue;
        };

        for py in y0..y1 {
            let fy = py as f32 + 0.5;
            let row = py * width;
            for px in x0..x1 {
                let fx = px as f32 + 0.5;
                // Barycentric coordinates; dividing by the signed area keeps
                // them correct whichever way the triangle winds.
                let l0 = ((sy[1] - sy[2]) * (fx - sx[2]) + (sx[2] - sx[1]) * (fy - sy[2])) * inv;
                let l1 = ((sy[2] - sy[0]) * (fx - sx[2]) + (sx[0] - sx[2]) * (fy - sy[2])) * inv;
                let l2 = 1.0 - l0 - l1;
                if l0 < 0.0 || l1 < 0.0 || l2 < 0.0 {
                    continue;
                }

                // Orthographic projection is affine, so depth interpolates
                // linearly and exactly.
                let z = l0 * v[0][2] + l1 * v[1][2] + l2 * v[2][2];
                let idx = row + px;
                if z <= t.depth[idx] {
                    continue;
                }

                let mut n = [
                    l0 * nv[0][0] + l1 * nv[1][0] + l2 * nv[2][0],
                    l0 * nv[0][1] + l1 * nv[1][1] + l2 * nv[2][1],
                    l0 * nv[0][2] + l1 * nv[1][2] + l2 * nv[2][2],
                ];
                let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                if len < 1e-6 {
                    continue;
                }
                n = [n[0] / len, n[1] / len, n[2] / len];
                // A ribbon is seen from both sides; light whichever faces us.
                if n[2] < 0.0 {
                    n = [-n[0], -n[1], -n[2]];
                }
                t.plot(idx, z, n, tri.c);
            }
        }
    }

    let span = (zmax - zmin).max(1e-3);

    // Depth cue: wash distant geometry toward the paper. A drawn figure does
    // not fade with distance, so the outlines carry the depth instead.
    let fog_strength = match style {
        Style::Soft => 0.42,
        Style::Illustrative => 0.0,
    };
    for y in 0..height {
        let bg = background(y, height, bg_top, bg_bottom);
        for x in 0..width {
            let idx = y * width + x;
            if t.depth[idx] == f32::NEG_INFINITY {
                continue;
            }
            let d = ((t.depth[idx] - zmin) / span).clamp(0.0, 1.0);
            let fog = (1.0 - d) * fog_strength;
            let o = idx * 3;
            for k in 0..3 {
                t.color[o + k] += (bg[k] - t.color[o + k]) * fog;
            }
        }
    }

    // Outline. Soft darkens what sits behind a depth jump and the rim of the
    // silhouette; illustrative lays down real ink, on both sides of the
    // silhouette, thickened so it holds at any resolution.
    let jump = (span * 0.012).max(0.45);
    match style {
        Style::Soft => {
            let src = t.color.clone();
            for y in 0..height {
                for x in 0..width {
                    let idx = y * width + x;
                    let z = t.depth[idx];
                    if z == f32::NEG_INFINITY {
                        continue;
                    }
                    let mut behind = false;
                    let mut open = false;
                    for nb in neighbours(idx, x, y, width, height).into_iter().flatten() {
                        let nz = t.depth[nb];
                        if nz == f32::NEG_INFINITY {
                            open = true;
                        } else if nz - z > jump {
                            behind = true;
                        }
                    }
                    let k = if behind {
                        0.62
                    } else if open {
                        0.80
                    } else {
                        continue;
                    };
                    let o = idx * 3;
                    for i in 0..3 {
                        t.color[o + i] = src[o + i] * k;
                    }
                }
            }
        }
        Style::Illustrative => {
            const INK: [f32; 3] = [0.07, 0.07, 0.09];
            let mut edge = vec![false; n];
            for y in 0..height {
                for x in 0..width {
                    let idx = y * width + x;
                    let z = t.depth[idx];
                    for nb in neighbours(idx, x, y, width, height).into_iter().flatten() {
                        let nz = t.depth[nb];
                        // The silhouette is inked from both sides, so the line
                        // reads the same against the page and against the body.
                        let silhouette = (z == f32::NEG_INFINITY) != (nz == f32::NEG_INFINITY);
                        let overlap =
                            z != f32::NEG_INFINITY && nz != f32::NEG_INFINITY && nz - z > jump;
                        if silhouette || overlap {
                            edge[idx] = true;
                            break;
                        }
                    }
                }
            }

            // One pixel of line per ~700 of frame, so the weight is the same
            // whether or not the display is retina.
            let weight = (width.min(height) / 700).max(1);
            for _ in 1..weight {
                let src = edge.clone();
                for y in 0..height {
                    for x in 0..width {
                        let idx = y * width + x;
                        if src[idx] {
                            continue;
                        }
                        edge[idx] = neighbours(idx, x, y, width, height)
                            .into_iter()
                            .flatten()
                            .any(|nb| src[nb]);
                    }
                }
            }

            for (idx, _) in edge.iter().enumerate().filter(|(_, e)| **e) {
                let o = idx * 3;
                t.color[o] = INK[0];
                t.color[o + 1] = INK[1];
                t.color[o + 2] = INK[2];
            }
        }
    }

    to_bgra(t.color, width, height)
}

/// The four pixels sharing an edge with this one.
#[inline]
fn neighbours(
    idx: usize,
    x: usize,
    y: usize,
    width: usize,
    height: usize,
) -> [Option<usize>; 4] {
    [
        (x > 0).then(|| idx - 1),
        (x + 1 < width).then_some(idx + 1),
        (y > 0).then(|| idx - width),
        (y + 1 < height).then_some(idx + width),
    ]
}

/// Clamp a floating-point rectangle to the frame, or `None` if it falls
/// entirely outside.
#[inline]
fn bbox(
    left: f32,
    right: f32,
    top: f32,
    bottom: f32,
    width: usize,
    height: usize,
) -> Option<(usize, usize, usize, usize)> {
    let x0 = left.floor().max(0.0) as usize;
    let x1 = (right.ceil().max(0.0) as usize).min(width);
    let y0 = top.floor().max(0.0) as usize;
    let y1 = (bottom.ceil().max(0.0) as usize).min(height);
    if x0 >= x1 || y0 >= y1 {
        None
    } else {
        Some((x0, x1, y0, y1))
    }
}

fn to_bgra(color: Vec<f32>, width: usize, height: usize) -> Frame {
    let mut bgra = vec![0u8; width * height * 4];
    for i in 0..width * height {
        let o = i * 3;
        let q = i * 4;
        bgra[q] = (color[o + 2].clamp(0.0, 1.0) * 255.0 + 0.5) as u8;
        bgra[q + 1] = (color[o + 1].clamp(0.0, 1.0) * 255.0 + 0.5) as u8;
        bgra[q + 2] = (color[o].clamp(0.0, 1.0) * 255.0 + 0.5) as u8;
        bgra[q + 3] = 255;
    }
    Frame {
        width,
        height,
        bgra,
    }
}
