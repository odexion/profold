//! Ray-traced pictures, for figures. The same scene the viewport draws,
//! seen the same way, but each point lit by rays: shadows from the light,
//! softened by its size, and ambient occlusion, the sky a point can see,
//! so a pocket goes dark and a cleft reads as deep. The frame is drawn
//! several times over and averaged down, so edges are smooth.
//!
//! The view is the rasterizer's: points warped as in `Camera::warp` and seen
//! straight on, so every primary ray runs down the view axis and finds what
//! the viewport shows. Fog, ink and the soft rim are then laid over by
//! `render::finish`, so each style keeps its own look.

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use rayon::prelude::*;

use crate::render::{self, Camera, Frame, Scene, Style};

type V = [f32; 3];

fn add(a: V, b: V) -> V {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn sub(a: V, b: V) -> V {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn mul(a: V, k: f32) -> V {
    [a[0] * k, a[1] * k, a[2] * k]
}
fn dot(a: V, b: V) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn unit(a: V) -> V {
    let l = dot(a, a).sqrt();
    if l > 1e-12 {
        mul(a, 1.0 / l)
    } else {
        [0.0, 0.0, 1.0]
    }
}
fn mix(a: V, b: V, t: f32) -> V {
    add(a, mul(sub(b, a), t))
}

/// How a picture is made, and how its making is followed and stopped.
#[derive(Default)]
pub struct Control<'a> {
    /// Set to give up; the trace then returns nothing.
    pub cancel: Option<&'a AtomicBool>,
    /// Counts rows as they are finished, out of `Traced::rows`.
    pub done: Option<&'a AtomicUsize>,
}

#[derive(Clone, Copy)]
enum Shape {
    Sphere(V, f32),
    Capsule(V, V, f32),
    /// Corners, and their normals.
    Tri([V; 3], [V; 3]),
}

#[derive(Clone, Copy)]
struct Prim {
    shape: Shape,
    ca: V,
    cb: V,
    alpha: f32,
}

impl Prim {
    fn bounds(&self) -> (V, V) {
        match self.shape {
            Shape::Sphere(c, r) => (sub(c, [r; 3]), add(c, [r; 3])),
            Shape::Capsule(a, b, r) => {
                let lo = std::array::from_fn(|k| a[k].min(b[k]) - r);
                let hi = std::array::from_fn(|k| a[k].max(b[k]) + r);
                (lo, hi)
            }
            Shape::Tri(v, _) => {
                let lo = std::array::from_fn(|k| v[0][k].min(v[1][k]).min(v[2][k]));
                let hi = std::array::from_fn(|k| v[0][k].max(v[1][k]).max(v[2][k]));
                (lo, hi)
            }
        }
    }

    fn radius(&self) -> Option<f32> {
        match self.shape {
            Shape::Sphere(_, r) | Shape::Capsule(_, _, r) => Some(r),
            Shape::Tri(..) => None,
        }
    }
}

/// Where a vertical line through `(x, y)` enters a primitive and leaves it,
/// nearest first, with the normal where it enters and how far along a
/// capsule that is.
struct Span {
    front: f32,
    back: f32,
    n: V,
    f: f32,
}

fn column(p: &Prim, x: f32, y: f32) -> Option<Span> {
    match p.shape {
        Shape::Sphere(c, r) => {
            let (dx, dy) = (x - c[0], y - c[1]);
            let d2 = dx * dx + dy * dy;
            if d2 > r * r {
                return None;
            }
            let dz = (r * r - d2).sqrt();
            Some(Span {
                front: c[2] + dz,
                back: c[2] - dz,
                n: [dx / r, dy / r, dz / r],
                f: 0.0,
            })
        }
        Shape::Capsule(a, b, r) => {
            // A capsule is convex, so the line leaves it where a ray sent up
            // from below first meets it.
            let top = a[2].max(b[2]) + r + 1.0;
            let bottom = a[2].min(b[2]) - r - 1.0;
            let t = capsule_hit([x, y, top], [0.0, 0.0, -1.0], a, b, r)?;
            let u =
                capsule_hit([x, y, bottom], [0.0, 0.0, 1.0], a, b, r).unwrap_or(top - bottom - t);
            let front = top - t;
            let (n, f) = capsule_normal([x, y, front], a, b, r);
            Some(Span {
                front,
                back: bottom + u,
                n,
                f,
            })
        }
        Shape::Tri(v, nv) => {
            let d = (v[1][1] - v[2][1]) * (v[0][0] - v[2][0])
                + (v[2][0] - v[1][0]) * (v[0][1] - v[2][1]);
            if d.abs() < 1e-12 {
                return None;
            }
            let l0 =
                ((v[1][1] - v[2][1]) * (x - v[2][0]) + (v[2][0] - v[1][0]) * (y - v[2][1])) / d;
            let l1 =
                ((v[2][1] - v[0][1]) * (x - v[2][0]) + (v[0][0] - v[2][0]) * (y - v[2][1])) / d;
            let l2 = 1.0 - l0 - l1;
            if l0 < 0.0 || l1 < 0.0 || l2 < 0.0 {
                return None;
            }
            let z = l0 * v[0][2] + l1 * v[1][2] + l2 * v[2][2];
            let mut n = unit(add(add(mul(nv[0], l0), mul(nv[1], l1)), mul(nv[2], l2)));
            // Seen from both sides: lit on whichever faces us.
            if n[2] < 0.0 {
                n = mul(n, -1.0);
            }
            Some(Span {
                front: z,
                back: z,
                n,
                f: 0.0,
            })
        }
    }
}

fn capsule_normal(p: V, a: V, b: V, r: f32) -> (V, f32) {
    let ba = sub(b, a);
    let bb = dot(ba, ba);
    let f = if bb > 1e-12 {
        (dot(sub(p, a), ba) / bb).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (mul(sub(p, add(a, mul(ba, f))), 1.0 / r), f)
}

fn sphere_hit(o: V, d: V, c: V, r: f32) -> Option<f32> {
    let oc = sub(o, c);
    let b = dot(oc, d);
    let h = b * b - (dot(oc, oc) - r * r);
    if h < 0.0 {
        return None;
    }
    Some(-b - h.sqrt())
}

/// Where a ray from `o` along unit `d` first enters a capsule from outside,
/// after Íñigo Quílez.
fn capsule_hit(o: V, d: V, a: V, b: V, r: f32) -> Option<f32> {
    let ba = sub(b, a);
    let oa = sub(o, a);
    let baba = dot(ba, ba);
    let bard = dot(ba, d);
    let baoa = dot(ba, oa);
    let k2 = baba - bard * bard;
    if k2 > 1e-6 * baba.max(1e-12) {
        let k1 = baba * dot(oa, d) - baoa * bard;
        let k0 = baba * dot(oa, oa) - baoa * baoa - r * r * baba;
        let h = k1 * k1 - k2 * k0;
        if h < 0.0 {
            return None;
        }
        let t = (-k1 - h.sqrt()) / k2;
        let y = baoa + t * bard;
        if y > 0.0 && y < baba {
            return Some(t);
        }
        let end = if y <= 0.0 { a } else { b };
        return sphere_hit(o, d, end, r);
    }
    // Along the axis: the nearer cap.
    match (sphere_hit(o, d, a, r), sphere_hit(o, d, b, r)) {
        (Some(s), Some(t)) => Some(s.min(t)),
        (s, t) => s.or(t),
    }
}

fn tri_hit(o: V, d: V, v: &[V; 3]) -> Option<f32> {
    let e1 = sub(v[1], v[0]);
    let e2 = sub(v[2], v[0]);
    let p = cross(d, e2);
    let det = dot(e1, p);
    if det.abs() < 1e-12 {
        return None;
    }
    let inv = 1.0 / det;
    let s = sub(o, v[0]);
    let u = dot(s, p) * inv;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = cross(s, e1);
    let w = dot(d, q) * inv;
    if w < 0.0 || u + w > 1.0 {
        return None;
    }
    Some(dot(e2, q) * inv)
}

fn ray_hit(p: &Prim, o: V, d: V) -> Option<f32> {
    match p.shape {
        Shape::Sphere(c, r) => sphere_hit(o, d, c, r),
        Shape::Capsule(a, b, r) => capsule_hit(o, d, a, b, r),
        Shape::Tri(v, _) => tri_hit(o, d, &v),
    }
}

/// A box in the tree: an inner one's children are at `first` and the one
/// after; a leaf holds `count` primitives from `first`.
struct Node {
    lo: V,
    hi: V,
    first: u32,
    count: u32,
}

struct Bvh {
    nodes: Vec<Node>,
    prims: Vec<Prim>,
}

const LEAF: usize = 4;

impl Bvh {
    fn new(mut prims: Vec<Prim>) -> Self {
        let mut nodes = Vec::with_capacity(prims.len().div_ceil(LEAF) * 2 + 1);
        if !prims.is_empty() {
            let boxes: Vec<(V, V)> = prims.iter().map(Prim::bounds).collect();
            let mut order: Vec<u32> = (0..prims.len() as u32).collect();
            nodes.push(Node {
                lo: [0.0; 3],
                hi: [0.0; 3],
                first: 0,
                count: 0,
            });
            Self::split(&mut nodes, 0, &boxes, &mut order, 0);
            prims = order.iter().map(|&i| prims[i as usize]).collect();
        }
        Bvh { nodes, prims }
    }

    /// Fill node `at` with `order`, which starts at `first` in the whole,
    /// halving at the median along its longest spread of centres.
    fn split(nodes: &mut Vec<Node>, at: usize, boxes: &[(V, V)], order: &mut [u32], first: usize) {
        let mut lo = [f32::MAX; 3];
        let mut hi = [f32::MIN; 3];
        let mut clo = [f32::MAX; 3];
        let mut chi = [f32::MIN; 3];
        for &i in order.iter() {
            let (l, h) = boxes[i as usize];
            for k in 0..3 {
                lo[k] = lo[k].min(l[k]);
                hi[k] = hi[k].max(h[k]);
                let c = (l[k] + h[k]) * 0.5;
                clo[k] = clo[k].min(c);
                chi[k] = chi[k].max(c);
            }
        }
        nodes[at].lo = lo;
        nodes[at].hi = hi;
        let spread = sub(chi, clo);
        if order.len() <= LEAF || spread.iter().all(|&s| s <= 0.0) {
            nodes[at].first = first as u32;
            nodes[at].count = order.len() as u32;
            return;
        }
        let axis = if spread[0] >= spread[1] && spread[0] >= spread[2] {
            0
        } else if spread[1] >= spread[2] {
            1
        } else {
            2
        };
        let mid = order.len() / 2;
        let centre = |i: &u32| {
            let (l, h) = boxes[*i as usize];
            l[axis] + h[axis]
        };
        order.select_nth_unstable_by(mid, |a, b| centre(a).total_cmp(&centre(b)));
        let left = nodes.len();
        nodes[at].first = left as u32;
        nodes[at].count = 0;
        for _ in 0..2 {
            nodes.push(Node {
                lo: [0.0; 3],
                hi: [0.0; 3],
                first: 0,
                count: 0,
            });
        }
        let (a, b) = order.split_at_mut(mid);
        Self::split(nodes, left, boxes, a, first);
        Self::split(nodes, left + 1, boxes, b, first + mid);
    }

    /// Every primitive whose box the vertical line through `(x, y)` passes
    /// and that reaches above `floor`, handed to `f`, which returns a new
    /// floor: the line is searched top down and stops below it.
    fn down(&self, x: f32, y: f32, mut floor: f32, mut f: impl FnMut(&Prim, f32) -> f32) {
        if self.nodes.is_empty() {
            return;
        }
        let mut stack = [0u32; 64];
        let mut top = 1;
        while top > 0 {
            top -= 1;
            let node = &self.nodes[stack[top] as usize];
            if x < node.lo[0]
                || x > node.hi[0]
                || y < node.lo[1]
                || y > node.hi[1]
                || node.hi[2] < floor
            {
                continue;
            }
            if node.count > 0 {
                let first = node.first as usize;
                for p in &self.prims[first..first + node.count as usize] {
                    floor = f(p, floor);
                }
            } else {
                // The higher child last, so it is searched first and raises
                // the floor under the other.
                let (a, b) = (node.first, node.first + 1);
                let (low, high) = if self.nodes[a as usize].hi[2] > self.nodes[b as usize].hi[2] {
                    (b, a)
                } else {
                    (a, b)
                };
                stack[top] = low;
                stack[top + 1] = high;
                top += 2;
            }
        }
    }

    /// Whether a ray from `o` along unit `d` meets anything within `reach`
    /// whose point lies between the slab's depths.
    fn blocked(&self, o: V, d: V, reach: f32, near: f32, far: f32) -> bool {
        if self.nodes.is_empty() {
            return false;
        }
        let inv = d.map(|c| 1.0 / c);
        let mut stack = [0u32; 64];
        let mut top = 1;
        while top > 0 {
            top -= 1;
            let node = &self.nodes[stack[top] as usize];
            let mut t0: f32 = 0.0;
            let mut t1 = reach;
            for k in 0..3 {
                let a = (node.lo[k] - o[k]) * inv[k];
                let b = (node.hi[k] - o[k]) * inv[k];
                t0 = t0.max(a.min(b));
                t1 = t1.min(a.max(b));
            }
            if t0 > t1 {
                continue;
            }
            if node.count > 0 {
                let first = node.first as usize;
                for p in &self.prims[first..first + node.count as usize] {
                    if let Some(t) = ray_hit(p, o, d) {
                        let z = o[2] + d[2] * t;
                        if t > 0.0 && t < reach && z <= near && z >= far {
                            return true;
                        }
                    }
                }
            } else {
                stack[top] = node.first;
                stack[top + 1] = node.first + 1;
                top += 2;
            }
        }
        false
    }
}

/// What a pixel's line of sight comes to rest on.
struct Surface {
    z: f32,
    n: V,
    base: V,
}

/// The nearest opaque surface down the line through `(x, y)`, cut by the
/// slab as the rasterizer cuts it: a solid the near plane passes through
/// shows a flat face there.
fn first_surface(bvh: &Bvh, x: f32, y: f32, near: f32, far: f32) -> Option<Surface> {
    let mut best: Option<(f32, V, V)> = None;
    bvh.down(x, y, far, |p, floor| {
        let Some(s) = column(p, x, y) else {
            return floor;
        };
        let base = mix(p.ca, p.cb, s.f);
        let found = if s.front > near {
            (s.back <= near && !matches!(p.shape, Shape::Tri(..)))
                .then(|| (near, render::CAP, base.map(|c| c * render::CAP_SHADE)))
        } else {
            (s.front >= far).then_some((s.front, s.n, base))
        };
        match found {
            Some(hit) if best.is_none_or(|b| hit.0 > b.0) => {
                best = Some(hit);
                hit.0
            }
            _ => floor,
        }
    });
    best.map(|(z, n, base)| Surface { z, n, base })
}

/// A hash of a pixel, for decorrelating its samples from its neighbours'.
fn hash(x: u32, y: u32, salt: u32) -> u32 {
    let mut h =
        x.wrapping_mul(0x8da6_b343) ^ y.wrapping_mul(0xd816_3841) ^ salt.wrapping_mul(0xcb1a_b31f);
    h ^= h >> 16;
    h = h.wrapping_mul(0x7feb_352d);
    h ^= h >> 15;
    h = h.wrapping_mul(0x846c_a68b);
    h ^ (h >> 16)
}

/// The `i`th of a well-spread sequence over the unit square, shifted by a
/// pixel's own offset: Roberts' R2.
fn r2(i: u32, shift: u32) -> (f32, f32) {
    let s = (shift & 0xffff) as f32 / 65536.0;
    let t = (shift >> 16) as f32 / 65536.0;
    let a = (0.5 + i as f32 * 0.754_877_7 + s).fract();
    let b = (0.5 + i as f32 * 0.569_840_3 + t).fract();
    (a, b)
}

/// Two directions square to `n` and to each other.
fn frame(n: V) -> (V, V) {
    let up = if n[2].abs() < 0.9 {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let t = unit(cross(up, n));
    (t, cross(n, t))
}

/// How much the shadows and the occlusion darken, by style: a drawn figure
/// takes less of either, so the ink still leads.
fn strengths(style: Style) -> (f32, f32) {
    match style {
        Style::Studio => (0.85, 0.75),
        Style::Soft => (0.7, 0.6),
        Style::Illustrative => (0.55, 0.5),
        Style::Flat => (0.5, 0.35),
    }
}

/// The light's half-width, in radians: a lamp rather than a point, so the
/// shadows have soft edges.
const LAMP: f32 = 0.06;

/// The traced frame before it is laid on a page: shaded colour and depth at
/// `ss` times the size asked for each way, and the see-through primitives
/// still to be laid over it.
pub struct Traced {
    width: usize,
    height: usize,
    ss: usize,
    color: Vec<f32>,
    depth: Vec<f32>,
    glass: Bvh,
    style: Style,
    cx: f32,
    cy: f32,
    scale: f32,
    near: f32,
    far: f32,
    zmin: f32,
    span: f32,
}

/// The most pixels traced in one frame, oversampling included: about 400 MB
/// of buffers. A large picture is traced fewer times over per pixel, down to
/// once, which at its size is fine enough.
const BUDGET: usize = 25_000_000;

/// How many times over each way a frame of this size is traced.
pub fn oversampling(width: usize, height: usize) -> usize {
    (1..=3)
        .rev()
        .find(|s| width * height * s * s <= BUDGET)
        .unwrap_or(1)
}

/// The rows a frame of this size traces, for following its progress.
pub fn rows(height: usize, width: usize) -> usize {
    height * oversampling(width, height)
}

/// Trace `scene` as `cam` sees it into a frame `width` by `height`. Nothing
/// if `control` cancels it.
pub fn trace(
    scene: &Scene,
    cam: &Camera,
    width: usize,
    height: usize,
    style: Style,
    control: &Control,
) -> Option<Traced> {
    let ss = oversampling(width, height);
    let (w, h) = (width * ss, height * ss);
    let m = cam.basis();
    let to_view = |p: V| {
        let q = sub(p, cam.center);
        [dot(m[0], q), dot(m[1], q), dot(m[2], q)]
    };
    let rot = |v: V| [dot(m[0], v), dot(m[1], v), dot(m[2], v)];

    let mut opaque = Vec::new();
    let mut glass = Vec::new();
    for s in &scene.spheres {
        if s.alpha <= 0.0 || s.r <= 0.0 {
            continue;
        }
        let (c, k) = cam.warp(to_view(s.c));
        let p = Prim {
            shape: Shape::Sphere(c, s.r * k),
            ca: s.color,
            cb: s.color,
            alpha: s.alpha,
        };
        if s.alpha < render::OPAQUE {
            glass.push(p)
        } else {
            opaque.push(p)
        }
    }
    for c in &scene.capsules {
        if c.alpha <= 0.0 {
            continue;
        }
        let (a, ka) = cam.warp(to_view(c.a));
        let (b, kb) = cam.warp(to_view(c.b));
        let p = Prim {
            shape: Shape::Capsule(a, b, c.r * (ka + kb) * 0.5),
            ca: c.ca,
            cb: c.cb,
            alpha: c.alpha,
        };
        if c.alpha < render::OPAQUE {
            glass.push(p)
        } else {
            opaque.push(p)
        }
    }
    for t in &scene.tris {
        let v = t.p.map(|p| cam.warp(to_view(p)).0);
        opaque.push(Prim {
            shape: Shape::Tri(v, t.n.map(rot)),
            ca: t.c,
            cb: t.c,
            alpha: 1.0,
        });
    }

    let (mut zmin, mut zmax) = (f32::MAX, f32::MIN);
    for p in opaque.iter().chain(&glass) {
        let (lo, hi) = p.bounds();
        zmin = zmin.min(lo[2]);
        zmax = zmax.max(hi[2]);
    }
    let span = (zmax - zmin).max(1e-3);

    // How far a point looks for what shuts it in: a share of the whole, but
    // never less than a few atoms across, so a lone molecule still shades
    // where its atoms meet.
    let mut radii: Vec<f32> = opaque.iter().filter_map(Prim::radius).collect();
    let typical = if radii.is_empty() {
        1.5
    } else {
        let mid = radii.len() / 2;
        *radii.select_nth_unstable_by(mid, f32::total_cmp).1
    };
    let reach = (cam.radius * 0.2).max(typical * 4.0);
    let eps = (cam.radius * 1e-4).max(1e-4);
    let light = render::LIGHT;
    let (lt, lb) = frame(light);
    let (shadow_k, sky_k) = strengths(style);
    let per = (ss * ss) as u32;
    // Rays for each pixel of the picture, shared among its samples.
    let shadow_rays = 36u32.div_ceil(per);
    let sky_rays = 96u32.div_ceil(per);

    let (near, far) = cam.depths();
    let scale = (w.min(h) as f32) / (2.2 * cam.radius) * cam.zoom;
    let cx = w as f32 * 0.5 + cam.pan[0] * ss as f32;
    let cy = h as f32 * 0.5 + cam.pan[1] * ss as f32;
    let inv = 1.0 / scale;

    let opaque = Bvh::new(opaque);
    let glass = Bvh::new(glass);
    let mut color = vec![0.0f32; w * h * 3];
    let mut depth = vec![f32::NEG_INFINITY; w * h];
    let cancelled = || control.cancel.is_some_and(|c| c.load(Ordering::Relaxed));

    color
        .par_chunks_mut(w * 3)
        .zip(depth.par_chunks_mut(w))
        .enumerate()
        .for_each(|(py, (row, zrow))| {
            if cancelled() {
                return;
            }
            let y = -(py as f32 + 0.5 - cy) * inv;
            for px in 0..w {
                let x = (px as f32 + 0.5 - cx) * inv;
                let Some(s) = first_surface(&opaque, x, y, near, far) else {
                    continue;
                };
                let o = add([x, y, s.z], mul(s.n, eps));
                // One sequence for the pixel of the picture this sample is
                // part of, each sample taking its own stretch of it, so the
                // samples averaged together spread evenly between them.
                let (ox, oy) = ((px / ss) as u32, (py / ss) as u32);
                let part = ((py % ss) * ss + px % ss) as u32;
                let seed = hash(ox, oy, 0);

                let lit = if dot(s.n, light) > 0.0 {
                    let mut through = 0;
                    for i in 0..shadow_rays {
                        let (u, v) = r2(part * shadow_rays + i, seed);
                        let (rr, phi) = (u.sqrt() * LAMP, v * std::f32::consts::TAU);
                        let d = unit(add(
                            light,
                            add(mul(lt, rr * phi.cos()), mul(lb, rr * phi.sin())),
                        ));
                        if !opaque.blocked(o, d, f32::MAX, near, far) {
                            through += 1;
                        }
                    }
                    1.0 - shadow_k * (1.0 - through as f32 / shadow_rays as f32)
                } else {
                    1.0
                };

                let (t, b) = frame(s.n);
                let seed = hash(ox, oy, 1);
                let mut shut = 0;
                for i in 0..sky_rays {
                    let (u, v) = r2(part * sky_rays + i, seed);
                    let (rr, phi) = (u.sqrt(), v * std::f32::consts::TAU);
                    let up = (1.0 - u).max(0.0).sqrt();
                    let d = add(
                        add(mul(t, rr * phi.cos()), mul(b, rr * phi.sin())),
                        mul(s.n, up),
                    );
                    if opaque.blocked(o, d, reach, near, far) {
                        shut += 1;
                    }
                }
                let open = 1.0 - sky_k * shut as f32 / sky_rays as f32;

                let c = render::shade_lit(s.n, s.base, style, lit, open);
                row[px * 3..px * 3 + 3].copy_from_slice(&c);
                zrow[px] = s.z;
            }
            if let Some(done) = control.done {
                done.fetch_add(1, Ordering::Relaxed);
            }
        });
    if cancelled() {
        return None;
    }

    Some(Traced {
        width,
        height,
        ss,
        color,
        depth,
        glass,
        style,
        cx,
        cy,
        scale,
        near,
        far,
        zmin,
        span,
    })
}

impl Traced {
    /// The frame on a page of colour `page`: fog, ink and the see-through
    /// layers laid on, then averaged down to the size asked for.
    pub fn develop(&self, page: V) -> Frame {
        let (w, h, ss) = (self.width * self.ss, self.height * self.ss, self.ss);
        let style = self.style;
        let mut color = self.color.clone();
        color.par_chunks_mut(3).zip(&self.depth).for_each(|(c, z)| {
            if *z == f32::NEG_INFINITY {
                c.copy_from_slice(&page);
            }
        });
        render::finish(
            &mut color,
            &self.depth,
            w,
            h,
            page,
            style,
            self.zmin,
            self.span,
            self.scale,
            false,
            ss,
        );

        if !self.glass.prims.is_empty() {
            let fog = style.fog();
            let rim = ((w.min(h) / ss / 700).max(1) * ss) as f32 + 0.5;
            let inv = 1.0 / self.scale;
            color
                .par_chunks_mut(w * 3)
                .enumerate()
                .for_each_init(Vec::new, |layers, (py, row)| {
                    let y = -(py as f32 + 0.5 - self.cy) * inv;
                    for px in 0..w {
                        let x = (px as f32 + 0.5 - self.cx) * inv;
                        let solid = self.depth[py * w + px];
                        layers.clear();
                        self.glass.down(x, y, solid.max(self.far), |p, floor| {
                            if let Some(s) = column(p, x, y) {
                                if s.front > solid && s.front <= self.near && s.front >= self.far {
                                    let edge = match p.shape {
                                        Shape::Sphere(c, r) => {
                                            let d = ((x - c[0]).powi(2) + (y - c[1]).powi(2))
                                                .sqrt()
                                                * self.scale;
                                            d > (r * self.scale - rim).max(0.0)
                                        }
                                        _ => false,
                                    };
                                    layers.push((
                                        s.front,
                                        s.n,
                                        mix(p.ca, p.cb, s.f),
                                        p.alpha,
                                        edge,
                                    ));
                                }
                            }
                            floor
                        });
                        layers.sort_by(|a, b| a.0.total_cmp(&b.0));
                        let out = &mut row[px * 3..px * 3 + 3];
                        for &(z, n, base, alpha, edge) in layers.iter() {
                            let mut c = render::shade_lit(n, base, style, 1.0, 1.0);
                            let d = ((z - self.zmin) / self.span).clamp(0.0, 1.0);
                            let wash = (1.0 - d) * fog;
                            for k in 0..3 {
                                c[k] += (page[k] - c[k]) * wash;
                            }
                            let (c, alpha) = match (edge, style) {
                                (true, Style::Illustrative | Style::Flat) => {
                                    ([0.07, 0.07, 0.09], alpha.max(0.45))
                                }
                                (true, Style::Soft) => (c.map(|x| x * 0.8), alpha),
                                _ => (c, alpha),
                            };
                            for k in 0..3 {
                                out[k] = c[k].min(1.0) * alpha + out[k] * (1.0 - alpha);
                            }
                        }
                    }
                });
        }

        // Averaged down, each sample clamped as the screen would show it.
        let (fw, fh) = (self.width, self.height);
        let mut small = vec![0.0f32; fw * fh * 3];
        let k = 1.0 / (ss * ss) as f32;
        small
            .par_chunks_mut(fw * 3)
            .enumerate()
            .for_each(|(y, row)| {
                for x in 0..fw {
                    let mut sum = [0.0f32; 3];
                    for sy in 0..ss {
                        for sx in 0..ss {
                            let i = ((y * ss + sy) * w + x * ss + sx) * 3;
                            for c in 0..3 {
                                sum[c] += color[i + c].clamp(0.0, 1.0);
                            }
                        }
                    }
                    for c in 0..3 {
                        row[x * 3 + c] = sum[c] * k;
                    }
                }
            });
        render::to_bgra(&small, fw, fh)
    }

    /// The frame on no page, as straight RGBA.
    pub fn transparent(&self) -> Vec<u8> {
        render::unmix(&self.develop([1.0; 3]), &self.develop([0.0; 3]))
    }
}

/// A traced frame on a page, or nothing if cancelled.
pub fn render(
    scene: &Scene,
    cam: &Camera,
    width: usize,
    height: usize,
    page: V,
    style: Style,
    control: &Control,
) -> Option<Frame> {
    trace(scene, cam, width, height, style, control).map(|t| t.develop(page))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::render::{build_scene, Appearance, Coloring, Palette, Repr};

    fn luma(f: &Frame, x: usize, y: usize) -> f32 {
        let p = &f.bgra[(y * f.width + x) * 4..][..3];
        (p[0] as f32 + p[1] as f32 + p[2] as f32) / 3.0
    }

    #[test]
    fn a_sphere_behind_another_sits_in_its_shadow() {
        // Two balls, the front one up and to the left of the other, the way
        // the light comes: the back one is darker where the front one shades
        // it than on its open side.
        let ball = |c: V| render::Sphere {
            c,
            r: 1.0,
            color: [0.8, 0.8, 0.8],
            alpha: 1.0,
        };
        let scene = Scene {
            spheres: vec![ball([-0.9, 1.2, 2.2]), ball([0.0, 0.0, 0.0])],
            ..Scene::default()
        };
        let cam = Camera {
            orientation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            zoom: 1.0,
            pan: [0.0, 0.0],
            center: [0.0, 0.0, 0.0],
            radius: 3.0,
            perspective: false,
            slab: None,
        };
        let page = [1.0; 3];
        let traced = render(
            &scene,
            &cam,
            200,
            200,
            page,
            Style::Soft,
            &Control::default(),
        )
        .unwrap();
        let flat = render::render_cpu(&scene, &cam, 200, 200, page, Style::Soft);
        // A point on the back ball just below and right of the front one,
        // against one the same way out on the far side.
        let scale = 200.0 / (2.2 * 3.0);
        let at = |x: f32, y: f32| ((100.0 + x * scale) as usize, (100.0 - y * scale) as usize);
        let (sx, sy) = at(-0.1, 0.35);
        let (ox, oy) = at(0.6, -0.6);
        let shaded = luma(&traced, sx, sy) / luma(&flat, sx, sy);
        let open = luma(&traced, ox, oy) / luma(&flat, ox, oy);
        assert!(shaded < open - 0.1, "shaded {shaded}, open {open}");
    }

    #[test]
    fn a_traced_protein_matches_the_viewport_in_outline() {
        let s = crate::formats::load_bytes(
            "1crn.pdb.gz",
            include_bytes!("../samples/1crn.pdb.gz").to_vec(),
        )
        .unwrap();
        for repr in [Repr::Cartoon, Repr::Spacefill, Repr::BallAndStick] {
            for style in Style::ALL {
                let scene = build_scene(
                    &s,
                    repr,
                    Coloring::Chain,
                    1.0,
                    style,
                    Palette::Spectrum,
                    false,
                    false,
                    &Appearance::default(),
                );
                for cut in [false, true] {
                    let mut cam = Camera::frame(&s);
                    cam.rotate(0.3, 0.2);
                    // In perspective through a slab, the cut faces and all.
                    if cut {
                        cam.perspective = true;
                        cam.slab = Some([0.4 * cam.radius, -0.3 * cam.radius]);
                    }
                    let page = [0.97, 0.96, 0.94];
                    let traced =
                        render(&scene, &cam, 160, 120, page, style, &Control::default()).unwrap();
                    let flat = render::render_cpu(&scene, &cam, 160, 120, page, style);
                    // The same pixels are covered, give or take the smoothed edge.
                    let bg = |f: &Frame, i: usize| {
                        let p = &f.bgra[i * 4..][..3];
                        p.iter()
                            .zip([page[2], page[1], page[0]])
                            .all(|(&a, b)| a.abs_diff((b * 255.0 + 0.5) as u8) <= 2)
                    };
                    let differ = (0..160 * 120)
                        .filter(|&i| bg(&traced, i) != bg(&flat, i))
                        .count();
                    assert!(
                        differ < 160 * 120 / 25,
                        "{repr:?} {style:?}: {differ} pixels differ in coverage"
                    );
                }
            }
        }
    }

    #[test]
    fn glass_blends_over_what_is_behind_it_as_in_the_viewport() {
        let ball = |c: V, color: V, alpha: f32| render::Sphere { c, r: 1.0, color, alpha };
        let scene = Scene {
            spheres: vec![ball([0.0, 0.0, 0.0], [0.2, 0.3, 0.9], 1.0), ball([0.3, 0.2, 2.0], [0.9, 0.3, 0.2], 0.5)],
            ..Scene::default()
        };
        let cam = Camera {
            orientation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            zoom: 1.0,
            pan: [0.0, 0.0],
            center: [0.0, 0.0, 0.0],
            radius: 3.0,
            perspective: false,
            slab: None,
        };
        for style in Style::ALL {
            let page = [0.97, 0.96, 0.94];
            let traced = render(&scene, &cam, 120, 120, page, style, &Control::default()).unwrap();
            let flat = render::render_cpu(&scene, &cam, 120, 120, page, style);
            // Through the glass, at the middle of both.
            let i = (60 * 120 + 62) * 4;
            for k in 0..3 {
                let (a, b) = (traced.bgra[i + k], flat.bgra[i + k]);
                assert!(a.abs_diff(b) <= 12, "{style:?}: {a} against {b}");
            }
        }
    }

    #[test]
    fn cancelling_stops_the_trace() {
        let s = crate::formats::load_bytes(
            "caffeine.sdf",
            include_bytes!("../samples/caffeine.sdf").to_vec(),
        )
        .unwrap();
        let scene = build_scene(
            &s,
            Repr::Spacefill,
            Coloring::Element,
            1.0,
            Style::Studio,
            Palette::Spectrum,
            false,
            true,
            &Appearance::default(),
        );
        let cam = Camera::frame(&s);
        let stop = AtomicBool::new(true);
        let control = Control {
            cancel: Some(&stop),
            done: None,
        };
        assert!(trace(&scene, &cam, 64, 64, Style::Studio, &control).is_none());
    }

    /// `PROFOLD_TRACE_DIR=<dir> cargo test --release trace_gallery -- --ignored --nocapture`
    /// writes each sample drawn both ways, and how long the trace took.
    #[test]
    #[ignore]
    fn trace_gallery() {
        let dir = std::path::PathBuf::from(
            std::env::var("PROFOLD_TRACE_DIR").expect("PROFOLD_TRACE_DIR"),
        );
        let file = std::env::var("PROFOLD_TRACE_FILE").ok();
        let mut cases: Vec<(String, Vec<u8>, Repr)> = vec![
            (
                "4hhb.pdb.gz".into(),
                include_bytes!("../samples/4hhb.pdb.gz").to_vec(),
                Repr::Cartoon,
            ),
            (
                "1crn.pdb.gz".into(),
                include_bytes!("../samples/1crn.pdb.gz").to_vec(),
                Repr::Spacefill,
            ),
            (
                "1ema.pdb.gz".into(),
                include_bytes!("../samples/1ema.pdb.gz").to_vec(),
                Repr::Spacefill,
            ),
            (
                "caffeine.sdf".into(),
                include_bytes!("../samples/caffeine.sdf").to_vec(),
                Repr::BallAndStick,
            ),
        ];
        if let Some(f) = file {
            cases = vec![(f.clone(), std::fs::read(&f).unwrap(), Repr::Spacefill)];
        }
        let page = [0.97, 0.96, 0.94];
        for (name, bytes, repr) in cases {
            let s = crate::formats::load_bytes(&name, bytes).unwrap();
            let stem = std::path::Path::new(&name)
                .file_name()
                .unwrap()
                .to_string_lossy()
                .split('.')
                .next()
                .unwrap()
                .to_string();
            for style in Style::ALL {
                let scene = build_scene(
                    &s,
                    repr,
                    Coloring::Chain,
                    1.0,
                    style,
                    Palette::Spectrum,
                    false,
                    true,
                    &Appearance::default(),
                );
                let mut cam = Camera::frame(&s);
                cam.rotate(0.3, 0.2);
                let (w, h) = (1200, 900);
                let t = std::time::Instant::now();
                let traced = render(&scene, &cam, w, h, page, style, &Control::default()).unwrap();
                eprintln!(
                    "{stem} {style:?}: {:.2?} at {}x",
                    t.elapsed(),
                    oversampling(w, h)
                );
                let flat = render::render_cpu(&scene, &cam, w, h, page, style);
                for (f, kind) in [(traced, "traced"), (flat, "flat")] {
                    let mut rgba = f.bgra;
                    for p in rgba.chunks_exact_mut(4) {
                        p.swap(0, 2);
                    }
                    let path = dir.join(format!(
                        "{stem}-{}-{kind}.png",
                        style.label().to_lowercase()
                    ));
                    image::save_buffer(path, &rgba, w as u32, h as u32, image::ColorType::Rgba8)
                        .unwrap();
                }
            }
        }
    }
}
