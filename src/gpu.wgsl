// The GPU half of the rasterizer. The same picture as the CPU path in
// render.rs: orthographic view, spheres and capsules ray-traced per pixel,
// ribbons as plain triangles, then fog and ink as a screen-space pass.
//
// View space has +z toward the viewer, so a larger z is nearer.

struct U {
    rot0: vec4f,
    rot1: vec4f,
    rot2: vec4f,
    center: vec4f,
    // cx, cy, scale
    view: vec4f,
    // width, height, whether the ink pass finds its own edges
    size: vec4f,
    // zfar, depth range, zmin, span
    depth: vec4f,
    // page colour, fog strength
    page: vec4f,
    // jump, style (0 soft, 1 illustrative), rim, line weight
    params: vec4f,
}

@group(0) @binding(0) var<uniform> u: U;

const LIGHT = vec3f(-0.4082, 0.5585, 0.7226);
const HALF = vec3f(-0.2225, 0.3044, 0.9262);
const INK = vec3f(0.07, 0.07, 0.09);
// The z target's clear value, read as "nothing here".
const EMPTY = -1e30;

fn rot(v: vec3f) -> vec3f {
    return vec3f(dot(u.rot0.xyz, v), dot(u.rot1.xyz, v), dot(u.rot2.xyz, v));
}

fn to_view(p: vec3f) -> vec3f {
    return rot(p - u.center.xyz);
}

fn clip(v: vec3f) -> vec4f {
    let sx = u.view.x + v.x * u.view.z;
    let sy = u.view.y - v.y * u.view.z;
    return vec4f(sx / u.size.x * 2.0 - 1.0, 1.0 - sy / u.size.y * 2.0, to_depth(v.z), 1.0);
}

fn to_depth(z: f32) -> f32 {
    return clamp((u.depth.x - z) / u.depth.y, 0.0, 1.0);
}

// The view-space point under a pixel centre.
fn view_xy(frag: vec2f) -> vec2f {
    return vec2f((frag.x - u.view.x) / u.view.z, -(frag.y - u.view.y) / u.view.z);
}

fn shade(n: vec3f, base: vec3f) -> vec3f {
    let diff = max(dot(n, LIGHT), 0.0);
    let nh = max(dot(n, HALF), 0.0);
    let n2 = nh * nh;
    let n8 = n2 * n2 * n2 * n2;
    var k: f32;
    var spec: f32;
    if u.params.y < 0.5 {
        k = 0.34 + 0.70 * diff;
        spec = n8 * n8 * n8 * 0.28;
    } else {
        let banded = floor(diff * 4.0) / 3.0;
        let mixed = diff * 0.35 + min(banded, 1.0) * 0.65;
        k = 0.62 + 0.42 * mixed;
        spec = n8 * n8 * 0.10;
    }
    return base * k + vec3f(spec);
}

// -- spheres and capsules ----------------------------------------------------

struct Hit {
    ok: bool,
    z: f32,
    n: vec3f,
    // Along the capsule, 0 at a and 1 at b.
    f: f32,
}

fn sphere_hit(c: vec3f, r: f32, xy: vec2f) -> Hit {
    let d = xy - c.xy;
    let d2 = dot(d, d);
    let rr = r * r;
    if d2 > rr {
        return Hit(false, 0.0, vec3f(0.0), 0.0);
    }
    let dz = sqrt(rr - d2);
    return Hit(true, c.z + dz, vec3f(d.x, d.y, dz) / r, 0.0);
}

fn capsule_hit(a: vec3f, b: vec3f, r: f32, xy: vec2f) -> Hit {
    let d = b - a;
    let len = length(d);
    var ax = vec3f(0.0, 0.0, 1.0);
    if len > 1e-6 {
        ax = d / len;
    }
    let q = 1.0 - ax.z * ax.z;
    let w = xy - a.xy;
    if q > 1e-6 && len > 1e-6 {
        let k = w.x * ax.x + w.y * ax.y;
        let bq = -2.0 * k * ax.z;
        let cq = dot(w, w) - k * k - r * r;
        let disc = bq * bq - 4.0 * q * cq;
        if disc >= 0.0 {
            let cz = (-bq + sqrt(disc)) / (2.0 * q);
            let along = k + cz * ax.z;
            let s = along / len;
            if s >= 0.0 && s <= 1.0 {
                let n = (vec3f(w.x, w.y, cz) - ax * along) / r;
                return Hit(true, a.z + cz, n, s);
            }
        }
    }
    // Outside the cylinder's span, the caps take over.
    var best = sphere_hit(a, r, xy);
    var hb = sphere_hit(b, r, xy);
    hb.f = 1.0;
    if hb.ok && (!best.ok || hb.z > best.z) {
        best = hb;
    }
    return best;
}

struct ImpOut {
    @builtin(position) pos: vec4f,
    // View-space ends, with the radius and the kind (0 sphere, 1 capsule).
    @location(0) @interpolate(flat) a: vec4f,
    @location(1) @interpolate(flat) b: vec4f,
    @location(2) @interpolate(flat) ca: vec4f,
    @location(3) @interpolate(flat) cb: vec4f,
}

@vertex
fn imp_vs(
    @builtin(vertex_index) vi: u32,
    @location(0) ar: vec4f,
    @location(1) bk: vec4f,
    @location(2) ca: vec4f,
    @location(3) cb: vec4f,
) -> ImpOut {
    let a = to_view(ar.xyz);
    let b = to_view(bk.xyz);
    let r = ar.w;
    // A pixel of slack, so every pixel centre on the rim is reached.
    let pad = r + 1.5 / u.view.z;
    let d = b.xy - a.xy;
    let len = length(d);
    var dir = vec2f(1.0, 0.0);
    if len > 1e-6 {
        dir = d / len;
    }
    let perp = vec2f(-dir.y, dir.x);
    var corners = array<vec2f, 6>(
        vec2f(0.0, -1.0), vec2f(1.0, -1.0), vec2f(1.0, 1.0),
        vec2f(0.0, -1.0), vec2f(1.0, 1.0), vec2f(0.0, 1.0),
    );
    let c = corners[vi];
    let xy = mix(a.xy - dir * pad, b.xy + dir * pad, c.x) + perp * pad * c.y;

    var out: ImpOut;
    out.pos = clip(vec3f(xy, 0.0));
    out.pos.z = 0.5;
    out.a = vec4f(a, r);
    out.b = vec4f(b, bk.w);
    out.ca = ca;
    out.cb = cb;
    return out;
}

fn trace(in: ImpOut) -> Hit {
    let xy = view_xy(in.pos.xy);
    if in.b.w < 0.5 {
        return sphere_hit(in.a.xyz, in.a.w, xy);
    }
    return capsule_hit(in.a.xyz, in.b.xyz, in.a.w, xy);
}

struct Opaque {
    @location(0) color: vec4f,
    @builtin(frag_depth) depth: f32,
}

@fragment
fn imp_fs(in: ImpOut) -> Opaque {
    let h = trace(in);
    if !h.ok {
        discard;
    }
    let base = mix(in.ca.rgb, in.cb.rgb, h.f);
    return Opaque(vec4f(shade(h.n, base), 1.0), to_depth(h.z));
}

struct Blended {
    @location(0) color: vec4f,
    @builtin(frag_depth) depth: f32,
}

// Laid over the finished opaque frame, back to front, with the same fog and,
// in the illustrative style, an ink rim on spheres.
@fragment
fn glass_fs(in: ImpOut) -> Blended {
    let h = trace(in);
    if !h.ok {
        discard;
    }
    var c = shade(h.n, mix(in.ca.rgb, in.cb.rgb, h.f));
    let d = clamp((h.z - u.depth.z) / u.depth.w, 0.0, 1.0);
    c += (u.page.rgb - c) * (1.0 - d) * u.page.w;
    var alpha = in.ca.w;
    let sphere = in.b.w < 0.5;
    let inner = max(in.a.w * u.view.z - u.params.z, 0.0);
    let from_centre = length(view_xy(in.pos.xy) - in.a.xy) * u.view.z;
    if sphere && from_centre > inner {
        if u.params.y < 0.5 {
            c *= 0.8;
        } else {
            c = INK;
            alpha = max(alpha, 0.45);
        }
    }
    return Blended(vec4f(c, alpha), to_depth(h.z));
}

// -- ribbons -------------------------------------------------------------

struct TriOut {
    @builtin(position) pos: vec4f,
    @location(0) n: vec3f,
    @location(1) @interpolate(flat) c: vec3f,
}

@vertex
fn tri_vs(@location(0) p: vec3f, @location(1) n: vec3f, @location(2) c: vec3f) -> TriOut {
    let v = to_view(p);
    var out: TriOut;
    out.pos = clip(v);
    out.n = rot(n);
    out.c = c;
    return out;
}

// Depth comes from the rasterizer here, so the hardware can reject hidden
// fragments before shading them.
@fragment
fn tri_fs(in: TriOut) -> @location(0) vec4f {
    let len = length(in.n);
    if len < 1e-6 {
        discard;
    }
    var n = in.n / len;
    // A ribbon is seen from both sides; light whichever faces us.
    if n.z < 0.0 {
        n = -n;
    }
    return vec4f(shade(n, in.c), 1.0);
}

// -- fog and ink ---------------------------------------------------------

@group(0) @binding(1) var color_tex: texture_2d<f32>;
@group(0) @binding(2) var depth_tex: texture_depth_2d;
// Where the ink goes, before it is thickened.
@group(0) @binding(3) var edge_tex: texture_2d<f32>;

@vertex
fn post_vs(@builtin(vertex_index) vi: u32) -> @builtin(position) vec4f {
    let x = f32((vi << 1u) & 2u);
    let y = f32(vi & 2u);
    return vec4f(x * 2.0 - 1.0, y * 2.0 - 1.0, 0.0, 1.0);
}

fn inside(p: vec2i) -> bool {
    return p.x >= 0 && p.y >= 0 && p.x < i32(u.size.x) && p.y < i32(u.size.y);
}

// The view-space depth at a pixel, or EMPTY where nothing was drawn. The
// depth buffer holds it mapped linearly into 0..1, and 1 where it is clear.
fn z_at(p: vec2i) -> f32 {
    let d = textureLoad(depth_tex, p, 0);
    if d >= 1.0 {
        return -3e38;
    }
    return u.depth.x - d * u.depth.y;
}

const OFFSETS = array<vec2i, 4>(vec2i(-1, 0), vec2i(1, 0), vec2i(0, -1), vec2i(0, 1));

// A silhouette, inked from both sides, or a jump in depth.
@fragment
fn edge_fs(@builtin(position) pos: vec4f) -> @location(0) vec4f {
    let p = vec2i(pos.xy);
    let z = z_at(p);
    for (var i = 0; i < 4; i++) {
        let q = p + OFFSETS[i];
        if !inside(q) {
            continue;
        }
        let nz = z_at(q);
        let silhouette = (z < EMPTY) != (nz < EMPTY);
        let overlap = z >= EMPTY && nz >= EMPTY && nz - z > u.params.x;
        if silhouette || overlap {
            return vec4f(1.0);
        }
    }
    return vec4f(0.0);
}

// Outside the frame, which a neighbour test skips.
const OUTSIDE = 3e38;

fn z_or_outside(p: vec2i) -> f32 {
    if inside(p) {
        return z_at(p);
    }
    return OUTSIDE;
}

// Whether a pixel of depth `z` is on a line, given one of its neighbours.
fn crosses(z: f32, nz: f32) -> bool {
    if nz == OUTSIDE {
        return false;
    }
    let silhouette = (z < EMPTY) != (nz < EMPTY);
    let overlap = z >= EMPTY && nz >= EMPTY && nz - z > u.params.x;
    return silhouette || overlap;
}

fn edge4(z: f32, a: f32, b: f32, c: f32, d: f32) -> bool {
    return z != OUTSIDE && (crosses(z, a) || crosses(z, b) || crosses(z, c) || crosses(z, d));
}

@fragment
fn post_fs(@builtin(position) pos: vec4f) -> @location(0) vec4f {
    let p = vec2i(pos.xy);
    // Lines of the usual weights, one or two pixels, are found here from the
    // depth around the pixel, which saves a pass; heavier ones, in a large
    // export, are read from `edge_fs`.
    let w = i32(u.params.w) - 1;
    if u.params.y >= 0.5 && u.size.z > 0.5 {
        // The pixel, its four neighbours, and for a two-pixel line theirs:
        // each depth read once, held in registers.
        let c = z_at(p);
        let l = z_or_outside(p + vec2i(-1, 0));
        let r = z_or_outside(p + vec2i(1, 0));
        let t = z_or_outside(p + vec2i(0, -1));
        let b = z_or_outside(p + vec2i(0, 1));
        if edge4(c, l, r, t, b) {
            return vec4f(INK, 1.0);
        }
        if w == 1 {
            let tl = z_or_outside(p + vec2i(-1, -1));
            let tr = z_or_outside(p + vec2i(1, -1));
            let bl = z_or_outside(p + vec2i(-1, 1));
            let br = z_or_outside(p + vec2i(1, 1));
            if edge4(l, z_or_outside(p + vec2i(-2, 0)), c, tl, bl)
                || edge4(r, z_or_outside(p + vec2i(2, 0)), c, tr, br)
                || edge4(t, z_or_outside(p + vec2i(0, -2)), c, tl, tr)
                || edge4(b, z_or_outside(p + vec2i(0, 2)), c, bl, br) {
                return vec4f(INK, 1.0);
            }
        }
    } else if u.params.y >= 0.5 {
        // The line, thickened to the frame's weight.
        for (var dy = -w; dy <= w; dy++) {
            let reach = w - abs(dy);
            for (var dx = -reach; dx <= reach; dx++) {
                let q = p + vec2i(dx, dy);
                if inside(q) && textureLoad(edge_tex, q, 0).r > 0.5 {
                    return vec4f(INK, 1.0);
                }
            }
        }
    }

    let z = z_at(p);
    if z < EMPTY {
        return vec4f(u.page.rgb, 1.0);
    }
    var c = textureLoad(color_tex, p, 0).rgb;
    let d = clamp((z - u.depth.z) / u.depth.w, 0.0, 1.0);
    c += (u.page.rgb - c) * (1.0 - d) * u.page.w;

    if u.params.y < 0.5 {
        var behind = false;
        var open = false;
        for (var i = 0; i < 4; i++) {
            let q = p + OFFSETS[i];
            if !inside(q) {
                continue;
            }
            let nz = z_at(q);
            if nz < EMPTY {
                open = true;
            } else if nz - z > u.params.x {
                behind = true;
            }
        }
        if behind {
            c *= 0.62;
        } else if open {
            c *= 0.80;
        }
    }
    return vec4f(c, 1.0);
}
