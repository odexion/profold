// The ray tracer's pixels on the graphics card: the same as the loop in
// trace.rs, step for step. For each sample, the first opaque surface down
// the view axis (cut by the slab as the rasterizer cuts it), then rays to
// the light and over the sky above it, through the same tree, with the same
// sequence of directions. Fog, ink, glass and the averaging down are left
// to the CPU, as they are cheap.
//
// View space has +z toward the viewer, so a larger z is nearer.

struct Prim {
    // Sphere: centre in v0, radius v1.w. Capsule: ends v0 and v1, radius
    // v1.w. Triangle: corners v0, v1, v2 and their normals n0, n1, n2.
    // The kind is v0.w: 0 sphere, 1 capsule, 2 triangle.
    v0: vec4f,
    v1: vec4f,
    v2: vec4f,
    n0: vec4f,
    n1: vec4f,
    n2: vec4f,
    ca: vec4f,
    cb: vec4f,
}

struct Node {
    lo: vec3f,
    // A leaf's first primitive, or an inner node's first child.
    first: u32,
    hi: vec3f,
    // Primitives in a leaf; none in an inner node.
    count: u32,
}

struct Params {
    // cx, cy, inv (view units per sample), eps
    view: vec4f,
    // near, far, reach, shadow strength
    cut: vec4f,
    // sky strength, style, unused, unused
    look: vec4f,
    light: vec4f,
    lt: vec4f,
    lb: vec4f,
    // samples across, first row of the band, rows in the band, oversampling
    size: vec4u,
    // shadow rays, sky rays, unused, unused
    rays: vec4u,
}

@group(0) @binding(0) var<uniform> p: Params;
@group(0) @binding(1) var<storage, read> nodes: array<Node>;
@group(0) @binding(2) var<storage, read> prims: array<Prim>;
// Colour, and depth in w: EMPTY where nothing was met.
@group(0) @binding(3) var<storage, read_write> out: array<vec4f>;

const EMPTY = -3.0e38;
const TAU = 6.28318530718;
const LAMP = 0.06;
const HALF = vec3f(-0.2225, 0.3044, 0.9262);
const CAP = vec3f(0.0, 0.0, 1.0);
const CAP_SHADE = 0.82;
const STACK = 144u;

fn unit(v: vec3f) -> vec3f {
    let l = length(v);
    if l > 1e-12 {
        return v / l;
    }
    return vec3f(0.0, 0.0, 1.0);
}

// -- hits --------------------------------------------------------------

// Where a ray first enters a sphere from outside: t, or a miss in y.
fn sphere_hit(o: vec3f, d: vec3f, c: vec3f, r: f32) -> vec2f {
    let oc = o - c;
    let b = dot(oc, d);
    let h = b * b - (dot(oc, oc) - r * r);
    if h < 0.0 {
        return vec2f(0.0, 0.0);
    }
    return vec2f(-b - sqrt(h), 1.0);
}

fn capsule_hit(o: vec3f, d: vec3f, a: vec3f, b: vec3f, r: f32) -> vec2f {
    let ba = b - a;
    let oa = o - a;
    let baba = dot(ba, ba);
    let bard = dot(ba, d);
    let baoa = dot(ba, oa);
    let k2 = baba - bard * bard;
    if k2 > 1e-6 * max(baba, 1e-12) {
        let k1 = baba * dot(oa, d) - baoa * bard;
        let k0 = baba * dot(oa, oa) - baoa * baoa - r * r * baba;
        let h = k1 * k1 - k2 * k0;
        if h < 0.0 {
            return vec2f(0.0, 0.0);
        }
        let t = (-k1 - sqrt(h)) / k2;
        let y = baoa + t * bard;
        if y > 0.0 && y < baba {
            return vec2f(t, 1.0);
        }
        if y <= 0.0 {
            return sphere_hit(o, d, a, r);
        }
        return sphere_hit(o, d, b, r);
    }
    // Along the axis: the nearer cap.
    let s = sphere_hit(o, d, a, r);
    let u = sphere_hit(o, d, b, r);
    if s.y > 0.5 && u.y > 0.5 {
        return vec2f(min(s.x, u.x), 1.0);
    }
    if s.y > 0.5 {
        return s;
    }
    return u;
}

fn tri_hit(o: vec3f, d: vec3f, v0: vec3f, v1: vec3f, v2: vec3f) -> vec2f {
    let e1 = v1 - v0;
    let e2 = v2 - v0;
    let q = cross(d, e2);
    let det = dot(e1, q);
    if abs(det) < 1e-12 {
        return vec2f(0.0, 0.0);
    }
    let inv = 1.0 / det;
    let s = o - v0;
    let u = dot(s, q) * inv;
    if u < 0.0 || u > 1.0 {
        return vec2f(0.0, 0.0);
    }
    let qq = cross(s, e1);
    let w = dot(d, qq) * inv;
    if w < 0.0 || u + w > 1.0 {
        return vec2f(0.0, 0.0);
    }
    return vec2f(dot(e2, qq) * inv, 1.0);
}

fn ray_hit(q: Prim, o: vec3f, d: vec3f) -> vec2f {
    let kind = q.v0.w;
    if kind < 0.5 {
        return sphere_hit(o, d, q.v0.xyz, q.v1.w);
    }
    if kind < 1.5 {
        return capsule_hit(o, d, q.v0.xyz, q.v1.xyz, q.v1.w);
    }
    return tri_hit(o, d, q.v0.xyz, q.v1.xyz, q.v2.xyz);
}

// -- the line of sight ---------------------------------------------------

struct Span {
    hit: bool,
    front: f32,
    back: f32,
    n: vec3f,
    f: f32,
}

fn capsule_normal(q: vec3f, a: vec3f, b: vec3f, r: f32) -> vec4f {
    let ba = b - a;
    let bb = dot(ba, ba);
    var f = 0.0;
    if bb > 1e-12 {
        f = clamp(dot(q - a, ba) / bb, 0.0, 1.0);
    }
    return vec4f((q - (a + ba * f)) / r, f);
}

// Where the vertical line through (x, y) enters a primitive and leaves it.
fn column(q: Prim, x: f32, y: f32) -> Span {
    var s = Span(false, 0.0, 0.0, vec3f(0.0), 0.0);
    let kind = q.v0.w;
    if kind < 0.5 {
        let c = q.v0.xyz;
        let r = q.v1.w;
        let dx = x - c.x;
        let dy = y - c.y;
        let d2 = dx * dx + dy * dy;
        if d2 > r * r {
            return s;
        }
        let dz = sqrt(r * r - d2);
        return Span(true, c.z + dz, c.z - dz, vec3f(dx / r, dy / r, dz / r), 0.0);
    }
    if kind < 1.5 {
        let a = q.v0.xyz;
        let b = q.v1.xyz;
        let r = q.v1.w;
        let top = max(a.z, b.z) + r + 1.0;
        let bottom = min(a.z, b.z) - r - 1.0;
        let t = capsule_hit(vec3f(x, y, top), vec3f(0.0, 0.0, -1.0), a, b, r);
        if t.y < 0.5 {
            return s;
        }
        let u = capsule_hit(vec3f(x, y, bottom), vec3f(0.0, 0.0, 1.0), a, b, r);
        var up = top - bottom - t.x;
        if u.y > 0.5 {
            up = u.x;
        }
        let front = top - t.x;
        let nf = capsule_normal(vec3f(x, y, front), a, b, r);
        return Span(true, front, bottom + up, nf.xyz, nf.w);
    }
    let v0 = q.v0.xyz;
    let v1 = q.v1.xyz;
    let v2 = q.v2.xyz;
    let d = (v1.y - v2.y) * (v0.x - v2.x) + (v2.x - v1.x) * (v0.y - v2.y);
    if abs(d) < 1e-12 {
        return s;
    }
    let l0 = ((v1.y - v2.y) * (x - v2.x) + (v2.x - v1.x) * (y - v2.y)) / d;
    let l1 = ((v2.y - v0.y) * (x - v2.x) + (v0.x - v2.x) * (y - v2.y)) / d;
    let l2 = 1.0 - l0 - l1;
    if l0 < 0.0 || l1 < 0.0 || l2 < 0.0 {
        return s;
    }
    let z = l0 * v0.z + l1 * v1.z + l2 * v2.z;
    var n = unit(q.n0.xyz * l0 + q.n1.xyz * l1 + q.n2.xyz * l2);
    // Seen from both sides: lit on whichever faces us.
    if n.z < 0.0 {
        n = -n;
    }
    return Span(true, z, z, n, 0.0);
}

struct Surface {
    hit: bool,
    z: f32,
    n: vec3f,
    base: vec3f,
}

// The nearest opaque surface down the line through (x, y), cut by the slab.
fn first_surface(x: f32, y: f32) -> Surface {
    var best = Surface(false, 0.0, vec3f(0.0), vec3f(0.0));
    var stack: array<u32, STACK>;
    var top = 1u;
    stack[0] = 0u;
    var floor = p.cut.y;
    let near = p.cut.x;
    let far = p.cut.y;
    while top > 0u {
        top -= 1u;
        let node = nodes[stack[top]];
        if x < node.lo.x || x > node.hi.x || y < node.lo.y || y > node.hi.y || node.hi.z < floor {
            continue;
        }
        if node.count > 0u {
            for (var i = node.first; i < node.first + node.count; i++) {
                let q = prims[i];
                let s = column(q, x, y);
                if !s.hit {
                    continue;
                }
                let base = mix(q.ca.xyz, q.cb.xyz, s.f);
                var found = false;
                var z = 0.0;
                var n = vec3f(0.0);
                var b = base;
                if s.front > near {
                    if s.back <= near && q.v0.w < 1.5 {
                        found = true;
                        z = near;
                        n = CAP;
                        b = base * CAP_SHADE;
                    }
                } else if s.front >= far {
                    found = true;
                    z = s.front;
                    n = s.n;
                }
                if found && (!best.hit || z > best.z) {
                    best = Surface(true, z, n, b);
                    floor = z;
                }
            }
        } else if top + 2u <= STACK {
            // The higher child last, so it is searched first and raises the
            // floor under the other.
            let a = node.first;
            let c = node.first + 1u;
            if nodes[a].hi.z > nodes[c].hi.z {
                stack[top] = c;
                stack[top + 1u] = a;
            } else {
                stack[top] = a;
                stack[top + 1u] = c;
            }
            top += 2u;
        }
    }
    return best;
}

// Whether a ray from o along unit d meets anything within reach whose point
// lies between the slab's depths.
fn blocked(o: vec3f, d: vec3f, reach: f32) -> bool {
    let inv = 1.0 / d;
    var stack: array<u32, STACK>;
    var top = 1u;
    stack[0] = 0u;
    let near = p.cut.x;
    let far = p.cut.y;
    while top > 0u {
        top -= 1u;
        let node = nodes[stack[top]];
        let a = (node.lo - o) * inv;
        let b = (node.hi - o) * inv;
        let lo = min(a, b);
        let hi = max(a, b);
        let t0 = max(max(max(lo.x, lo.y), lo.z), 0.0);
        let t1 = min(min(min(hi.x, hi.y), hi.z), reach);
        if t0 > t1 {
            continue;
        }
        if node.count > 0u {
            for (var i = node.first; i < node.first + node.count; i++) {
                let h = ray_hit(prims[i], o, d);
                if h.y > 0.5 {
                    let z = o.z + d.z * h.x;
                    if h.x > 0.0 && h.x < reach && z <= near && z >= far {
                        return true;
                    }
                }
            }
        } else if top + 2u <= STACK {
            stack[top] = node.first;
            stack[top + 1u] = node.first + 1u;
            top += 2u;
        }
    }
    return false;
}

// -- sampling and light ---------------------------------------------------

fn hash(x: u32, y: u32, salt: u32) -> u32 {
    var h = (x * 0x8da6b343u) ^ (y * 0xd8163841u) ^ (salt * 0xcb1ab31fu);
    h ^= h >> 16u;
    h *= 0x7feb352du;
    h ^= h >> 15u;
    h *= 0x846ca68bu;
    return h ^ (h >> 16u);
}

fn r2(i: u32, shift: u32) -> vec2f {
    let s = f32(shift & 0xffffu) / 65536.0;
    let t = f32(shift >> 16u) / 65536.0;
    return vec2f(fract(0.5 + f32(i) * 0.7548777 + s), fract(0.5 + f32(i) * 0.5698403 + t));
}

fn shade_lit(n: vec3f, base: vec3f, style: f32, lit: f32, open: f32) -> vec3f {
    let diff = max(dot(n, p.light.xyz), 0.0) * lit;
    let nh = max(dot(n, HALF), 0.0);
    let n2 = nh * nh;
    let n8 = n2 * n2 * n2 * n2;
    var k: f32;
    var spec: f32;
    if style < 0.5 {
        k = 0.34 + 0.70 * diff;
        spec = n8 * n8 * n8 * 0.28;
    } else if style < 1.5 {
        let banded = floor(diff * 4.0) / 3.0;
        let mixed = diff * 0.35 + min(banded, 1.0) * 0.65;
        k = 0.62 + 0.42 * mixed;
        spec = n8 * n8 * 0.10;
    } else if style < 2.5 {
        k = 0.76 + 0.26 * diff;
        spec = 0.0;
    } else {
        let n32 = n8 * n8 * n8 * n8;
        k = 0.26 + 0.80 * diff;
        spec = n32 * 0.5 + n8 * 0.04;
    }
    k *= open;
    spec *= lit;
    return base * k + vec3f(spec);
}

@compute @workgroup_size(8, 8)
fn trace(@builtin(global_invocation_id) id: vec3u) {
    let width = p.size.x;
    let rows = p.size.z;
    if id.x >= width || id.y >= rows {
        return;
    }
    let px = id.x;
    let py = p.size.y + id.y;
    let slot = id.y * width + id.x;
    let x = (f32(px) + 0.5 - p.view.x) * p.view.z;
    let y = -(f32(py) + 0.5 - p.view.y) * p.view.z;
    let s = first_surface(x, y);
    if !s.hit {
        out[slot] = vec4f(0.0, 0.0, 0.0, EMPTY);
        return;
    }
    let o = vec3f(x, y, s.z) + s.n * p.view.w;
    let ss = p.size.w;
    let ox = px / ss;
    let oy = py / ss;
    let part = (py % ss) * ss + px % ss;
    let light = p.light.xyz;

    var lit = 1.0;
    if dot(s.n, light) > 0.0 {
        let seed = hash(ox, oy, 0u);
        let count = p.rays.x;
        var through = 0u;
        for (var i = 0u; i < count; i++) {
            let uv = r2(part * count + i, seed);
            let rr = sqrt(uv.x) * LAMP;
            let phi = uv.y * TAU;
            let d = unit(light + p.lt.xyz * (rr * cos(phi)) + p.lb.xyz * (rr * sin(phi)));
            if !blocked(o, d, 3.0e38) {
                through += 1u;
            }
        }
        lit = 1.0 - p.cut.w * (1.0 - f32(through) / f32(count));
    }

    var up = vec3f(0.0, 0.0, 1.0);
    if abs(s.n.z) >= 0.9 {
        up = vec3f(1.0, 0.0, 0.0);
    }
    let t = unit(cross(up, s.n));
    let b = cross(s.n, t);
    let seed = hash(ox, oy, 1u);
    let count = p.rays.y;
    var shut = 0u;
    for (var i = 0u; i < count; i++) {
        let uv = r2(part * count + i, seed);
        let rr = sqrt(uv.x);
        let phi = uv.y * TAU;
        let rise = sqrt(max(1.0 - uv.x, 0.0));
        let d = t * (rr * cos(phi)) + b * (rr * sin(phi)) + s.n * rise;
        if blocked(o, d, p.cut.z) {
            shut += 1u;
        }
    }
    let open = 1.0 - p.look.x * f32(shut) / f32(count);
    out[slot] = vec4f(shade_lit(s.n, s.base, p.look.y, lit, open), s.z);
}
