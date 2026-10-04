# Draws assets/logo.svg: `python3 assets/logo.py assets/logo.svg`.
import math, colorsys, sys

# The Spectrum ramp from src/render.rs: hue 250 -> 0, s .62, l .52
def spectrum(t):
    return colorsys.hls_to_rgb(((1-t)*250)/360, 0.52, 0.62)
def hexc(c, k=1.0):
    return '#%02x%02x%02x' % tuple(max(0,min(255,round(v*k*255))) for v in c)

# Backbone in 3D: an unfolded strand that coils into an alpha helix.
turns, R, H = 2.2, 1.0, 2.85
def helix(s):
    th = s*turns*2*math.pi - math.pi/2
    return (R*math.cos(th), s*H, R*math.sin(th))
raw = []
start = helix(0.0)
strand = [(-3.1, -0.25, 0.4), (-2.0, -0.55, 0.3), (-0.9, -0.55, 0.0), start]
def bez(p, t):  # Catmull-Rom through strand points
    return p
# smooth incoming strand with a cubic Bezier into the helix tangent
th0 = -math.pi/2
tan = (-R*math.sin(th0)*turns*2*math.pi, H, R*math.cos(th0)*turns*2*math.pi)
L = math.sqrt(sum(v*v for v in tan)); tan = tuple(v/L for v in tan)
P0, P3 = (-3.0, -0.6, 0.6), start
P1 = (-2.0, -1.2, 0.5)
P2 = tuple(P3[k] - tan[k]*1.4 for k in range(3))
for i in range(90):
    t = i/90
    raw.append(tuple((1-t)**3*P0[k] + 3*(1-t)**2*t*P1[k] + 3*(1-t)*t*t*P2[k] + t**3*P3[k] for k in range(3)))
for i in range(601):
    raw.append(helix(i/600))

# View: tilt the helix axis toward upper-right and look slightly from above.
def rot(p):
    x,y,z = p
    a = math.radians(24)          # pitch: see into the coils
    y, z = y*math.cos(a) - z*math.sin(a), y*math.sin(a) + z*math.cos(a)
    b = math.radians(-30)         # roll: lean the axis right
    x, y = x*math.cos(b) - y*math.sin(b), x*math.sin(b) + y*math.cos(b)
    return (x, -y, z)             # screen y grows downward
P = [rot(p) for p in raw]
xs = [p[0] for p in P]; ys = [p[1] for p in P]
W_units = 0.5                    # tube width in model units
pad = W_units/2 + 0.05
x0, x1, y0, y1 = min(xs)-pad, max(xs)+pad, min(ys)-pad, max(ys)+pad
box = 680.0
sc = box/max(x1-x0, y1-y0)
cx, cy = 512 - (x0+x1)/2*sc, 520 - (y0+y1)/2*sc
pts = [(p[0]*sc+cx, p[1]*sc+cy, p[2]) for p in P]
W = W_units*sc
h = W*0.11

n = len(pts)
zmin = min(p[2] for p in pts); zmax = max(p[2] for p in pts)
# Runs of the backbone behind and in front of the helix axis. Back runs are
# drawn first, so front turns cross over them; a run's ends are butt so the
# joins between runs leave no seam.
runs, s0 = [], 0
for k in range(1, n):
    if (pts[k][2] >= 0) != (pts[s0][2] >= 0) or k == n-1:
        runs.append((s0, k)); s0 = k
# fold slivers into their neighbour
while any(e - s < 12 for s, e in runs) and len(runs) > 1:
    i = next(i for i, (s, e) in enumerate(runs) if e - s < 12)
    j = i - 1 if i > 0 else i + 1
    a, b = sorted((i, j))
    runs[a:b+1] = [(runs[a][0], runs[b][1])]
runs.sort(key=lambda r: (pts[r[0]][2] + pts[r[1]][2] + pts[(r[0]+r[1])//2][2]) >= 0)

def pl(a, b, dx=0, dy=0):
    return ' '.join(f'{pts[k][0]+dx:.1f},{pts[k][1]+dy:.1f}' for k in range(a, b+1))
def col(k):
    depth = (pts[k][2]-zmin)/(zmax-zmin)
    return hexc(spectrum(k/(n-1)), 0.82 + 0.18*depth)

INK = '#16171d'
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1024 1024" width="1024" height="1024">',
 '<defs>',
 '<linearGradient id="bg" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#fffdf8"/><stop offset="1" stop-color="#f3eee3"/></linearGradient>',
 '<filter id="sh" x="-20%" y="-20%" width="140%" height="140%"><feDropShadow dx="0" dy="12" stdDeviation="16" flood-color="#1f1d17" flood-opacity="0.26"/></filter>',
 '</defs>',
 '<rect x="100" y="100" width="824" height="824" rx="185" fill="url(#bg)" filter="url(#sh)"/>',
 '<rect x="100.5" y="100.5" width="823" height="823" rx="184.5" fill="none" stroke="#e4ddcf" stroke-width="1"/>',
 '<g fill="none" stroke-linejoin="round" stroke-linecap="butt">']
SUB = 3
def tube(s, e, ext, caps):
    # outline, then the fill in short solid pieces that follow the ramp
    o.append(f'<polyline points="{pl(s,e)}" stroke="{INK}" stroke-width="{W+2*h:.1f}" stroke-linecap="{caps}"/>')
    fs, fe = max(0, s-ext), min(n-1, e+ext)
    for k in (fs, fe):
        if k in (0, n-1):
            o.append(f'<circle cx="{pts[k][0]:.1f}" cy="{pts[k][1]:.1f}" r="{W/2:.1f}" fill="{col(k)}"/>')
    for a in range(fs, fe, SUB):
        b = min(fe, a+SUB+1)
        o.append(f'<polyline points="{pl(a,b)}" stroke="{col((a+b)//2)}" stroke-width="{W:.1f}"/>')
    o.append(f'<polyline points="{pl(fs+3,fe-3,-W*0.12,-W*0.15)}" stroke="#fff" stroke-opacity="0.28" stroke-width="{W*0.2:.1f}" stroke-linecap="round"/>')
# The whole chain once, then the turns in front of the axis over it. Every
# join lands on the tube already drawn, so none shows.
tube(0, n-1, 0, 'round')
for s, e in runs:
    if pts[(s+e)//2][2] > 0:
        tube(s, e, 16, 'butt')
o.append('</g></svg>')
open(sys.argv[1], 'w').write('\n'.join(o))
