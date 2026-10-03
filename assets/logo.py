# Draws assets/logo.svg: `python3 assets/logo.py assets/logo.svg`.
#
# A DNA double helix in ball and stick, as the viewer draws it: one colour,
# shaded spheres, a dark hairline. Two backbones of beads, one rung per base
# pair; what lies further back fades a little toward the paper.
import math, sys

PAPER = '#fdfdfc'
LIGHT, MID, BASE, LINE = (0x4b, 0x37, 0x69), (0x3a, 0x28, 0x56), (0x2b, 0x1b, 0x40), (0x13, 0x13, 0x1b)
PAIRS = 12                        # base pairs drawn
TWIST = 36.0                      # degrees per pair
RISE = 0.48                       # axis rise per pair, in helix radii
GROOVE = 180.0                    # angle between the two strands' beads
BEAD = 0.21                       # sphere radius, in helix radii
STICK = 0.065
EDGE = 0.016
FADE = 0.28                       # how far the back fades toward the paper
TILT = -34.0                      # lean of the axis on the tile
PITCH = 0.0                       # look down the axis a little, into the rungs
BOX = 690.0

pts, sticks = [], []
for i in range(PAIRS):
    th = math.radians(i*TWIST + TWIST/2)
    y = (i - (PAIRS-1)/2)*RISE
    for s, off in ((0, 0.0), (1, math.radians(GROOVE))):
        pts.append((math.cos(th+off), y, math.sin(th+off)))
    a, b = 2*i, 2*i+1
    sticks.append((a, b))                      # the base pair
    if i:
        sticks += [(a-2, a), (b-2, b)]         # the backbones

def view(p):
    x, y, z = p
    q = math.radians(PITCH)
    y, z = y*math.cos(q) - z*math.sin(q), y*math.sin(q) + z*math.cos(q)
    t = math.radians(TILT)
    x, y = x*math.cos(t) - y*math.sin(t), x*math.sin(t) + y*math.cos(t)
    return (x, -y, z)
P = [view(p) for p in pts]
x0 = min(p[0] for p in P)-BEAD; x1 = max(p[0] for p in P)+BEAD
y0 = min(p[1] for p in P)-BEAD; y1 = max(p[1] for p in P)+BEAD
sc = BOX/max(x1-x0, y1-y0)
xy = [((p[0]-(x0+x1)/2)*sc + 512, (p[1]-(y0+y1)/2)*sc + 512) for p in P]
zs = [p[2] for p in P]

PAP = tuple(int(PAPER[k:k+2], 16) for k in (1, 3, 5))
def tint(c, z):
    f = FADE*(1 - (z + 1)/2)                  # 0 at the front, FADE at the back
    return '#%02x%02x%02x' % tuple(round(c[k] + (PAP[k]-c[k])*f) for k in range(3))

o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1024 1024" width="1024" height="1024">',
 '<defs><filter id="sh" x="-20%" y="-20%" width="140%" height="140%">'
 '<feDropShadow dx="0" dy="10" stdDeviation="14" flood-color="#1f1d17" flood-opacity="0.18"/></filter>']
for i, z in enumerate(zs):
    o.append(f'<radialGradient id="b{i}" cx="0.42" cy="0.38" r="0.62">'
             f'<stop offset="0" stop-color="{tint(LIGHT, z)}"/><stop offset="0.55" stop-color="{tint(MID, z)}"/>'
             f'<stop offset="1" stop-color="{tint(BASE, z)}"/></radialGradient>')
o += ['</defs>', f'<rect x="100" y="100" width="824" height="824" rx="185" fill="{PAPER}" filter="url(#sh)"/>']

# Back to front. A stick sits just behind the nearer of its two beads.
items = [(z, 1, i) for i, z in enumerate(zs)]
items += [(min(zs[a], zs[b]) - 1e-3, 0, k) for k, (a, b) in enumerate(sticks)]
items.sort()
e, r, w = EDGE*sc, BEAD*sc, STICK*2*sc
for z, kind, k in items:
    if kind:
        x, y = xy[k]
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="url(#b{k})" stroke="{tint(LINE, z)}" stroke-width="{e:.1f}"/>')
    else:
        a, b = sticks[k]
        (ax, ay), (bx, by) = xy[a], xy[b]
        seg = f'x1="{ax:.1f}" y1="{ay:.1f}" x2="{bx:.1f}" y2="{by:.1f}"'
        o.append(f'<line {seg} stroke="{tint(LINE, z)}" stroke-width="{w+2*e:.1f}"/>')
        o.append(f'<line {seg} stroke="{tint(BASE, z)}" stroke-width="{w:.1f}"/>')
        d = math.hypot(bx-ax, by-ay) or 1.0; nx, ny = -(by-ay)/d, (bx-ax)/d
        if nx + ny > 0: nx, ny = -nx, -ny
        o.append(f'<line x1="{ax+nx*w*0.14:.1f}" y1="{ay+ny*w*0.14:.1f}" x2="{bx+nx*w*0.14:.1f}" y2="{by+ny*w*0.14:.1f}" '
                 f'stroke="{tint(MID, z)}" stroke-width="{w*0.5:.1f}"/>')
o.append('</svg>')
open(sys.argv[1], 'w').write('\n'.join(o))
