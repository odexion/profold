# Draws assets/logo.svg: `python3 assets/logo.py assets/logo.svg`.
#
# Caffeine, the first of the built-in samples, in ball and stick as the
# viewer draws it: one colour, shaded spheres, a dark hairline. Coordinates
# come from samples/caffeine.sdf and are seen down z, as the app opens it.
import math, os, sys

PAPER = '#fdfdfc'
LIGHT, BASE, LINE = '#4b3769', '#2b1b40', '#13131b'
BOX = 700.0                       # the drawing's longest side, in px
RADIUS = {'H': 0.32}              # sphere radius in Å; heavy atoms get HEAVY
HEAVY = 0.45
STICK = 0.13                      # stick radius in Å
EDGE = 0.022                      # outline width in Å

sdf = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'samples', 'caffeine.sdf')
lines = open(sdf).read().splitlines()
na, nb = int(lines[3][:3]), int(lines[3][3:6])
atoms = [(float(l.split()[0]), float(l.split()[1]), float(l.split()[2]), l.split()[3])
         for l in lines[4:4+na]]
bonds = [(int(l[:3])-1, int(l[3:6])-1) for l in lines[4+na:4+na+nb]]
rad = [RADIUS.get(a[3], HEAVY) for a in atoms]

xs = [a[0] for a in atoms]; ys = [a[1] for a in atoms]
x0 = min(a[0]-r for a, r in zip(atoms, rad)); x1 = max(a[0]+r for a, r in zip(atoms, rad))
y0 = min(a[1]-r for a, r in zip(atoms, rad)); y1 = max(a[1]+r for a, r in zip(atoms, rad))
sc = BOX/max(x1-x0, y1-y0)
def at(i):
    a = atoms[i]
    return ((a[0]-(x0+x1)/2)*sc + 512, -(a[1]-(y0+y1)/2)*sc + 512)

o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1024 1024" width="1024" height="1024">',
 '<defs>',
 '<filter id="sh" x="-20%" y="-20%" width="140%" height="140%">'
 '<feDropShadow dx="0" dy="10" stdDeviation="14" flood-color="#1f1d17" flood-opacity="0.18"/></filter>',
 # a sphere lit from the upper left
 f'<radialGradient id="ball" cx="0.42" cy="0.38" r="0.62"><stop offset="0" stop-color="{LIGHT}"/>'
 f'<stop offset="0.55" stop-color="#3a2856"/><stop offset="1" stop-color="{BASE}"/></radialGradient>',
 '</defs>',
 f'<rect x="100" y="100" width="824" height="824" rx="185" fill="{PAPER}" filter="url(#sh)"/>']

# Back to front. A stick sits just behind the nearer of its two atoms.
items = [(atoms[i][2], 1, i) for i in range(na)]
items += [(min(atoms[a][2], atoms[b][2]) - 1e-3, 0, k) for k, (a, b) in enumerate(bonds)]
items.sort()
e = EDGE*sc
for z, kind, k in items:
    if kind == 1:
        x, y = at(k); r = rad[k]*sc
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="url(#ball)" stroke="{LINE}" stroke-width="{e:.1f}"/>')
    else:
        a, b = bonds[k]
        (ax, ay), (bx, by) = at(a), at(b)
        w = STICK*2*sc
        seg = f'x1="{ax:.1f}" y1="{ay:.1f}" x2="{bx:.1f}" y2="{by:.1f}"'
        o.append(f'<line {seg} stroke="{LINE}" stroke-width="{w+2*e:.1f}"/>')
        o.append(f'<line {seg} stroke="{BASE}" stroke-width="{w:.1f}"/>')
        # the round of the stick: a lighter band toward the light
        d = math.hypot(bx-ax, by-ay); nx, ny = -(by-ay)/d, (bx-ax)/d
        if nx + ny > 0: nx, ny = -nx, -ny
        off = w*0.14
        o.append(f'<line x1="{ax+nx*off:.1f}" y1="{ay+ny*off:.1f}" x2="{bx+nx*off:.1f}" y2="{by+ny*off:.1f}" '
                 f'stroke="#3d2a5a" stroke-width="{w*0.5:.1f}"/>')
o.append('</svg>')
open(sys.argv[1], 'w').write('\n'.join(o))
