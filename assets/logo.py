# Draws assets/logo.svg: `python3 assets/logo.py assets/logo.svg`.
#
# Caffeine, the first of the built-in samples, as a skeletal formula in one
# weight of ink on paper. Coordinates come from samples/caffeine.sdf.
import math, os, sys

INK, PAPER = '#16171d', '#fdfcf9'
FONT = "'Helvetica Neue', Helvetica, Arial, sans-serif"
BOX = 700.0                       # the drawing's longest side, in px
LINE = 0.10                       # stroke width, as a fraction of a bond
LABEL = 0.50                      # cap height of N and O, as a fraction of a bond

sdf = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'samples', 'caffeine.sdf')
lines = open(sdf).read().splitlines()
na, nb = int(lines[3][:3]), int(lines[3][3:6])
atoms = []
for l in lines[4:4+na]:
    f = l.split()
    atoms.append((float(f[0]), float(f[1]), f[3]))
bonds = []
for l in lines[4+na:4+na+nb]:
    a, b, order = int(l[:3])-1, int(l[3:6])-1, int(l[6:9])
    if atoms[a][2] != 'H' and atoms[b][2] != 'H':
        bonds.append((a, b, order))
heavy = sorted({i for b in bonds for i in b[:2]})

# rings, for which side a ring double bond's second line goes
adj = {i: set() for i in heavy}
for a, b, _ in bonds:
    adj[a].add(b); adj[b].add(a)
def ring_center(a, b):
    # the smallest ring through the bond a-b, by breadth-first search
    from collections import deque
    q = deque([(b, [b])])
    best = None
    while q:
        node, path = q.popleft()
        if len(path) > 7: break
        for m in adj[node]:
            if m == a and len(path) > 2:
                ring = [a] + path
                if best is None or len(ring) < len(best): best = ring
            elif m not in path and m != a:
                q.append((m, path + [m]))
    if not best: return None
    return (sum(P[i][0] for i in best)/len(best), sum(P[i][1] for i in best)/len(best), len(best))

# Turn it so the bond the two rings share stands upright, the way caffeine
# is usually drawn, with the six-membered ring on the left.
fused = next((a, b) for a, b, _ in bonds
             if len(adj[a] & adj[b]) == 0 and len(adj[a]) == 3 and len(adj[b]) == 3
             and all(atoms[i][2] == 'C' for i in (a, b)))
(ax, ay), (bx, by) = atoms[fused[0]][:2], atoms[fused[1]][:2]
r = -math.pi/2 - math.atan2(by-ay, bx-ax)
P = {i: (atoms[i][0]*math.cos(r) - atoms[i][1]*math.sin(r),
         -(atoms[i][0]*math.sin(r) + atoms[i][1]*math.cos(r))) for i in heavy}
bond = sum(math.dist(P[a], P[b]) for a, b, _ in bonds)/len(bonds)
xs = [p[0] for p in P.values()]; ys = [p[1] for p in P.values()]
pad = 0.45*bond
x0, x1, y0, y1 = min(xs)-pad, max(xs)+pad, min(ys)-pad, max(ys)+pad
sc = BOX/max(x1-x0, y1-y0)
P = {i: ((p[0]-(x0+x1)/2)*sc + 512, (p[1]-(y0+y1)/2)*sc + 512) for i, p in P.items()}
B = bond*sc
W = LINE*B
cap = LABEL*B
labelled = {i for i in heavy if atoms[i][2] != 'C'}
clear = 0.78*cap                   # bonds stop this far from a letter

def seg(p, q, a=None, b=None):
    (x1, y1), (x2, y2) = p, q
    d = math.dist(p, q); ux, uy = (x2-x1)/d, (y2-y1)/d
    if a in labelled: x1, y1 = x1+ux*clear, y1+uy*clear
    if b in labelled: x2, y2 = x2-ux*clear, y2-uy*clear
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}"/>'

o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1024 1024" width="1024" height="1024">',
 '<defs><filter id="sh" x="-20%" y="-20%" width="140%" height="140%">'
 '<feDropShadow dx="0" dy="10" stdDeviation="14" flood-color="#1f1d17" flood-opacity="0.18"/></filter></defs>',
 f'<rect x="100" y="100" width="824" height="824" rx="185" fill="{PAPER}" filter="url(#sh)"/>',
 f'<g stroke="{INK}" stroke-width="{W:.1f}" stroke-linecap="round">']
gap = 0.20*B
for a, b, order in bonds:
    p, q = P[a], P[b]
    if order == 1:
        o.append(seg(p, q, a, b)); continue
    d = math.dist(p, q); nx, ny = -(q[1]-p[1])/d, (q[0]-p[0])/d
    rc = ring_center(a, b)
    if rc:
        # ring bond: the main line, and a shorter one inside the ring
        o.append(seg(p, q, a, b))
        mx, my = (p[0]+q[0])/2, (p[1]+q[1])/2
        if (rc[0]-mx)*nx + (rc[1]-my)*ny < 0: nx, ny = -nx, -ny
        k = 0.16
        p2 = (p[0]+(q[0]-p[0])*k + nx*gap, p[1]+(q[1]-p[1])*k + ny*gap)
        q2 = (q[0]+(p[0]-q[0])*k + nx*gap, q[1]+(p[1]-q[1])*k + ny*gap)
        o.append(seg(p2, q2, a, b))
    else:
        # exocyclic, as two lines either side of the axis
        h = gap/2
        o.append(seg((p[0]+nx*h, p[1]+ny*h), (q[0]+nx*h, q[1]+ny*h), a, b))
        o.append(seg((p[0]-nx*h, p[1]-ny*h), (q[0]-nx*h, q[1]-ny*h), a, b))
o.append('</g>')
for i in sorted(labelled):
    x, y = P[i]
    o.append(f'<text x="{x:.1f}" y="{y + cap*0.5:.1f}" font-family="{FONT}" font-weight="500" '
             f'font-size="{cap*1.4:.1f}" text-anchor="middle" fill="{INK}">{atoms[i][2]}</text>')
o.append('</svg>')
open(sys.argv[1], 'w').write('\n'.join(o))
