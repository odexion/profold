# Draws assets/logo.svg: `python3 assets/logo.py assets/logo.svg`.
#
# One ink line on paper: a strand that folds into a helix. Where the chain
# passes behind itself the back stroke is broken, as in a knot diagram.
import math, sys

INK, PAPER = '#16171d', '#fdfcf9'

turns, R, H = 2.5, 1.0, 3.4
def helix(s):
    th = s*turns*2*math.pi - math.pi/2
    return (R*math.cos(th), s*H, R*math.sin(th))

# incoming strand, a cubic Bezier that meets the helix along its tangent
th0 = -math.pi/2
tan = (-R*math.sin(th0)*turns*2*math.pi, H, R*math.cos(th0)*turns*2*math.pi)
L = math.sqrt(sum(v*v for v in tan)); tan = tuple(v/L for v in tan)
P0, P1, P3 = (-0.9, -2.0, 0.5), (-0.9, -1.3, 0.6), helix(0.0)
P2 = tuple(P3[k] - tan[k]*1.4 for k in range(3))
raw = []
for i in range(90):
    t = i/90
    raw.append(tuple((1-t)**3*P0[k] + 3*(1-t)**2*t*P1[k] + 3*(1-t)*t*t*P2[k] + t**3*P3[k] for k in range(3)))
for i in range(601):
    raw.append(helix(i/600))

def view(p):
    x, y, z = p
    a = math.radians(38)          # pitch: see into the coils
    y, z = y*math.cos(a) - z*math.sin(a), y*math.sin(a) + z*math.cos(a)
    b = math.radians(-40)         # roll: lean the axis right
    x, y = x*math.cos(b) - y*math.sin(b), x*math.sin(b) + y*math.cos(b)
    return (x, -y, z)
P = [view(p) for p in raw]

LINE = 0.17                       # stroke width, model units
GAP = 0.12                        # paper either side of a stroke that crosses over
xs = [p[0] for p in P]; ys = [p[1] for p in P]
x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
sc = 680.0/max(x1-x0, y1-y0)
cx, cy = 512 - (x0+x1)/2*sc, 512 - (y0+y1)/2*sc
pts = [(p[0]*sc+cx, p[1]*sc+cy, p[2]) for p in P]
n = len(pts)
W, G = LINE*sc, GAP*sc

# A point is hidden where another part of the chain, nearer the eye, passes
# over it on the page. What is left is drawn as separate strokes.
reach = W + G
hidden = [any(abs(j-k) > 25 and pts[j][2] > pts[k][2]
              and math.hypot(pts[j][0]-pts[k][0], pts[j][1]-pts[k][1]) < reach
              for j in range(n)) for k in range(n)]
strokes, cur = [], []
for k in range(n):
    if hidden[k]:
        if len(cur) > 1: strokes.append(cur)
        cur = []
    else:
        cur.append(k)
if len(cur) > 1: strokes.append(cur)

o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1024 1024" width="1024" height="1024">',
 '<defs><filter id="sh" x="-20%" y="-20%" width="140%" height="140%">'
 '<feDropShadow dx="0" dy="10" stdDeviation="14" flood-color="#1f1d17" flood-opacity="0.18"/></filter></defs>',
 f'<rect x="100" y="100" width="824" height="824" rx="185" fill="{PAPER}" filter="url(#sh)"/>',
 f'<g fill="none" stroke="{INK}" stroke-width="{W:.1f}" stroke-linejoin="round" stroke-linecap="round">']
for st in strokes:
    o.append('<polyline points="' + ' '.join(f'{pts[k][0]:.1f},{pts[k][1]:.1f}' for k in st) + '"/>')
o.append('</g></svg>')
open(sys.argv[1], 'w').write('\n'.join(o))
