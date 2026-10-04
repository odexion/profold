# Draws assets/logo.svg: `python3 assets/logo.py assets/logo.svg`.
#
# A strand that folds into a helix, coloured along its length from the
# N-terminus to the C. No outline and no gloss: where the chain passes behind
# itself it simply breaks, and the paper shows through. The strand swells
# where it comes toward the eye and thins where it turns away.
import math, sys

PAPER = '#ffffff'
# the ramp along the chain, start to end
RAMP = [(0x43, 0x38, 0xca), (0x7c, 0x3a, 0xed), (0xc0, 0x26, 0xd3), (0xf4, 0x3f, 0x5e), (0xfb, 0x92, 0x3c)]
TURNS, H = 2.18, 4.1              # helix turns and length, in radii
LINE = 0.36                       # strand width, in radii
THIN = 0.45                       # width at the far side, as a fraction of the near
TAPER = 0.06                      # the ends narrow over this fraction of the chain
GAP = 0.08                        # paper either side of a strand that passes over
BOX = 690.0

def ramp(t):
    t = max(0.0, min(1.0, t))*(len(RAMP)-1)
    i = min(int(t), len(RAMP)-2); f = t - i
    a, b = RAMP[i], RAMP[i+1]
    return '#%02x%02x%02x' % tuple(round(a[k] + (b[k]-a[k])*f) for k in range(3))

def helix(s):
    th = s*TURNS*2*math.pi - math.pi/2
    return (math.cos(th), s*H, math.sin(th))
# the incoming strand, a cubic Bezier that meets the helix along its tangent
th0 = -math.pi/2
tan = (-math.sin(th0)*TURNS*2*math.pi, H, math.cos(th0)*TURNS*2*math.pi)
L = math.sqrt(sum(v*v for v in tan)); tan = tuple(v/L for v in tan)
P0, P1, P3 = (-0.9, -2.0, 0.5), (-0.9, -1.3, 0.6), helix(0.0)
P2 = tuple(P3[k] - tan[k]*1.4 for k in range(3))
raw = []
for i in range(120):
    t = i/120
    raw.append(tuple((1-t)**3*P0[k] + 3*(1-t)**2*t*P1[k] + 3*(1-t)*t*t*P2[k] + t**3*P3[k] for k in range(3)))
for i in range(801):
    raw.append(helix(i/800))

def view(p):
    x, y, z = p
    a = math.radians(38)          # pitch: see into the coils
    y, z = y*math.cos(a) - z*math.sin(a), y*math.sin(a) + z*math.cos(a)
    b = math.radians(-40)         # roll: lay the axis on the diagonal
    x, y = x*math.cos(b) - y*math.sin(b), x*math.sin(b) + y*math.cos(b)
    return (x, -y, z)
P = [view(p) for p in raw]
x0 = min(p[0] for p in P) - LINE/2; x1 = max(p[0] for p in P) + LINE/2
y0 = min(p[1] for p in P) - LINE/2; y1 = max(p[1] for p in P) + LINE/2
sc = BOX/max(x1-x0, y1-y0)
pts = [((p[0]-(x0+x1)/2)*sc + 512, (p[1]-(y0+y1)/2)*sc + 512, p[2]) for p in P]
n = len(pts)
W, G = LINE*sc, GAP*sc
zmin = min(p[2] for p in pts); zmax = max(p[2] for p in pts)
def width(k):
    t = (pts[k][2] - zmin)/(zmax - zmin)          # 0 far, 1 near
    e = min(1.0, min(k, n-1-k)/((n-1)*TAPER))     # 0 at an end
    e = e*e*(3 - 2*e)                             # eased, so the taper has no corner
    return W*(THIN + (1-THIN)*t)*(0.6 + 0.4*e)
wd = [width(k) for k in range(n)]

# A point is hidden where another part of the chain, nearer the eye, passes
# over it on the page.
arc = [0.0]
for k in range(1, n):
    arc.append(arc[-1] + math.hypot(pts[k][0]-pts[k-1][0], pts[k][1]-pts[k-1][1]))
hidden = [any(abs(arc[j]-arc[k]) > 3*W and pts[j][2] > pts[k][2]
              and math.hypot(pts[j][0]-pts[k][0], pts[j][1]-pts[k][1]) < (wd[j]+wd[k])/2 + G
              for j in range(0, n, 2)) for k in range(n)]

o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1024 1024" width="1024" height="1024">',
 '<defs><filter id="sh" x="-20%" y="-20%" width="140%" height="140%">'
 '<feDropShadow dx="0" dy="8" stdDeviation="16" flood-color="#0f0a1e" flood-opacity="0.14"/></filter></defs>',
 f'<rect x="100" y="100" width="824" height="824" rx="185" fill="{PAPER}" filter="url(#sh)"/>',
 '<g fill="none" stroke-linecap="round" stroke-linejoin="round">']
# short pieces, each its own colour along the ramp; round caps join them
STEP = 1
for a in range(0, n-1, STEP):
    b = min(n-1, a+STEP)
    seg = [k for k in range(a, b+1) if not hidden[k]]
    if len(seg) < 2 or seg[-1]-seg[0] != len(seg)-1: 
        # split around hidden points
        runs, cur = [], []
        for k in range(a, b+1):
            if hidden[k]:
                if len(cur) > 1: runs.append(cur)
                cur = []
            else: cur.append(k)
        if len(cur) > 1: runs.append(cur)
    else:
        runs = [seg]
    for r in runs:
        o.append(f'<polyline points="' + ' '.join(f'{pts[k][0]:.1f},{pts[k][1]:.1f}' for k in r)
                 + f'" stroke="{ramp((r[0]+r[-1])/2/(n-1))}" stroke-width="{wd[(r[0]+r[-1])//2]:.1f}"/>')
o.append('</g></svg>')
open(sys.argv[1], 'w').write('\n'.join(o))
