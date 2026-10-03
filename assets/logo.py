# Draws assets/logo.svg: `python3 assets/logo.py assets/logo.svg`.
#
# A DNA double helix in cartoon: two backbone tubes winding round each other,
# the base pairs as rungs between them. One violet, a dark hairline, the far
# side a little faded toward the paper.
import math, sys

PAPER = '#fdfdfc'
INK = (0x13, 0x13, 0x1b)
TUBE = (0x3a, 0x28, 0x56)
RUNG = (0x5c, 0x48, 0x7c)
SHEEN = (0x6d, 0x58, 0x8f)
TURNS = 1.4                       # how much helix is shown
PITCH_LEN = 3.9                   # axis length of one turn, in helix radii
PAIRS_PER_TURN = 10
W = 0.36                          # backbone width, in helix radii
RW = 0.17                         # rung width
EDGE = 0.045                      # outline
FADE = 0.30                       # how far the back fades toward the paper
TILT = -34.0                      # lean of the axis on the tile
BOX = 690.0

N = 900
def strand(s, phase):
    th = s*TURNS*2*math.pi + phase
    return (math.cos(th), (s - 0.5)*TURNS*PITCH_LEN, math.sin(th))
def view(p):
    x, y, z = p
    t = math.radians(TILT)
    return (x*math.cos(t) - y*math.sin(t), -(x*math.sin(t) + y*math.cos(t)), z)

PH = math.radians(108)         # ends just short of where the strands cross
A = [view(strand(i/N, PH)) for i in range(N+1)]
B = [view(strand(i/N, PH + math.pi)) for i in range(N+1)]
pairs = int(TURNS*PAIRS_PER_TURN)
rungs = []
for j in range(pairs):
    s = (j + 0.5)/pairs
    rungs.append((view(strand(s, PH)), view(strand(s, PH + math.pi))))

allp = A + B
x0 = min(p[0] for p in allp) - W/2 - EDGE; x1 = max(p[0] for p in allp) + W/2 + EDGE
y0 = min(p[1] for p in allp) - W/2 - EDGE; y1 = max(p[1] for p in allp) + W/2 + EDGE
sc = BOX/max(x1-x0, y1-y0)
def px(p):
    return ((p[0]-(x0+x1)/2)*sc + 512, (p[1]-(y0+y1)/2)*sc + 512, p[2])
A = [px(p) for p in A]; B = [px(p) for p in B]
rungs = [(px(a), px(b)) for a, b in rungs]
w, rw, e = W*sc, RW*sc, EDGE*sc

PAP = (0xfd, 0xfd, 0xfc)
def tint(c, z):
    f = FADE*(1 - (z + 1)/2)
    return '#%02x%02x%02x' % tuple(round(c[k] + (PAP[k]-c[k])*f) for k in range(3))
def pl(P, a, b, dx=0, dy=0):
    return ' '.join(f'{P[k][0]+dx:.1f},{P[k][1]+dy:.1f}' for k in range(a, b+1))

o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1024 1024" width="1024" height="1024">',
 '<defs><filter id="sh" x="-20%" y="-20%" width="140%" height="140%">'
 '<feDropShadow dx="0" dy="10" stdDeviation="14" flood-color="#1f1d17" flood-opacity="0.18"/></filter></defs>',
 f'<rect x="100" y="100" width="824" height="824" rx="185" fill="{PAPER}" filter="url(#sh)"/>',
 '<g fill="none" stroke-linejoin="round" stroke-linecap="butt">']

SUB = 4
def tube(P, s, t, ext, caps):
    o.append(f'<polyline points="{pl(P, s, t)}" stroke="{tint(INK, P[(s+t)//2][2])}" stroke-width="{w+2*e:.1f}" stroke-linecap="{caps}"/>')
    fs, ft = max(0, s-ext), min(len(P)-1, t+ext)
    for k in (fs, ft):
        if k in (0, len(P)-1):
            o.append(f'<circle cx="{P[k][0]:.1f}" cy="{P[k][1]:.1f}" r="{w/2:.1f}" fill="{tint(TUBE, P[k][2])}"/>')
    for a in range(fs, ft, SUB):
        b = min(ft, a+SUB+1)
        o.append(f'<polyline points="{pl(P, a, b)}" stroke="{tint(TUBE, P[(a+b)//2][2])}" stroke-width="{w:.1f}"/>')
    for a in range(fs+2, ft-2, SUB):
        b = min(ft-2, a+SUB+1)
        o.append(f'<polyline points="{pl(P, a, b, -w*0.13, -w*0.16)}" stroke="{tint(SHEEN, P[(a+b)//2][2])}" stroke-width="{w*0.22:.1f}"/>')

def runs(P):
    out, s0 = [], 0
    for k in range(1, len(P)):
        if (P[k][2] >= 0) != (P[s0][2] >= 0) or k == len(P)-1:
            out.append((s0, k)); s0 = k
    return [(s, t) for s, t in out if P[(s+t)//2][2] > 0]

# Both backbones whole, then the rungs, then each backbone's near half again
# over the top. Every join lands on tube already drawn, so none shows.
tube(A, 0, N, 0, 'round')
tube(B, 0, N, 0, 'round')
for a, b in rungs:
    d = math.hypot(b[0]-a[0], b[1]-a[1])
    if d < w*1.6: continue
    z = (a[2]+b[2])/2
    # stop at the tubes' edges, so a rung tucks in under the far strand too
    k = (w/2 + e)/d
    ax, ay = a[0] + (b[0]-a[0])*k, a[1] + (b[1]-a[1])*k
    bx, by = b[0] - (b[0]-a[0])*k, b[1] - (b[1]-a[1])*k
    seg = f'x1="{ax:.1f}" y1="{ay:.1f}" x2="{bx:.1f}" y2="{by:.1f}"'
    o.append(f'<line {seg} stroke="{tint(INK, z)}" stroke-width="{rw+2*e:.1f}" stroke-linecap="butt"/>')
    o.append(f'<line {seg} stroke="{tint(RUNG, z)}" stroke-width="{rw:.1f}" stroke-linecap="butt"/>')
for P in (A, B):
    for s, t in runs(P):
        tube(P, s, t, 14, 'butt')
o.append('</g></svg>')
open(sys.argv[1], 'w').write('\n'.join(o))
