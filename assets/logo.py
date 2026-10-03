# Draws assets/logo.svg: `python3 assets/logo.py assets/logo.svg`.
#
# Two strands winding round each other, and nothing else. One violet on
# paper, the strands thinning and paling as they turn away.
import math, sys

PAPER = '#fdfdfc'
INK = '#2b1b40'
TURNS = 1.5                       # how much helix is shown
PITCH_LEN = 4.6                   # axis length of one turn, in helix radii
W = 0.42                          # strand width, in helix radii
THIN = 0.35                       # width of the far side, as a fraction of the near
FADE = 0.55                       # how far the far side pales toward the paper
TILT = -24.0                      # lean of the axis on the tile
PH = 0.0                          # each strand ends at the widest point of a turn
BOX = 690.0

N = 900
def strand(s, phase):
    th = s*TURNS*2*math.pi + phase
    return (math.cos(th), (s - 0.5)*TURNS*PITCH_LEN, math.sin(th))
def view(p):
    x, y, z = p
    t = math.radians(TILT)
    return (x*math.cos(t) - y*math.sin(t), -(x*math.sin(t) + y*math.cos(t)), z)

S = [[view(strand(i/N, PH + k*math.pi)) for i in range(N+1)] for k in (0, 1)]
allp = S[0] + S[1]
x0 = min(p[0] for p in allp) - W/2; x1 = max(p[0] for p in allp) + W/2
y0 = min(p[1] for p in allp) - W/2; y1 = max(p[1] for p in allp) + W/2
sc = BOX/max(x1-x0, y1-y0)
S = [[((p[0]-(x0+x1)/2)*sc + 512, (p[1]-(y0+y1)/2)*sc + 512, p[2]) for p in P] for P in S]
w = W*sc

# Depth is the only shading: a strand swells and darkens as it comes toward
# the eye and thins and pales as it goes away. Short round-capped pieces,
# drawn far to near, so the near strand passes over the far one.
PAP, INKC = (0xfd, 0xfd, 0xfc), tuple(int(INK[k:k+2], 16) for k in (1, 3, 5))
def look(z):
    t = (z + 1)/2                                  # 0 far, 1 near
    width = w*(THIN + (1 - THIN)*t)
    f = FADE*(1 - t)
    c = '#%02x%02x%02x' % tuple(round(INKC[k] + (PAP[k]-INKC[k])*f) for k in range(3))
    return width, c
pieces = []
for P in S:
    for i in range(len(P)-1):
        pieces.append(((P[i][2]+P[i+1][2])/2, P[i], P[i+1]))
pieces.sort(key=lambda q: q[0])

o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1024 1024" width="1024" height="1024">',
 '<defs><filter id="sh" x="-20%" y="-20%" width="140%" height="140%">'
 '<feDropShadow dx="0" dy="10" stdDeviation="14" flood-color="#1f1d17" flood-opacity="0.18"/></filter></defs>',
 f'<rect x="100" y="100" width="824" height="824" rx="185" fill="{PAPER}" filter="url(#sh)"/>',
 '<g fill="none" stroke-linecap="round">']
for z, p, q in pieces:
    width, c = look(z)
    o.append(f'<line x1="{p[0]:.1f}" y1="{p[1]:.1f}" x2="{q[0]:.1f}" y2="{q[1]:.1f}" stroke="{c}" stroke-width="{width:.1f}"/>')
o.append('</g></svg>')
open(sys.argv[1], 'w').write('\n'.join(o))
