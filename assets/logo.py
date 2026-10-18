# Draws the logo: `python3 assets/logo.py assets/logo.svg`, and the glyph for
# the smallest icon sizes with `python3 assets/logo.py --small assets/logo-small.svg`.
# `--bare` leaves off the tile and crops to the strand, for the app's own
# masthead: `python3 assets/logo.py --small --bare assets/mark.svg`.
# scripts/icons.sh draws all three and renders the icons from them.
#
# A protein cartoon: a ribbon wound into a helix that leaves as a strand and
# ends in an arrowhead, coloured along its length from the N-terminus to the C.
# The ribbon's outer face is lit, its inner face is shaded, and a dark edge
# keeps one turn apart from the next where they cross. The small glyph is the
# same fold in a wider ribbon and a heavier line, so it still reads at 16 pixels.
import math, sys

SMALL = '--small' in sys.argv
BARE = '--bare' in sys.argv
OUT = [s for s in sys.argv[1:] if not s.startswith('--')][0]

TILE = ('#ffffff', '#f1f2f6')     # top and bottom of the paper
INK = (8, 6, 20)                  # what the inner face and the edge darken toward
# the ramp along the chain, start to end
RAMP = [(0x43, 0x38, 0xca), (0x7c, 0x3a, 0xed), (0xc0, 0x26, 0xd3), (0xf4, 0x3f, 0x5e), (0xfb, 0x92, 0x3c)]
TURNS = 2.5
RISE = 2.5 if SMALL else 2.1      # along the axis per turn, in radii
W = 1.45 if SMALL else 1.1        # ribbon width, in radii
EDGE = 9.0 if SMALL else 5.0      # edge line, in pixels
BOX = 660.0 if SMALL else 620.0
PHASE = 0.7                       # where the helix starts, in half turns
EXIT_ANGLE = 5.0                  # the strand, in degrees above the page's horizontal
EXIT_LEN = 1.2                    # and its length, in radii
HANDLE = 0.8                      # how gently it bends out of the helix

def ramp(t):
    t = max(0.0, min(1.0, t))*(len(RAMP)-1)
    i = min(int(t), len(RAMP)-2); f = t - i
    return tuple(RAMP[i][k] + (RAMP[i+1][k]-RAMP[i][k])*f for k in range(3))
def mix(c, d, f): return tuple(c[k] + (d[k]-c[k])*f for k in range(3))
def hexc(c): return '#%02x%02x%02x' % tuple(max(0, min(255, round(v))) for v in c)
def add(p, q): return tuple(p[k]+q[k] for k in range(3))
def sub(p, q): return tuple(p[k]-q[k] for k in range(3))
def mul(p, s): return tuple(v*s for v in p)
def dot(p, q): return sum(p[k]*q[k] for k in range(3))
def cross(p, q): return (p[1]*q[2]-p[2]*q[1], p[2]*q[0]-p[0]*q[2], p[0]*q[1]-p[1]*q[0])
def unit(p): return mul(p, 1/math.sqrt(dot(p, p)))

def view(p):
    x, y, z = p
    a = math.radians(24)     # tip the axis toward the eye
    y, z = y*math.cos(a) - z*math.sin(a), y*math.sin(a) + z*math.cos(a)
    b = math.radians(-26)                     # lay it on the diagonal
    x, y = x*math.cos(b) - y*math.sin(b), x*math.sin(b) + y*math.cos(b)
    return (x, -y, z)
# the page's right, up and toward the eye, in the chain's own frame
BASIS = [view(v) for v in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
RIGHT = tuple(v[0] for v in BASIS); UP = tuple(-v[1] for v in BASIS); EYE = tuple(v[2] for v in BASIS)

# The centreline, with the direction the ribbon's width lies along. On the
# helix that is the axis, kept on the cylinder; past it, the ribbon turns to
# lie flat to the eye, so the strand and its arrowhead are seen face on.
H = TURNS*RISE
PH = math.pi*PHASE
C, across = [], []
NH = int(260*TURNS)
for i in range(NH+1):
    th = i/NH*TURNS*2*math.pi + PH
    C.append((math.cos(th), i/NH*H, -math.sin(th)))   # right-handed, as an α-helix is
    t = unit((-math.sin(th)*TURNS*2*math.pi, H, -math.cos(th)*TURNS*2*math.pi))
    across.append(unit(sub((0, 1, 0), mul(t, t[1]))))
# out of the helix along its tangent, and round into a straight strand
tH = unit(sub(C[-1], C[-2]))
# the strand points up and to the right on the page, flat to the eye
EXIT = math.radians(EXIT_ANGLE)
tS = unit(add(mul(RIGHT, math.cos(EXIT)), mul(UP, math.sin(EXIT))))
P0 = C[-1]; P1 = add(P0, mul(tH, HANDLE))
P3 = add(P0, mul(tS, EXIT_LEN))
P2 = sub(P3, mul(tS, HANDLE))
for i in range(1, 71):
    u = i/70
    C.append(tuple((1-u)**3*P0[k] + 3*(1-u)**2*u*P1[k] + 3*(1-u)*u*u*P2[k] + u**3*P3[k] for k in range(3)))
    across.append(None)
SL = 0.0
for i in range(1, 51 if SL else 1):
    C.append(add(P3, mul(tS, SL*i/50))); across.append(None)
last = across[NH]
for k in range(NH+1, len(C)):
    t = unit(sub(C[k], C[k-1]))
    flat = unit(cross(t, EYE))
    if dot(flat, last) < 0: flat = mul(flat, -1)
    e = min(1.0, (k-NH)/60); e = e*e*(3-2*e)
    b = add(mul(last, 1-e), mul(flat, e))
    across[k] = unit(sub(b, mul(t, dot(b, t))))
body = len(C)
# the arrowhead: the same point twice, narrow then wide, makes its shoulders
AL, NA = 0.8, 30
tE = unit(sub(C[-1], C[-2])); E = C[-1]
wd = [W]*body
C.append(E); across.append(across[-1]); wd.append(1.6*W)
for i in range(1, NA+1):
    C.append(add(E, mul(tE, AL*i/NA))); across.append(across[body-1]); wd.append(1.6*W*(1 - i/NA) + 1e-3)
n = len(C)

L = [view(add(C[k], mul(across[k], wd[k]/2))) for k in range(n)]
R = [view(sub(C[k], mul(across[k], wd[k]/2))) for k in range(n)]
x0 = min(p[0] for p in L+R); x1 = max(p[0] for p in L+R)
y0 = min(p[1] for p in L+R); y1 = max(p[1] for p in L+R)
sc = BOX/max(x1-x0, y1-y0)
def page(p): return ((p[0]-(x0+x1)/2)*sc + 512, (p[1]-(y0+y1)/2)*sc + 518, p[2])
L = [page(p) for p in L]; R = [page(p) for p in R]

# One quad per step along the chain, painted far to near.
LIGHT = unit((-0.4, 0.6, 1.0))    # from the upper left, toward the eye
quads = []
for k in range(n-1):
    l0, l1, r1, r0 = L[k], L[k+1], R[k+1], R[k]
    u = sub(l1, l0) if l1 != l0 else sub(r1, r0)
    nn = cross(u, sub(r0, l0))
    nn = unit(nn) if dot(nn, nn) > 1e-9 else (0, 0, 1)
    base = ramp(min(k, body-1)/(body-1))
    if nn[2] >= 0:                # the outer face
        lit = abs(dot(nn, LIGHT))
        col = mix(mix(base, (255, 255, 255), 0.22*lit**6), (0, 0, 0), 0.10*(1-lit))
    else:                         # the inner face
        col = mix(base, INK, 0.38)
    quads.append(((l0[2]+l1[2]+r0[2]+r1[2])/4, k, col, mix(base, INK, 0.55)))
quads.sort()

V = (512 - BOX/2 - EDGE, 518 - BOX/2 - EDGE, BOX + 2*EDGE) if BARE else (0, 0, 1024)
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{V[0]:.0f} {V[1]:.0f} {V[2]:.0f} {V[2]:.0f}" width="1024" height="1024">',
 '<defs><linearGradient id="bg" x1="0" y1="0" x2="0" y2="1">'
 f'<stop offset="0" stop-color="{TILE[0]}"/><stop offset="1" stop-color="{TILE[1]}"/></linearGradient></defs>']
# the shadow as rings of faint tile, not a filter, so every renderer draws it
for i in range(0 if BARE else 12, 0, -1):
    g = i*2.5
    o.append(f'<rect x="{100-g:.1f}" y="{100-g+8:.1f}" width="{824+2*g:.1f}" height="{824+2*g:.1f}" '
             f'rx="{185+g:.1f}" fill="#0f0a1e" fill-opacity="0.009"/>')
if not BARE:
    o.append('<rect x="100" y="100" width="824" height="824" rx="185" fill="url(#bg)"/>')
o += ['<g stroke-linejoin="round" stroke-linecap="round">']
xy = lambda p: f'{p[0]:.1f},{p[1]:.1f}'
for _, k, col, edge in quads:
    l0, l1, r1, r0 = L[k], L[k+1], R[k+1], R[k]
    c = hexc(col)
    # stroked in its own colour so neighbouring quads leave no seam
    o.append(f'<polygon points="{xy(l0)} {xy(l1)} {xy(r1)} {xy(r0)}" fill="{c}" stroke="{c}" stroke-width="1.2"/>')
    d = f'M{xy(l0)}L{xy(l1)}M{xy(r0)}L{xy(r1)}'
    if k == 0 or l0[:2] == l1[:2]:  # the N-terminal end, and the arrow's shoulders
        d += f'M{xy(l0)}L{xy(r0)}'
    o.append(f'<path d="{d}" fill="none" stroke="{hexc(edge)}" stroke-width="{EDGE}"/>')
o.append('</g></svg>')
open(OUT, 'w').write('\n'.join(o))
