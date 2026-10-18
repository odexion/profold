#!/usr/bin/env bash
# Draw the logo and every icon made from it, on macOS:
#
#   scripts/icons.sh
#
# From assets/logo.py:
#   logo.svg, logo-small.svg  the logo on its tile, and the glyph for small sizes
#   mark.svg, mark.png        the bare glyph, for the app's masthead
#   icon.png                  the Dock icon for an unbundled build, and Linux's
#   profold.icns              the bundle icon
#   profold.ico               the Windows executable's icon
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root/assets"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

python3 logo.py logo.svg
python3 logo.py --small logo-small.svg
python3 logo.py --small --bare mark.svg
swiftc -O "$root/scripts/raster.swift" -o "$work/raster"
r="$work/raster"

"$r" logo.svg icon.png 1024
"$r" mark.svg mark.png 64      # twice the masthead's 32 points, so it is drawn 1:1

set="$work/profold.iconset"
mkdir "$set"
for s in 16 32; do
    "$r" logo-small.svg "$set/icon_${s}x${s}.png" "$s"
    "$r" logo-small.svg "$set/icon_${s}x${s}@2x.png" $((s*2))
done
for s in 128 256 512; do
    "$r" logo.svg "$set/icon_${s}x${s}.png" "$s"
    "$r" logo.svg "$set/icon_${s}x${s}@2x.png" $((s*2))
done
iconutil -c icns "$set" -o profold.icns

# an .ico of PNGs, smallest first
"$r" logo.svg "$work/48.png" 48
python3 - "$set/icon_16x16.png" "$set/icon_32x32.png" "$work/48.png" "$set/icon_256x256.png" <<'EOF'
import struct, sys
pngs = [open(f, 'rb').read() for f in sys.argv[1:]]
off = 6 + 16*len(pngs); entries = b''
for s, d in zip((16, 32, 48, 256), pngs):
    entries += struct.pack('<BBBBHHII', s % 256, s % 256, 0, 0, 1, 32, len(d), off); off += len(d)
open('profold.ico', 'wb').write(struct.pack('<HHH', 0, 1, len(pngs)) + entries + b''.join(pngs))
EOF
