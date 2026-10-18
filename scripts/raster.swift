// Renders an SVG to a square PNG with its transparency: `raster in.svg out.png size`.
// Used by scripts/icons.sh; macOS draws the SVG itself, so nothing else is needed.
import AppKit

let a = CommandLine.arguments
guard a.count == 4, let s = Int(a[3]), let svg = NSImage(contentsOfFile: a[1]) else {
    print("usage: raster in.svg out.png size"); exit(1)
}
let rep = NSBitmapImageRep(bitmapDataPlanes: nil, pixelsWide: s, pixelsHigh: s, bitsPerSample: 8,
                           samplesPerPixel: 4, hasAlpha: true, isPlanar: false,
                           colorSpaceName: .deviceRGB, bytesPerRow: 0, bitsPerPixel: 0)!
rep.size = NSSize(width: s, height: s)
NSGraphicsContext.saveGraphicsState()
NSGraphicsContext.current = NSGraphicsContext(bitmapImageRep: rep)
NSGraphicsContext.current!.imageInterpolation = .high
svg.draw(in: NSRect(x: 0, y: 0, width: s, height: s))
NSGraphicsContext.restoreGraphicsState()
try! rep.representation(using: .png, properties: [:])!.write(to: URL(fileURLWithPath: a[2]))
