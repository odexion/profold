#!/usr/bin/env bash
# Build a release archive for one target: dist/profold-<target>.tar.gz and
# its .sha256. Run by the release workflow on each platform, and by hand:
#
#   scripts/package.sh aarch64-apple-darwin
#
# The archive holds what the installers and the updater put in place:
#   macOS    profold.app
#   Linux    profold, and profold.png for the launcher
#   Windows  profold.exe
set -euo pipefail

target=${1:?usage: scripts/package.sh <rust target>}
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"

version=$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -n1)
dist="$root/dist"
stage="$dist/stage-$target"
rm -rf "$stage"
mkdir -p "$stage"

# The target, baked in, is what lets a build replace itself (see update.rs).
PROFOLD_TARGET="$target" cargo build --release --locked --target "$target"
out="target/$target/release"

case "$target" in
    *apple-darwin)
        app="$stage/profold.app"
        mkdir -p "$app/Contents/MacOS" "$app/Contents/Resources"
        cp "$out/profold" "$app/Contents/MacOS/profold"
        cp assets/profold.icns "$app/Contents/Resources/profold.icns"
        cat >"$app/Contents/Info.plist" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleName</key><string>profold</string>
    <key>CFBundleDisplayName</key><string>profold</string>
    <key>CFBundleIdentifier</key><string>dev.profold</string>
    <key>CFBundleExecutable</key><string>profold</string>
    <key>CFBundleIconFile</key><string>profold.icns</string>
    <key>CFBundlePackageType</key><string>APPL</string>
    <key>CFBundleVersion</key><string>$version</string>
    <key>CFBundleShortVersionString</key><string>$version</string>
    <key>CFBundleInfoDictionaryVersion</key><string>6.0</string>
    <key>LSMinimumSystemVersion</key><string>11.0</string>
    <key>LSApplicationCategoryType</key><string>public.app-category.education</string>
    <key>NSHighResolutionCapable</key><true/>
    <key>NSHumanReadableCopyright</key><string>MIT License</string>
</dict>
</plist>
EOF
        # Ad hoc: enough for Apple silicon to run it. Notarizing needs a
        # Developer ID, which can be added here later.
        codesign --force --deep --sign - "$app"
        ;;
    *windows*)
        cp "$out/profold.exe" "$stage/profold.exe"
        ;;
    *)
        cp "$out/profold" "$stage/profold"
        if command -v convert >/dev/null 2>&1; then
            convert assets/icon.png -resize 256x256 "$stage/profold.png"
        else
            cp assets/icon.png "$stage/profold.png"
        fi
        ;;
esac

archive="profold-$target.tar.gz"
# No `._` AppleDouble files from macOS tar.
COPYFILE_DISABLE=1 tar -czf "$dist/$archive" -C "$stage" .
rm -rf "$stage"

cd "$dist"
if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$archive" >"$archive.sha256"
else
    shasum -a 256 "$archive" >"$archive.sha256"
fi
echo "dist/$archive ($version)"
cat "$archive.sha256"
