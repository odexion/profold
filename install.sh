#!/bin/sh
# The profold installer for macOS and Linux.
#
#   curl -LsSf https://github.com/odexion/profold/releases/latest/download/install.sh | sh
#
# Installs for the current user, with no sudo:
#   macOS  ~/Applications/profold.app, and `profold` in ~/.local/bin
#   Linux  ~/.local/bin/profold, with a launcher in the applications menu
#
# Settings, all optional:
#   PROFOLD_VERSION=0.2.0     a particular release rather than the latest
#   PROFOLD_BIN_DIR=DIR       where the `profold` command goes
#   PROFOLD_APP_DIR=DIR       macOS: where profold.app goes
#   PROFOLD_NO_MODIFY_PATH=1  leave shell profiles alone
#   PROFOLD_DOWNLOAD_BASE=URL fetch the archives from a mirror
#   NO_COLOR=1                plain output

set -eu

REPO="odexion/profold"

# --- output -----------------------------------------------------------------

if [ -t 2 ] && [ -z "${NO_COLOR:-}" ] && [ "${TERM:-}" != "dumb" ]; then
    BOLD=$(printf '\033[1m'); DIM=$(printf '\033[2m'); RESET=$(printf '\033[0m')
    RED=$(printf '\033[31m'); GREEN=$(printf '\033[32m'); YELLOW=$(printf '\033[33m')
    MAGENTA=$(printf '\033[35m'); CYAN=$(printf '\033[36m')
    FANCY=1
else
    BOLD=""; DIM=""; RESET=""; RED=""; GREEN=""; YELLOW=""; MAGENTA=""; CYAN=""
    FANCY=""
fi

say() { printf '%s\n' "$*" >&2; }
step() { say "${BOLD}${CYAN}==>${RESET} ${BOLD}$*${RESET}"; }
ok() { say "  ${GREEN}✓${RESET} $*"; }
warn() { say "  ${YELLOW}!${RESET} $*"; }
fail() {
    say "${RED}${BOLD}error${RESET}${BOLD}:${RESET} $*"
    exit 1
}

need() {
    command -v "$1" >/dev/null 2>&1 || fail "this installer needs \`$1\`, which was not found"
}

# --- what to install ----------------------------------------------------------

detect_target() {
    os=$(uname -s)
    arch=$(uname -m)
    case "$os" in
        Darwin)
            # A shell under Rosetta reports x86_64 on an Apple silicon Mac.
            if [ "$arch" = "x86_64" ] && [ "$(sysctl -n sysctl.proc_translated 2>/dev/null || echo 0)" = "1" ]; then
                arch=arm64
            fi
            case "$arch" in
                arm64 | aarch64) echo "aarch64-apple-darwin" ;;
                x86_64) echo "x86_64-apple-darwin" ;;
                *) fail "unsupported Mac architecture: $arch" ;;
            esac
            ;;
        Linux)
            if ldd --version 2>&1 | grep -qi musl; then
                fail "profold is built for glibc Linux; musl systems such as Alpine are not supported yet"
            fi
            case "$arch" in
                x86_64 | amd64) echo "x86_64-unknown-linux-gnu" ;;
                aarch64 | arm64) echo "aarch64-unknown-linux-gnu" ;;
                *) fail "unsupported Linux architecture: $arch" ;;
            esac
            ;;
        MINGW* | MSYS* | CYGWIN*)
            fail "on Windows, run this in PowerShell instead:
    powershell -ExecutionPolicy ByPass -c \"irm https://github.com/$REPO/releases/latest/download/install.ps1 | iex\""
            ;;
        *) fail "unsupported system: $os" ;;
    esac
}

# --- downloading, with a bar ----------------------------------------------------

# bar HAVE TOTAL ELAPSED LABEL
bar() {
    awk -v have="$1" -v total="$2" -v secs="$3" -v label="$4" \
        -v on="$MAGENTA" -v off="$DIM" -v reset="$RESET" -v bold="$BOLD" '
    function human(b,   u, i) {
        split("B KiB MiB GiB", u, " "); i = 1
        while (b >= 1024 && i < 4) { b /= 1024; i++ }
        return sprintf(i == 1 ? "%d %s" : "%.1f %s", b, u[i])
    }
    BEGIN {
        width = 30
        rate = secs > 0 ? have / secs : 0
        speed = rate > 0 ? "  " human(rate) "/s" : ""
        if (total > 0) {
            f = have / total; if (f > 1) f = 1
            n = int(f * width)
            s = ""
            for (i = 0; i < n; i++) s = s "━"
            tail = ""
            if (n < width) { tail = "╺"; for (i = n + 1; i < width; i++) tail = tail "━" }
            printf "\r\033[2K  %s%s%s  %s%s%s%s%s%s  %3d%%  %s / %s%s", bold, label, reset, on, s, reset, off, tail, reset, f * 100, human(have), human(total), speed
        } else {
            printf "\r\033[2K  %s%s%s  %s%s", bold, label, reset, human(have), speed
        }
    }' >&2
}

# download URL FILE LABEL
download() {
    url=$1; out=$2; label=$3
    if [ -z "$FANCY" ]; then
        say "  $label…"
        curl -fsSL --retry 3 -o "$out" "$url" </dev/null || return 1
        return 0
    fi
    headers="$out.headers"
    : >"$out"
    : >"$headers"
    curl -fsSL --retry 3 -D "$headers" -o "$out" "$url" </dev/null &
    pid=$!
    started=$(date +%s)
    while kill -0 "$pid" 2>/dev/null; do
        # The last Content-Length is the file's own, after the redirects.
        total=$(tr -d '\r' <"$headers" | awk 'tolower($1) == "content-length:" { n = $2 } END { print n + 0 }')
        have=$(wc -c <"$out" | tr -d ' ')
        bar "$have" "$total" "$(( $(date +%s) - started ))" "$label"
        sleep 0.1
    done
    if ! wait "$pid"; then
        printf '\r\033[2K' >&2
        rm -f "$headers"
        return 1
    fi
    have=$(wc -c <"$out" | tr -d ' ')
    bar "$have" "$have" "$(( $(date +%s) - started ))" "$label"
    printf '\n' >&2
    rm -f "$headers"
}

sha256() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{ print $1 }'
    else
        shasum -a 256 "$1" | awk '{ print $1 }'
    fi
}

# --- the shell's PATH -------------------------------------------------------------

add_to_path() {
    dir=$1
    case ":$PATH:" in *":$dir:"*) return 0 ;; esac
    if [ -n "${PROFOLD_NO_MODIFY_PATH:-}" ]; then
        warn "$dir is not on your PATH; add it to run \`profold\` from a terminal"
        return 0
    fi
    line="export PATH=\"$dir:\$PATH\""
    touched=""
    for rc in "$HOME/.profile" "$HOME/.bashrc" "$HOME/.zshrc"; do
        case "$rc" in
            */.bashrc) [ -f "$rc" ] || continue ;;
            */.zshrc) [ -f "$rc" ] || [ "$(basename "${SHELL:-}")" = "zsh" ] || continue ;;
        esac
        if ! grep -qsF "$line" "$rc"; then
            printf '\n# Added by the profold installer\n%s\n' "$line" >>"$rc"
            touched="$touched ${rc#"$HOME"/}"
        fi
    done
    if [ -n "$touched" ]; then
        ok "added ${BOLD}$dir${RESET} to PATH in ~/${touched# }"
        warn "open a new terminal, or run: ${BOLD}$line${RESET}"
    fi
}

# --- installing -------------------------------------------------------------------

install_macos() {
    unpacked=$1
    app_dir=${PROFOLD_APP_DIR:-$HOME/Applications}
    bin_dir=${PROFOLD_BIN_DIR:-$HOME/.local/bin}
    [ -d "$unpacked/profold.app" ] || fail "the download has no profold.app in it"
    mkdir -p "$app_dir" "$bin_dir"

    app="$app_dir/profold.app"
    if [ -e "$app" ]; then
        rm -rf "$app.old"
        mv "$app" "$app.old"
    fi
    mv "$unpacked/profold.app" "$app"
    rm -rf "$app.old"
    # Downloaded with curl, so there is no quarantine flag for Gatekeeper to
    # act on; clear it anyway in case something added one.
    xattr -dr com.apple.quarantine "$app" 2>/dev/null || true
    # So Spotlight and Launchpad find it straight away.
    lsregister=/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister
    if [ -x "$lsregister" ]; then
        "$lsregister" -f "$app" >/dev/null 2>&1 || true
    fi
    ok "installed ${BOLD}$app${RESET}"

    cat >"$bin_dir/profold" <<EOF
#!/bin/sh
# Starts profold, installed at $app.
exec "$app/Contents/MacOS/profold" "\$@"
EOF
    chmod 755 "$bin_dir/profold"
    ok "linked the ${BOLD}profold${RESET} command in $bin_dir"
    add_to_path "$bin_dir"
    LAUNCH="open it from Applications or Spotlight, or run ${BOLD}profold${RESET}"
}

install_linux() {
    unpacked=$1
    bin_dir=${PROFOLD_BIN_DIR:-$HOME/.local/bin}
    data=${XDG_DATA_HOME:-$HOME/.local/share}
    [ -f "$unpacked/profold" ] || fail "the download has no profold binary in it"
    mkdir -p "$bin_dir"

    # Renamed into place, so a running copy keeps its file until it exits.
    cp "$unpacked/profold" "$bin_dir/.profold.new"
    chmod 755 "$bin_dir/.profold.new"
    mv -f "$bin_dir/.profold.new" "$bin_dir/profold"
    ok "installed ${BOLD}$bin_dir/profold${RESET}"

    if [ -f "$unpacked/profold.png" ]; then
        mkdir -p "$data/icons/hicolor/256x256/apps"
        cp "$unpacked/profold.png" "$data/icons/hicolor/256x256/apps/dev.profold.png"
    fi
    mkdir -p "$data/applications"
    # Named for the app id the window reports, so the launcher and the
    # window are one in the dock.
    cat >"$data/applications/dev.profold.desktop" <<EOF
[Desktop Entry]
Type=Application
Name=profold
GenericName=Molecular viewer
Comment=View proteins, molecules, crystals and simulations
Exec="$bin_dir/profold" %F
Icon=dev.profold
Terminal=false
Categories=Science;Chemistry;Biology;Education;
Keywords=protein;molecule;PDB;structure;crystal;trajectory;
StartupWMClass=dev.profold
EOF
    if command -v update-desktop-database >/dev/null 2>&1; then
        update-desktop-database -q "$data/applications" 2>/dev/null || true
    fi
    ok "added profold to the applications menu"

    # GPUI loads these at run time; say which are missing rather than
    # leaving a window that never opens.
    if command -v ldd >/dev/null 2>&1; then
        missing=$(ldd "$bin_dir/profold" 2>/dev/null | awk '/not found/ { print $1 }' | tr '\n' ' ')
        if [ -n "$missing" ]; then
            warn "missing system libraries: ${BOLD}$missing${RESET}"
            warn "Debian/Ubuntu: sudo apt install libxkbcommon-x11-0 libvulkan1 libwayland-client0 libxcb1 libfontconfig1"
            warn "Fedora:        sudo dnf install libxkbcommon-x11 vulkan-loader wayland libxcb fontconfig"
        fi
    fi
    add_to_path "$bin_dir"
    LAUNCH="open it from your applications menu, or run ${BOLD}profold${RESET}"
}

main() {
    need curl
    need tar
    need uname

    say ""
    say "  ${BOLD}${MAGENTA}profold${RESET}${BOLD} installer${RESET}  ${DIM}a molecular viewer · github.com/$REPO${RESET}"
    say ""

    target=$(detect_target)
    if [ -n "${PROFOLD_VERSION:-}" ]; then
        version=${PROFOLD_VERSION#v}
        base="https://github.com/$REPO/releases/download/v$version"
    else
        version=""
        base="https://github.com/$REPO/releases/latest/download"
    fi
    base=${PROFOLD_DOWNLOAD_BASE:-$base}
    archive="profold-$target.tar.gz"

    step "Installing profold${version:+ $version} for ${BOLD}$target${RESET}"
    tmp=$(mktemp -d 2>/dev/null || mktemp -d -t profold)
    trap 'rm -rf "$tmp"' EXIT INT TERM

    download "$base/$archive" "$tmp/$archive" "Downloading" ||
        fail "could not download $base/$archive
  Is there a release for $target? See https://github.com/$REPO/releases"
    curl -fsSL --retry 3 -o "$tmp/$archive.sha256" "$base/$archive.sha256" </dev/null ||
        fail "could not download the checksum for $archive"

    want=$(awk '{ print $1; exit }' "$tmp/$archive.sha256")
    got=$(sha256 "$tmp/$archive")
    [ "$want" = "$got" ] || fail "the download is damaged (checksum mismatch); try again"
    ok "checksum verified"

    mkdir -p "$tmp/unpacked"
    tar -xzf "$tmp/$archive" -C "$tmp/unpacked"

    case "$target" in
        *apple-darwin) install_macos "$tmp/unpacked" ;;
        *) install_linux "$tmp/unpacked" ;;
    esac

    say ""
    say "  ${GREEN}${BOLD}profold is installed!${RESET} To start, $LAUNCH."
    say "  ${DIM}It keeps itself up to date: new versions arrive with a \"Restart to update\" button.${RESET}"
    say ""
}

main "$@"
