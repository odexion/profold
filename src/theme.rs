//! The palette.
//!
//! Paper white, ink black, one weight of hairline. The only saturated colour in
//! the app is the molecule itself.

use gpui::Hsla;

/// Window background.
pub const PAPER: Hsla = Hsla { h: 0.1250, s: 0.3333, l: 0.9765, a: 1.00 };

/// Cards and panels.
pub const SURFACE: Hsla = Hsla { h: 0.0000, s: 0.0000, l: 1.0000, a: 1.00 };

/// Input wells and inactive chips.
pub const SUNKEN: Hsla = Hsla { h: 0.1333, s: 0.1852, l: 0.9471, a: 1.00 };

/// Hairlines.
pub const BORDER: Hsla = Hsla { h: 0.1111, s: 0.1579, l: 0.8882, a: 1.00 };

/// Focused field outline.
pub const BORDER_FOCUS: Hsla = Hsla { h: 0.6667, s: 0.0833, l: 0.0941, a: 1.00 };

/// Primary text, and the fill of primary buttons.
pub const INK: Hsla = Hsla { h: 0.6667, s: 0.0833, l: 0.0941, a: 1.00 };

/// Body text.
pub const TEXT: Hsla = Hsla { h: 0.6667, s: 0.0617, l: 0.1588, a: 1.00 };

/// Secondary text.
pub const MUTED: Hsla = Hsla { h: 0.1061, s: 0.0435, l: 0.5039, a: 1.00 };

/// Placeholders and disabled text.
pub const FAINT: Hsla = Hsla { h: 0.1111, s: 0.0750, l: 0.6863, a: 1.00 };

/// Text caret.
pub const CARET: Hsla = Hsla { h: 0.6667, s: 0.0833, l: 0.0941, a: 1.00 };

/// Text selection.
pub const SELECTION: Hsla = Hsla { h: 0.6667, s: 0.0833, l: 0.0941, a: 0.14 };

/// Errors.
pub const DANGER: Hsla = Hsla { h: 0.0092, s: 0.7005, l: 0.4059, a: 1.00 };

/// Text on an ink-filled surface.
pub const ON_INK: Hsla = Hsla { h: 0.1250, s: 0.3333, l: 0.9765, a: 1.00 };


/// The page behind the molecule, in 0..1 RGB for the rasterizer. Both styles
/// sit on the same flat paper, so switching between them changes the drawing
/// and nothing else.
pub const VIEW: [f32; 3] = [0.992, 0.992, 0.988];

/// Monospace family, used for sequences and identifiers: one that ships
/// with each system.
#[cfg(target_os = "macos")]
pub const MONO: &str = "Menlo";
#[cfg(windows)]
pub const MONO: &str = "Consolas";
#[cfg(not(any(target_os = "macos", windows)))]
pub const MONO: &str = "DejaVu Sans Mono";

/// Display serif, used only for the name of the structure on screen.
#[cfg(target_os = "macos")]
pub const SERIF: &str = "Charter";
#[cfg(windows)]
pub const SERIF: &str = "Cambria";
#[cfg(not(any(target_os = "macos", windows)))]
pub const SERIF: &str = "DejaVu Serif";

/// A shortcut as this system writes it: `⇧⌘T` on a Mac, `Ctrl+Shift+T`
/// elsewhere. Written the Mac way, with ⇧ before ⌘.
pub fn keys(mac: &str) -> String {
    if cfg!(target_os = "macos") {
        mac.to_string()
    } else {
        spelled_out(mac)
    }
}

#[cfg_attr(target_os = "macos", allow(dead_code))]
fn spelled_out(mac: &str) -> String {
    let mut parts = Vec::new();
    let mut rest = mac;
    for (glyph, name) in [("⌃", "Ctrl"), ("⌥", "Alt"), ("⇧", "Shift"), ("⌘", "Ctrl")] {
        if let Some(r) = rest.strip_prefix(glyph) {
            parts.push(name);
            rest = r;
        }
    }
    // ⌘ comes last on a Mac but Ctrl comes first elsewhere.
    if parts.last() == Some(&"Ctrl") && parts.len() > 1 {
        parts.pop();
        parts.insert(0, "Ctrl");
    }
    parts.push(rest);
    parts.join("+")
}

#[cfg(test)]
mod tests {
    #[test]
    fn shortcuts_are_spelled_out_off_a_mac() {
        assert_eq!(super::spelled_out("⌘T"), "Ctrl+T");
        assert_eq!(super::spelled_out("⇧⌘T"), "Ctrl+Shift+T");
        assert_eq!(super::spelled_out("⌘O"), "Ctrl+O");
    }
}
