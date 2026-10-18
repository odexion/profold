//! The palette.
//!
//! Paper white, ink black, one weight of hairline; or, dark, the same turned
//! over. The only saturated colour in the app is the molecule itself.

use std::sync::atomic::{AtomicBool, Ordering};

use gpui::Hsla;

/// Whether the window is drawn dark. One switch for the whole app, read by
/// every colour below and by the rasterizer's page. Light unless asked,
/// whatever the system is set to.
#[cfg(not(test))]
static DARK: AtomicBool = AtomicBool::new(false);

// Tests run side by side, and one that turns the window dark must not tint
// another's render halfway through: each test thread keeps its own switch.
#[cfg(test)]
thread_local! {
    static DARK: AtomicBool = const { AtomicBool::new(false) };
}

#[cfg(not(test))]
pub fn dark() -> bool {
    DARK.load(Ordering::Relaxed)
}

#[cfg(not(test))]
pub fn set_dark(dark: bool) {
    DARK.store(dark, Ordering::Relaxed);
}

#[cfg(test)]
pub fn dark() -> bool {
    DARK.with(|d| d.load(Ordering::Relaxed))
}

#[cfg(test)]
pub fn set_dark(dark: bool) {
    DARK.with(|d| d.store(dark, Ordering::Relaxed));
}

const fn hsla(h: f32, s: f32, l: f32, a: f32) -> Hsla {
    Hsla { h, s, l, a }
}

/// Each colour as a light and a dark pair. The dark side is the same warm
/// grey turned over: charcoal paper, chalk ink.
macro_rules! palette {
    ($($(#[$doc:meta])* $name:ident: $light:expr, $dark:expr;)*) => {
        $(
            $(#[$doc])*
            pub fn $name() -> Hsla {
                if dark() { $dark } else { $light }
            }
        )*
    };
}

palette! {
    /// Window background.
    paper: hsla(0.1250, 0.3333, 0.9765, 1.00), hsla(0.1111, 0.0600, 0.0950, 1.00);

    /// Cards and panels.
    surface: hsla(0.0000, 0.0000, 1.0000, 1.00), hsla(0.1111, 0.0500, 0.1250, 1.00);

    /// Input wells and inactive chips.
    sunken: hsla(0.1333, 0.1852, 0.9471, 1.00), hsla(0.1111, 0.0600, 0.0700, 1.00);

    /// Hairlines.
    border: hsla(0.1111, 0.1579, 0.8882, 1.00), hsla(0.1111, 0.0500, 0.2050, 1.00);

    /// Focused field outline.
    border_focus: hsla(0.6667, 0.0833, 0.0941, 1.00), hsla(0.1111, 0.1500, 0.9200, 1.00);

    /// Primary text, and the fill of primary buttons.
    ink: hsla(0.6667, 0.0833, 0.0941, 1.00), hsla(0.1111, 0.1500, 0.9200, 1.00);

    /// Body text.
    text: hsla(0.6667, 0.0617, 0.1588, 1.00), hsla(0.1111, 0.0800, 0.8300, 1.00);

    /// Secondary text.
    muted: hsla(0.1061, 0.0435, 0.5039, 1.00), hsla(0.1061, 0.0435, 0.5800, 1.00);

    /// Placeholders and disabled text.
    faint: hsla(0.1111, 0.0750, 0.6863, 1.00), hsla(0.1111, 0.0500, 0.4200, 1.00);

    /// Text caret.
    caret: hsla(0.6667, 0.0833, 0.0941, 1.00), hsla(0.1111, 0.1500, 0.9200, 1.00);

    /// Text selection.
    selection: hsla(0.6667, 0.0833, 0.0941, 0.14), hsla(0.1111, 0.1500, 0.9200, 0.20);

    /// Errors.
    danger: hsla(0.0092, 0.7005, 0.4059, 1.00), hsla(0.0092, 0.7500, 0.6200, 1.00);

    /// Text on an ink-filled surface.
    on_ink: hsla(0.1250, 0.3333, 0.9765, 1.00), hsla(0.1111, 0.0600, 0.0950, 1.00);
}

/// The page behind the molecule, in 0..1 RGB for the rasterizer. Every style
/// sits on the same flat paper, so switching between them changes the drawing
/// and nothing else.
pub fn view() -> [f32; 3] {
    if dark() {
        [0.110, 0.108, 0.104]
    } else {
        [0.992, 0.992, 0.988]
    }
}

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
