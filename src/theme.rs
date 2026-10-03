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

/// Monospace family, used for sequences and identifiers.
pub const MONO: &str = "Menlo";

/// Display serif, used only for the name of the structure on screen.
pub const SERIF: &str = "Charter";
