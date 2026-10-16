//! Viewer preferences that outlive a session: momentum, style, palette, and
//! the size, opacity and colour set for each element or atom type.
//!
//! Kept as JSON in `~/.profold/settings.json`, beside SimpleFold's
//! environment. A file that is missing or unreadable means the defaults; a
//! setting the app no longer knows is ignored rather than refused.

use std::path::PathBuf;

use serde_json::{json, Map, Value};

use crate::elements::Element;
use crate::render::{Appearance, Look, Palette, Species, Style};

#[derive(Clone)]
pub struct Settings {
    pub momentum: bool,
    pub style: Style,
    pub palette: Palette,
    /// Elements in their standard colours rather than the palette's. Only
    /// elements: the palette still colours chains and ramps.
    pub standard_colors: bool,
    pub appearance: Appearance,
    /// Drawn dark rather than light.
    pub dark: bool,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            momentum: true,
            style: Style::Illustrative,
            palette: Palette::Spectrum,
            standard_colors: true,
            appearance: Appearance::default(),
            dark: false,
        }
    }
}

fn path() -> Option<PathBuf> {
    file("settings.json")
}

/// A file of profold's own under `~/.profold`.
pub fn file(name: &str) -> Option<PathBuf> {
    // Tests run the real app, and must not read or overwrite the real file.
    if cfg!(test) {
        return None;
    }
    Some(dirs::home_dir()?.join(".profold").join(name))
}

impl Settings {
    pub fn load() -> Settings {
        path()
            .and_then(|p| std::fs::read_to_string(p).ok())
            .and_then(|text| serde_json::from_str::<Value>(&text).ok())
            .map(|v| Settings::from_json(&v))
            .unwrap_or_default()
    }

    /// Write the settings out, quietly: losing a preference is not worth
    /// interrupting anyone over.
    pub fn save(json: &str) {
        let Some(p) = path() else { return };
        if let Some(dir) = p.parent() {
            let _ = std::fs::create_dir_all(dir);
        }
        let _ = std::fs::write(p, json);
    }

    pub fn to_json(&self) -> String {
        let looks: Vec<Value> = {
            let mut entries: Vec<(&Species, &Look)> = self.appearance.entries().collect();
            // A stable order, so an unchanged setting writes an unchanged file.
            entries.sort_by_key(|(s, _)| match s {
                Species::Element(z) => (0, *z as usize, String::new()),
                Species::Kind(k) => (1, 0, k.clone()),
            });
            entries
                .into_iter()
                .map(|(species, look)| {
                    let mut m = Map::new();
                    match species {
                        Species::Element(z) => {
                            m.insert("element".into(), json!(Element(*z).symbol()));
                        }
                        Species::Kind(k) => {
                            m.insert("type".into(), json!(k));
                        }
                    }
                    m.insert("size".into(), json!(look.size));
                    m.insert("opacity".into(), json!(look.alpha));
                    if let Some(c) = look.color {
                        m.insert("color".into(), json!(hex(c)));
                    }
                    Value::Object(m)
                })
                .collect()
        };
        let mut v = json!({
            "momentum": self.momentum,
            "style": self.style.label(),
            "palette": self.palette.label(),
            "element_colours": if self.standard_colors { "standard" } else { "palette" },
            "atoms": looks,
        });
        if self.dark {
            v["theme"] = json!("dark");
        }
        serde_json::to_string_pretty(&v).unwrap_or_default()
    }

    fn from_json(v: &Value) -> Settings {
        let mut s = Settings::default();
        if let Some(m) = v.get("momentum").and_then(Value::as_bool) {
            s.momentum = m;
        }
        if let Some(name) = v.get("style").and_then(Value::as_str) {
            if let Some(st) = Style::ALL.iter().find(|st| st.label() == name) {
                s.style = *st;
            }
        }
        s.dark = v.get("theme").and_then(Value::as_str) == Some("dark");
        // The old "fixed_colours" is not read: it was off unless chosen, and
        // standard colours are now the default.
        if let Some(c) = v.get("element_colours").and_then(Value::as_str) {
            s.standard_colors = c != "palette";
        }
        if let Some(name) = v.get("palette").and_then(Value::as_str) {
            if let Some(p) = Palette::ALL.iter().find(|p| p.label() == name) {
                s.palette = *p;
            }
        }
        for entry in v.get("atoms").and_then(Value::as_array).into_iter().flatten() {
            let species = if let Some(sym) = entry.get("element").and_then(Value::as_str) {
                let e = Element::from_symbol(sym);
                if !e.is_known() {
                    continue;
                }
                Species::Element(e.0)
            } else if let Some(k) = entry.get("type").and_then(Value::as_str) {
                Species::Kind(k.to_string())
            } else {
                continue;
            };
            let number = |key: &str, default: f32| {
                entry
                    .get(key)
                    .and_then(Value::as_f64)
                    .map(|x| x as f32)
                    .unwrap_or(default)
            };
            let look = Look {
                size: number("size", 1.0).clamp(0.0, MAX_SIZE),
                alpha: number("opacity", 1.0).clamp(0.0, 1.0),
                color: entry.get("color").and_then(Value::as_str).and_then(parse_hex),
            };
            s.appearance.set(species, look);
        }
        s
    }
}

// ---------------------------------------------------------------------------
// Values typed into fields
// ---------------------------------------------------------------------------

/// Largest size multiple that can be typed; the slider stops sooner.
pub const MAX_SIZE: f32 = 10.0;

fn plain_number(text: &str, suffixes: &[char]) -> Option<f32> {
    let t = text.trim().trim_end_matches(suffixes).trim();
    let v: f32 = t.parse().ok()?;
    v.is_finite().then_some(v)
}

/// A size multiple: "1", "1.5", "2x", "0.5×".
pub fn parse_size(text: &str) -> Option<f32> {
    plain_number(text, &['x', 'X', '×']).map(|v| v.clamp(0.0, MAX_SIZE))
}

/// An opacity, as a fraction or a percentage: "0.1", "10%" and "10" are all
/// one tenth. A bare number above one can only be a percentage.
pub fn parse_opacity(text: &str) -> Option<f32> {
    let t = text.trim();
    let v = if let Some(p) = t.strip_suffix('%') {
        p.trim().parse::<f32>().ok()? / 100.0
    } else {
        let v: f32 = t.parse().ok()?;
        if v > 1.0 {
            v / 100.0
        } else {
            v
        }
    };
    v.is_finite().then(|| v.clamp(0.0, 1.0))
}

/// An angle in degrees: "30", "-45°".
pub fn parse_degrees(text: &str) -> Option<f32> {
    plain_number(text, &['°'])
}

/// A zoom multiple: "2", "1.5x".
pub fn parse_zoom(text: &str) -> Option<f32> {
    plain_number(text, &['x', 'X', '×']).filter(|v| *v > 0.0)
}

/// "#ff8800", "ff8800" or the short "#f80".
pub fn parse_hex(text: &str) -> Option<[f32; 3]> {
    let t = text.trim().trim_start_matches('#');
    let full: String = match t.len() {
        3 => t.chars().flat_map(|c| [c, c]).collect(),
        6 => t.to_string(),
        _ => return None,
    };
    let byte = |i: usize| u8::from_str_radix(&full[i..i + 2], 16).ok();
    Some([
        byte(0)? as f32 / 255.0,
        byte(2)? as f32 / 255.0,
        byte(4)? as f32 / 255.0,
    ])
}

/// "255, 128, 0", "255 128 0" or "rgb(255, 128, 0)". Three numbers all at
/// most one, with a decimal point among them, are read as fractions.
pub fn parse_rgb(text: &str) -> Option<[f32; 3]> {
    let t = text.trim();
    let t = t
        .strip_prefix("rgb(")
        .or_else(|| t.strip_prefix("RGB("))
        .and_then(|r| r.strip_suffix(')'))
        .unwrap_or(t);
    let parts: Vec<&str> = t
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|p| !p.is_empty())
        .collect();
    if parts.len() != 3 {
        return None;
    }
    let v: Vec<f32> = parts.iter().map(|p| p.parse().ok()).collect::<Option<_>>()?;
    let fractions = v.iter().all(|x| *x <= 1.0) && parts.iter().any(|p| p.contains('.'));
    let scale = if fractions { 1.0 } else { 255.0 };
    Some([
        (v[0] / scale).clamp(0.0, 1.0),
        (v[1] / scale).clamp(0.0, 1.0),
        (v[2] / scale).clamp(0.0, 1.0),
    ])
}

pub fn hex(c: [f32; 3]) -> String {
    let b = |x: f32| (x.clamp(0.0, 1.0) * 255.0).round() as u8;
    format!("#{:02x}{:02x}{:02x}", b(c[0]), b(c[1]), b(c[2]))
}

pub fn rgb(c: [f32; 3]) -> String {
    let b = |x: f32| (x.clamp(0.0, 1.0) * 255.0).round() as u8;
    format!("{}, {}, {}", b(c[0]), b(c[1]), b(c[2]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn typed_values() {
        assert_eq!(parse_size("1"), Some(1.0));
        assert_eq!(parse_size("2.5×"), Some(2.5));
        assert_eq!(parse_size("abc"), None);
        assert_eq!(parse_opacity("0.1"), Some(0.1));
        assert_eq!(parse_opacity("10%"), Some(0.1));
        assert_eq!(parse_opacity("10"), Some(0.1));
        assert_eq!(parse_opacity("1"), Some(1.0));
        assert_eq!(parse_degrees("-45°"), Some(-45.0));
        assert_eq!(parse_zoom("0"), None);
        assert_eq!(parse_hex("#f80"), Some([1.0, 0x88 as f32 / 255.0, 0.0]));
        assert_eq!(parse_hex("ff8800").map(hex), Some("#ff8800".into()));
        assert_eq!(parse_rgb("rgb(255, 128, 0)").map(hex), Some("#ff8000".into()));
        assert_eq!(parse_rgb("255 128 0").map(hex), Some("#ff8000".into()));
        assert_eq!(parse_rgb("1.0, 0.5, 0").map(hex), Some("#ff8000".into()));
        assert_eq!(parse_rgb("1, 2"), None);
    }

    #[test]
    fn round_trip() {
        let mut s = Settings {
            momentum: false,
            style: Style::Soft,
            palette: Palette::Viridis,
            standard_colors: false,
            ..Settings::default()
        };
        s.appearance.set(
            Species::Element(1),
            Look {
                size: 0.0,
                alpha: 1.0,
                color: None,
            },
        );
        s.appearance.set(
            Species::Kind("2".into()),
            Look {
                size: 1.5,
                alpha: 0.25,
                color: Some([1.0, 0.0, 0.0]),
            },
        );
        let text = s.to_json();
        let back = Settings::from_json(&serde_json::from_str(&text).unwrap());
        assert_eq!(back.to_json(), text);
        assert!(!back.momentum);
        assert_eq!(back.style, Style::Soft);
        assert!(!back.standard_colors);
        assert_eq!(back.appearance.get(&Species::Kind("2".into())).alpha, 0.25);
    }

    #[test]
    fn every_style_and_palette_is_remembered() {
        for style in Style::ALL {
            for palette in Palette::ALL {
                let s = Settings {
                    style,
                    palette,
                    ..Settings::default()
                };
                let back = Settings::from_json(&serde_json::from_str(&s.to_json()).unwrap());
                assert_eq!((back.style, back.palette), (style, palette));
            }
        }
    }
}
