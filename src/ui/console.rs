//! The console: commands typed, or read from a script, that do what the
//! buttons and menus do, so a figure can be made again the same way.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context as _, Result};
use gpui::Context;

use super::{residue_pick, spread, Layer, Target, Turn, Viewer};
use crate::formats::write;
use crate::measure::{self, Measurement};
use crate::render::{self, Coloring, LayerStyle, Palette, Repr, Style};
use crate::select::Selection;

/// What `help` lists: each command and what it takes.
pub const HELP: &[(&str, &str)] = &[
    ("open <id, name or path>", "open a structure, as the search field does"),
    ("overlay <id or path>", "lay another structure over this one"),
    ("show cartoon|tube|sticks|licorice|spheres|surface", "draw it so"),
    ("color chain|rainbow|plddt|bfactor|secondary|hydrophobic|charge|element", "colour it so"),
    ("palette <name>", "Spectrum, Viridis, Magma and the rest"),
    ("style illustrative|flat|soft|studio", "the look of the drawing"),
    ("sticks|spheres|hide <selection>", "draw a selection its own way"),
    ("clear", "take every selection away"),
    ("distance|angle|dihedral <sel>, <sel>, …", "measure between the first atom of each"),
    ("focus <selection>", "fly to it"),
    ("site <ligand>", "look into a ligand's binding site"),
    ("assembly <n>|deposited", "show an assembly, or what was deposited"),
    ("turn x|y <degrees>, zoom <factor>, reset", "move the camera"),
    ("spin, rock, stop", "turn by itself, or stop"),
    ("perspective on|off, slab on|off", "the lens and the slab"),
    ("water on|off, hydrogens on|off", "show them or not"),
    ("frame <n>, play", "step a trajectory"),
    ("save <file.png|.pdb|.cif|.profold>", "write it out, no questions asked"),
    ("run <script>", "the commands in a file, one a line"),
    ("undo, redo, help", ""),
];

fn on_off(arg: &str) -> Result<bool> {
    match arg.trim() {
        "on" | "yes" | "true" | "1" | "" => Ok(true),
        "off" | "no" | "false" | "0" => Ok(false),
        other => bail!("\"{other}\" is neither on nor off"),
    }
}

impl Viewer {
    /// Run one command. What it says back, or what went wrong.
    pub(super) fn run_command(&mut self, line: &str, cx: &mut Context<Self>) -> Result<String> {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            return Ok(String::new());
        }
        let (verb, rest) = line.split_once(char::is_whitespace).unwrap_or((line, ""));
        let rest = rest.trim();
        let verb = verb.to_ascii_lowercase();
        let need = |what: &str| -> Result<()> {
            if rest.is_empty() {
                bail!("{verb} wants {what}");
            }
            Ok(())
        };
        let shown = || -> Result<Arc<crate::protein::Structure>> {
            self.structure.clone().ok_or_else(|| anyhow!("nothing is open"))
        };
        Ok(match verb.as_str() {
            "help" | "?" => HELP.iter().map(|(c, _)| *c).collect::<Vec<_>>().join(" · "),
            "open" | "load" | "fetch" => {
                need("something to open")?;
                let path = PathBuf::from(rest);
                if path.is_file() {
                    self.open_path(path, cx);
                } else {
                    self.open(rest, cx);
                }
                format!("opening {rest}")
            }
            "overlay" => {
                need("something to lay over")?;
                shown()?;
                self.overlay_next = true;
                let path = PathBuf::from(rest);
                if path.is_file() {
                    self.open_path(path, cx);
                } else {
                    self.open(rest, cx);
                }
                format!("laying {rest} over")
            }
            "show" | "as" => {
                let repr = match rest.to_ascii_lowercase().as_str() {
                    "cartoon" | "ribbon" => Repr::Cartoon,
                    "tube" | "trace" => Repr::Tube,
                    "sticks" | "ballandstick" | "ball&stick" | "balls" => Repr::BallAndStick,
                    "licorice" => Repr::Licorice,
                    "spheres" | "spacefill" | "cpk" => Repr::Spacefill,
                    "surface" => Repr::Surface,
                    other => bail!("\"{other}\" is not a way to draw it"),
                };
                let s = shown()?;
                if repr.needs_polymer() && !s.has_polymer() {
                    bail!("{} needs a backbone", repr.label());
                }
                self.set_repr(repr, cx);
                format!("drawn as {}", repr.label())
            }
            "color" | "colour" => {
                let want = rest.to_ascii_lowercase();
                let coloring = match want.as_str() {
                    "chain" | "chains" | "molecule" => Coloring::Chain,
                    "rainbow" => Coloring::Rainbow,
                    "plddt" | "confidence" => Coloring::Confidence,
                    "bfactor" | "b-factor" | "b" => Coloring::BFactor,
                    "secondary" | "ss" | "structure" => Coloring::Structure,
                    "hydrophobic" | "hydrophobicity" => Coloring::Hydrophobic,
                    "charge" | "electrostatic" | "electrostatics" => Coloring::Charge,
                    "element" | "elements" | "atom" => Coloring::Element,
                    _ => {
                        // A palette's name colours by it too.
                        return self.run_command(&format!("palette {rest}"), cx);
                    }
                };
                shown()?;
                self.coloring = coloring;
                self.rebuild(cx);
                format!("coloured by {}", coloring.label())
            }
            "palette" => {
                let p = Palette::ALL
                    .into_iter()
                    .find(|p| p.label().eq_ignore_ascii_case(rest))
                    .ok_or_else(|| anyhow!("\"{rest}\" is not a palette"))?;
                self.palette = p;
                self.rebuild(cx);
                format!("in {}", p.label())
            }
            "style" | "look" => {
                let st = Style::ALL
                    .into_iter()
                    .find(|st| st.label().eq_ignore_ascii_case(rest))
                    .ok_or_else(|| anyhow!("\"{rest}\" is not a style"))?;
                self.set_style(st, cx);
                format!("in the {} style", st.label())
            }
            "sticks" | "spheres" | "hide" | "select" => {
                need("a selection")?;
                shown()?;
                let selection = Selection::parse(rest)?;
                let style = match verb.as_str() {
                    "spheres" => LayerStyle::Spheres,
                    "hide" => LayerStyle::Hidden,
                    _ => LayerStyle::Sticks,
                };
                let n = selection.atoms(self.structure.as_ref().unwrap()).iter().filter(|a| **a).count();
                self.layers.push(Layer { selection, atoms: Arc::new(Vec::new()), style, color: None });
                self.rebuild(cx);
                format!("{n} atoms as {}", style.label().to_lowercase())
            }
            "clear" => {
                self.layers.clear();
                self.measurements.clear();
                self.rebuild(cx);
                "selections and measurements cleared".into()
            }
            "distance" | "angle" | "dihedral" => {
                let kind = match verb.as_str() {
                    "distance" => measure::Kind::Distance,
                    "angle" => measure::Kind::Angle,
                    _ => measure::Kind::Dihedral,
                };
                let s = shown()?;
                let parts: Vec<&str> = rest.split([',', ';']).map(str::trim).filter(|p| !p.is_empty()).collect();
                if parts.len() != kind.atoms() {
                    bail!("{verb} wants {} selections, a comma between each", kind.atoms());
                }
                let atoms: Vec<usize> = parts
                    .iter()
                    .map(|p| {
                        let mask = Selection::parse(p)?.atoms(&s);
                        mask.iter().position(|a| *a).ok_or_else(|| anyhow!("\"{p}\" takes no atom"))
                    })
                    .collect::<Result<_>>()?;
                let m = Measurement { kind, which: 0, atoms };
                let value = m.value(&s).map(|v| kind.format(v)).unwrap_or_default();
                self.measurements.push(m);
                self.compose(cx);
                format!("{verb} {value}")
            }
            "focus" | "zoomto" | "center" | "centre" => {
                need("a selection")?;
                let s = shown()?;
                let mask = Selection::parse(rest)?.atoms(&s);
                let points: Vec<[f32; 3]> = (0..s.atoms.len()).filter(|&i| mask[i]).map(|i| s.atoms[i].pos).collect();
                if points.is_empty() {
                    bail!("\"{rest}\" takes no atom");
                }
                let (centre, reach) = spread(points.iter().copied());
                self.fly_to(centre, reach * 2.0 + 8.0, None, cx);
                format!("{} atoms in view", points.len())
            }
            "site" => {
                need("a ligand's name")?;
                let s = shown()?;
                if s.copies(&rest.to_uppercase()).is_empty() {
                    bail!("there is no ligand called {rest}");
                }
                self.next_site(&rest.to_uppercase(), cx);
                format!("the site of {}", rest.to_uppercase())
            }
            "assembly" => {
                let k = match rest {
                    "" | "deposited" | "none" | "0" => None,
                    n => Some(n.parse::<usize>().context("assembly wants a number")?.checked_sub(1).ok_or_else(|| anyhow!("assemblies count from 1"))?),
                };
                let base = self.deposited.clone().or_else(|| self.structure.clone()).ok_or_else(|| anyhow!("nothing is open"))?;
                if let Some(k) = k {
                    if k >= base.assemblies.len() {
                        bail!("there are {} assemblies", base.assemblies.len());
                    }
                }
                self.show_assembly(k, cx);
                match k {
                    Some(k) => format!("assembly {}", k + 1),
                    None => "as deposited".into(),
                }
            }
            "turn" | "rotate" => {
                let mut it = rest.split_whitespace();
                let (axis, deg) = (it.next().unwrap_or(""), it.next().and_then(|d| d.parse::<f32>().ok()));
                let Some(deg) = deg else { bail!("turn wants an axis and degrees, as in \"turn y 90\"") };
                let r = deg.to_radians();
                match axis {
                    "y" => self.camera.rotate(r, 0.0),
                    "x" => self.camera.rotate(0.0, r),
                    _ => bail!("turn about x or y"),
                }
                cx.notify();
                format!("turned {deg}° about {axis}")
            }
            "zoom" => {
                let f: f32 = rest.parse().context("zoom wants a factor")?;
                self.camera.zoom = (self.camera.zoom * f).clamp(0.12, 24.0);
                cx.notify();
                format!("zoom {:.2}×", self.camera.zoom)
            }
            "reset" => {
                self.reset_view(cx);
                "view reset".into()
            }
            "spin" => {
                self.toggle_turn(Turn::Spin, cx);
                "spinning".into()
            }
            "rock" => {
                self.toggle_turn(Turn::Rock, cx);
                "rocking".into()
            }
            "stop" => {
                self.halt();
                self.stop_playing();
                cx.notify();
                "stopped".into()
            }
            "perspective" => {
                self.camera.perspective = on_off(rest)?;
                cx.notify();
                format!("perspective {}", if self.camera.perspective { "on" } else { "off" })
            }
            "slab" => {
                if on_off(rest)? != self.camera.slab.is_some() {
                    self.toggle_slab(cx);
                }
                format!("slab {}", if self.camera.slab.is_some() { "on" } else { "off" })
            }
            "water" => {
                self.appearance.hide_solvent = !on_off(rest)?;
                self.rebuild(cx);
                format!("water {}", if self.appearance.hide_solvent { "hidden" } else { "shown" })
            }
            "hydrogens" | "hydrogen" => {
                let on = on_off(rest)?;
                let h = render::Species::Element(1);
                let look = self.appearance.get(&h);
                self.appearance.set(h, render::Look { size: if on { 1.0 } else { 0.0 }, ..look });
                self.rebuild(cx);
                format!("hydrogens {}", if on { "shown" } else { "hidden" })
            }
            "frame" => {
                let n: usize = rest.parse().context("frame wants a number")?;
                let count = self.frame_count();
                if n == 0 || n > count {
                    bail!("there are {count} frames");
                }
                self.set_frame(n - 1, cx);
                format!("frame {n} of {count}")
            }
            "play" => {
                if !self.playing {
                    self.toggle_play(cx);
                }
                "playing".into()
            }
            "pause" => {
                self.stop_playing();
                cx.notify();
                "paused".into()
            }
            "save" | "write" => {
                need("a file to write")?;
                let path = Path::new(rest);
                let ext = path.extension().map(|e| e.to_string_lossy().to_ascii_lowercase()).unwrap_or_default();
                let s = shown()?;
                match ext.as_str() {
                    "pdb" => std::fs::write(path, write::pdb(&s, |_| true))?,
                    "cif" | "mmcif" => std::fs::write(path, write::mmcif(&s, |_| true))?,
                    "profold" => std::fs::write(path, self.session_json().unwrap_or_default())?,
                    "png" => {
                        let (w, h) = self.cache.borrow().size;
                        let (w, h) = (w.max(640), h.max(480));
                        let frame = render::render(&self.scene, &self.camera, w, h, crate::theme::view(), self.style);
                        let mut rgba = frame.bgra;
                        for p in rgba.chunks_exact_mut(4) {
                            p.swap(0, 2);
                        }
                        image::save_buffer(path, &rgba, w as u32, h as u32, image::ColorType::Rgba8)?;
                    }
                    other => bail!("save writes .png, .pdb, .cif or .profold, not .{other}"),
                }
                format!("saved {}", path.display())
            }
            "run" | "source" | "@" => {
                need("a script")?;
                let text = std::fs::read_to_string(rest).with_context(|| format!("could not read {rest}"))?;
                let lines: Vec<String> = text.lines().map(str::to_string).collect();
                self.run_script(lines, cx);
                format!("ran {rest}")
            }
            "undo" => {
                self.undo(false, cx);
                "undone".into()
            }
            "redo" => {
                self.undo(true, cx);
                "redone".into()
            }
            "pick" => {
                let s = shown()?;
                let mask = Selection::parse(rest)?.atoms(&s);
                let atom = mask.iter().position(|a| *a).ok_or_else(|| anyhow!("\"{rest}\" takes no atom"))?;
                let target = residue_pick(&s, s.atoms[atom].residue as usize).map(|picked| Target { which: 0, picked });
                self.select(target, cx);
                "picked".into()
            }
            other => bail!("\"{other}\" is not a command; help lists them"),
        })
    }

    /// Run lines one after another. Opening something waits: the lines after
    /// it run once it is open.
    pub(super) fn run_script(&mut self, lines: Vec<String>, cx: &mut Context<Self>) {
        let mut lines = lines.into_iter();
        while let Some(line) = lines.next() {
            let result = self.run_command(&line, cx);
            self.say(&line, result);
            // Whatever comes after an open waits for what it opens.
            if self.busy.is_some() {
                self.script = lines.collect();
                return;
            }
        }
    }

    /// Put a command and what came of it in the console's log.
    pub(super) fn say(&mut self, line: &str, result: Result<String>) {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            return;
        }
        let (ok, text) = match result {
            Ok(t) => (true, t),
            Err(e) => (false, format!("{e:#}")),
        };
        self.console_log.push((ok, format!("{line}  →  {text}")));
        let len = self.console_log.len();
        if len > 50 {
            self.console_log.drain(..len - 50);
        }
    }
}
