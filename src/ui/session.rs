//! A session: the structure on screen and how it is shown, in one file that
//! opens again as it was left.
//!
//! The structure goes in whole, written as mmCIF, so a session still opens
//! once the file it came from has moved, and an assembly or a download
//! comes back without the network. Around it: the camera, the style and
//! colours, each chain's and ligand's own colour, the measurements, what is
//! picked out, and the binding site looked into. A trajectory is kept at
//! the frame shown.

use anyhow::{anyhow, bail, Context as _, Result};
use gpui::Context;
use serde_json::{json, Value};

use super::{Target, Viewer};
use crate::formats::{self, write, Format};
use crate::measure::{self, Measurement};
use crate::protein::{LigandPaint, Source, Structure};
use crate::render::{Camera, Coloring, Palette, Picked, Repr, Style};

/// The extension a session goes by.
pub const EXTENSION: &str = "profold";

const VERSION: u64 = 1;

fn rgb(c: [f32; 3]) -> Value {
    json!(c)
}

fn read_rgb(v: &Value) -> Option<[f32; 3]> {
    let a = v.as_array()?;
    Some([a.first()?.as_f64()? as f32, a.get(1)?.as_f64()? as f32, a.get(2)?.as_f64()? as f32])
}

fn floats<const N: usize>(v: &Value) -> Option<[f32; N]> {
    let a = v.as_array()?;
    let mut out = [0.0; N];
    for (i, x) in out.iter_mut().enumerate() {
        *x = a.get(i)?.as_f64()? as f32;
    }
    Some(out)
}

impl Viewer {
    /// The session as JSON, or `None` with nothing on screen.
    pub(super) fn session_json(&self) -> Option<String> {
        let s = self.structure.as_ref()?;
        let c = &self.camera;
        let o = c.orientation;
        let chains: Vec<Value> = s
            .chains
            .iter()
            .filter_map(|ch| Some(json!([ch.id, rgb(ch.color?)])))
            .collect();
        let ligands: Vec<Value> = s
            .ligand_paint
            .iter()
            .map(|(name, paint)| match paint {
                LigandPaint::Solid(c) => json!([name, "solid", rgb(*c)]),
                LigandPaint::Carbons(c) => json!([name, "carbons", rgb(*c)]),
            })
            .collect();
        let measurements: Vec<Value> = self
            .measurements
            .iter()
            .filter(|m| m.which == 0)
            .map(|m| json!({ "kind": m.kind.label(), "atoms": m.atoms }))
            .collect();
        let selected = self
            .selected
            .filter(|t| t.which == 0)
            .map(|t| json!({ "atom": t.picked.atom, "residue": t.picked.residue }));
        let value = json!({
            "profold": VERSION,
            "id": s.id,
            "source": s.source.label(),
            "origin": self.origin,
            "structure": write::mmcif(s, |_| true),
            "camera": {
                "orientation": [o[0][0], o[0][1], o[0][2], o[1][0], o[1][1], o[1][2], o[2][0], o[2][1], o[2][2]],
                "zoom": c.zoom,
                "pan": c.pan,
                "center": c.center,
                "radius": c.radius,
                "perspective": c.perspective,
                "slab": c.slab,
            },
            "repr": self.repr.label(),
            "coloring": self.coloring.label(),
            "style": self.style.label(),
            "palette": self.palette.label(),
            "standard_colors": self.standard_colors,
            "show_cell": self.show_cell,
            "chain_colors": chains,
            "ligand_paint": ligands,
            "measurements": measurements,
            "selected": selected,
            "site": self.site.as_ref().map(|site| site.ligand),
        });
        serde_json::to_string_pretty(&value).ok()
    }

    /// Open a session saved by `session_json`.
    pub(super) fn open_session(&mut self, text: &str, cx: &mut Context<Self>) -> Result<()> {
        let v: Value = serde_json::from_str(text).context("reading the session")?;
        if v.get("profold").and_then(Value::as_u64).is_none_or(|n| n > VERSION) {
            bail!("this session was saved by a newer profold");
        }
        let cif = v["structure"].as_str().ok_or_else(|| anyhow!("the session holds no structure"))?;
        let id = v["id"].as_str().unwrap_or("session").to_string();
        let mut s: Structure = formats::parse(cif, Format::Cif, id, Source::File)?;
        // Where it came from, so it still says so.
        if let Some(label) = v["source"].as_str() {
            for source in [Source::Pdb, Source::AlphaFold, Source::PubChem, Source::SimpleFold, Source::File] {
                if source.label() == label {
                    s.source = source;
                }
            }
        }
        for c in v["chain_colors"].as_array().into_iter().flatten() {
            if let (Some(id), Some(color)) = (c[0].as_str(), read_rgb(&c[1])) {
                if let Some(ch) = s.chains.iter_mut().find(|ch| ch.id == id) {
                    ch.color = Some(color);
                }
            }
        }
        for l in v["ligand_paint"].as_array().into_iter().flatten() {
            if let (Some(name), Some(kind), Some(color)) = (l[0].as_str(), l[1].as_str(), read_rgb(&l[2])) {
                let paint = if kind == "carbons" { LigandPaint::Carbons(color) } else { LigandPaint::Solid(color) };
                s.ligand_paint.insert(name.to_string(), paint);
            }
        }
        let atoms = s.atoms.len();

        self.show(s, cx);
        let pick = |list: &[&'static str], v: &Value| v.as_str().and_then(|t| list.iter().position(|l| *l == t));
        if let Some(i) = pick(&Repr::ALL.map(Repr::label), &v["repr"]) {
            self.repr = Repr::ALL[i];
        }
        if let Some(i) = pick(&Coloring::ALL.map(Coloring::label), &v["coloring"]) {
            self.coloring = Coloring::ALL[i];
        }
        if let Some(i) = pick(&Style::ALL.map(Style::label), &v["style"]) {
            self.style = Style::ALL[i];
        }
        if let Some(i) = pick(&Palette::ALL.map(Palette::label), &v["palette"]) {
            self.palette = Palette::ALL[i];
        }
        if let Some(b) = v["standard_colors"].as_bool() {
            self.standard_colors = b;
        }
        if let Some(b) = v["show_cell"].as_bool() {
            self.show_cell = b;
        }
        let cam = &v["camera"];
        if let (Some(o), Some(pan), Some(center)) = (floats::<9>(&cam["orientation"]), floats::<2>(&cam["pan"]), floats::<3>(&cam["center"])) {
            self.camera = Camera {
                orientation: [[o[0], o[1], o[2]], [o[3], o[4], o[5]], [o[6], o[7], o[8]]],
                zoom: cam["zoom"].as_f64().unwrap_or(1.0) as f32,
                pan,
                center,
                radius: cam["radius"].as_f64().unwrap_or(self.camera.radius as f64) as f32,
                perspective: cam["perspective"].as_bool().unwrap_or(false),
                slab: floats::<2>(&cam["slab"]),
            };
        }
        let in_range = |i: &Value| i.as_u64().map(|i| i as usize).filter(|&i| i < atoms);
        self.measurements = v["measurements"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(|m| {
                let kind = measure::Kind::ALL.into_iter().find(|k| Some(k.label()) == m["kind"].as_str())?;
                let atoms: Vec<usize> = m["atoms"].as_array()?.iter().map(in_range).collect::<Option<_>>()?;
                (atoms.len() == kind.atoms()).then_some(Measurement { kind, which: 0, atoms })
            })
            .collect();
        self.selected = in_range(&v["selected"]["atom"]).map(|atom| Target {
            which: 0,
            picked: Picked { atom, residue: v["selected"]["residue"].as_bool().unwrap_or(false), depth: 0.0 },
        });
        self.rebuild(cx);
        if let Some(ligand) = v["site"].as_u64().map(|r| r as usize) {
            let camera = self.camera;
            if self.structure.as_ref().is_some_and(|s| ligand < s.residues.len()) {
                self.show_site(ligand, cx);
                // The site comes back as it was seen, without flying there.
                self.fly_task = None;
                self.camera = camera;
            }
        }
        cx.notify();
        Ok(())
    }
}
