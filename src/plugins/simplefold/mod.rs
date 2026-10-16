//! Structure prediction with Apple's SimpleFold.
//!
//! The viewer works without this; it is compiled in by the `simplefold`
//! feature, and appears only in a fold tab, which is opened from the + menu.
//! It contributes three things: jobs the tab runs (fold a sequence, fold an
//! entry's own chains, fold and compare against the deposited structure), the
//! settings row for them, and the Compare button. Everything else — running
//! jobs, showing structures, superimposed pairs — is the viewer's.

mod download;
mod engine;

use std::sync::mpsc::Sender;
use std::sync::Arc;

use gpui::{deferred, div, prelude::*, px, Context, SharedString};

use crate::align;
use crate::fetch::{self, Retrieval};
use crate::formats::parse_any;
use crate::job::{Event, Outcome};
use crate::protein::{Atom, Chain, Source, Structure};
use crate::theme;
use crate::ui::app::{label, pill, Viewer};
use engine::{model_label, Engine, MODELS};

/// Work for the app to run in the background, and what to say first.
pub struct Job {
    pub stage: &'static str,
    pub notes: Vec<String>,
    pub work: Box<dyn FnOnce(&Sender<Event>) -> anyhow::Result<Outcome> + Send>,
}

pub struct SimpleFold {
    engine: Arc<Engine>,
    model: usize,
    steps: u32,
    plddt: bool,
    menu_open: bool,
}

/// An example that shows off folding, for the list on the empty screen.
pub const EXAMPLE: (&str, &str) = (
    "GSHMKAEELKKKAEELLKKAEELAKKLAEELKKLAEELKKKAEELLKKLA",
    "a designed helix, folded here",
);

impl SimpleFold {
    pub fn new() -> Self {
        SimpleFold {
            engine: Arc::new(Engine::new()),
            model: 0,
            steps: 200,
            plddt: false,
            menu_open: false,
        }
    }

    /// What every job says before it starts, about downloads it will make.
    fn preamble(&self) -> Vec<String> {
        let mut notes = Vec::new();
        if !self.engine.is_installed() {
            notes.push("First run: installing SimpleFold. This takes a while.".into());
        }
        if self.plddt && !self.engine.has_plddt_weights() {
            notes.push("pLDDT scoring needs SimpleFold's 1.6B checkpoint (~6.5 GB), once.".into());
        }
        notes
    }

    pub fn fold(&self, sequence: String) -> Job {
        let engine = self.engine.clone();
        let model = MODELS[self.model].to_string();
        let label = model_label(MODELS[self.model]).to_string();
        let steps = self.steps;
        let plddt = self.plddt;
        let installing = !engine.is_installed();
        Job {
            stage: if installing {
                "Setting up the folding engine"
            } else {
                "Folding"
            },
            notes: self.preamble(),
            work: Box::new(move |tx| {
                if installing {
                    engine.install(tx)?;
                }
                let path = engine.fold("query", &sequence, &model, steps, plddt, tx)?;
                let _ = tx.send(Event::Stage("Reading structure".into()));
                let text = std::fs::read_to_string(&path)?;
                // The file is named after the FASTA record, which says nothing;
                // name it after what produced it.
                let mut s = parse_any(&text, label.clone(), Source::SimpleFold)?;
                s.title = format!("SimpleFold {label} · {steps} steps");
                Ok(Outcome::Structure(s))
            }),
        }
    }

    /// Fold an entry's own sequences and show the prediction on its own. Each
    /// distinct chain is folded in turn, since SimpleFold takes one sequence
    /// at a time; identical copies are folded once.
    pub fn fold_entry(&self, what: Retrieval, entry: String) -> Job {
        let engine = self.engine.clone();
        let model = MODELS[self.model].to_string();
        let label = model_label(MODELS[self.model]).to_string();
        let steps = self.steps;
        let plddt = self.plddt;
        let installing = !engine.is_installed();
        Job {
            stage: "Fetching",
            notes: self.preamble(),
            work: Box::new(move |tx| {
                // Only the sequence is wanted, so an accession is looked up at
                // UniProt rather than pulling the whole AlphaFold model down.
                let chains: Vec<(Vec<String>, String)> = match what {
                    Retrieval::AlphaFold(ref acc) => {
                        let sequence = fetch::fetch_uniprot_sequence(acc)?;
                        vec![(vec!["A".to_string()], sequence)]
                    }
                    Retrieval::Pdb(ref id) => fetch::fetch_pdb(id, tx)?.distinct_chains(),
                    Retrieval::PubChem(_) => anyhow::bail!("a compound has no protein chains to fold"),
                };

                let (foldable, short): (Vec<_>, Vec<_>) =
                    chains.into_iter().partition(|(_, seq)| seq.len() >= 12);
                for (ids, seq) in &short {
                    let _ = tx.send(Event::Log(format!(
                        "Chain {} is too short to fold ({} residues).",
                        ids.join("/"),
                        seq.len()
                    )));
                }
                if foldable.is_empty() {
                    anyhow::bail!("nothing long enough to fold");
                }

                let total = foldable.len();
                if total > 1 {
                    let _ = tx.send(Event::Log(format!(
                        "{total} distinct chains to fold, one at a time."
                    )));
                }
                if installing {
                    engine.install(tx)?;
                }

                let mut parts = Vec::with_capacity(total);
                for (i, (ids, sequence)) in foldable.into_iter().enumerate() {
                    // Copies of one chain share a row: "A/C" is folded once.
                    let name = ids.join("/");
                    if total > 1 {
                        let _ = tx.send(Event::Stage(format!(
                            "Folding chain {name} ({} of {total})",
                            i + 1
                        )));
                    }
                    let file: String = name.chars().filter(|c| c.is_ascii_alphanumeric()).collect();
                    let file = if file.is_empty() {
                        format!("chain{}", i + 1)
                    } else {
                        file
                    };
                    let path = engine.fold(&file, &sequence, &model, steps, plddt, tx)?;
                    let text = std::fs::read_to_string(&path)?;
                    parts.push((
                        name.clone(),
                        parse_any(&text, name, Source::SimpleFold)?,
                    ));
                }

                let _ = tx.send(Event::Stage("Reading structure".into()));
                Ok(Outcome::Structure(assemble(
                    parts,
                    entry.to_uppercase(),
                    format!("SimpleFold {label} · {steps} steps from {}", entry.to_uppercase()),
                )))
            }),
        }
    }

    /// Fetch the deposited structure, fold its own sequence, and lay one over
    /// the other.
    pub fn compare(&self, what: Retrieval) -> Job {
        let engine = self.engine.clone();
        let model = MODELS[self.model].to_string();
        let label = model_label(MODELS[self.model]).to_string();
        let steps = self.steps;
        let plddt = self.plddt;
        let installing = !engine.is_installed();
        Job {
            stage: "Fetching",
            notes: self.preamble(),
            work: Box::new(move |tx| {
                let mut reference = match what.run(tx)? {
                    Outcome::Structure(s) => s,
                    _ => anyhow::bail!("nothing to compare against"),
                };
                let (_, sequence) = reference.primary_chain();
                if sequence.len() < 12 {
                    anyhow::bail!("the longest chain is too short to fold");
                }
                if reference.chains.len() > 1 {
                    let _ = tx.send(Event::Log(format!(
                        "Comparing the longest chain only (1 of {} chains).",
                        reference.chains.len()
                    )));
                }

                if installing {
                    engine.install(tx)?;
                }
                let path = engine.fold("query", &sequence, &model, steps, plddt, tx)?;

                let _ = tx.send(Event::Stage("Superimposing".into()));
                let text = std::fs::read_to_string(&path)?;
                let mut prediction = parse_any(&text, label.clone(), Source::SimpleFold)?;
                prediction.title = format!("SimpleFold {label} · {steps} steps");
                let fit = align::superimpose(&mut reference, &mut prediction)?;
                Ok(Outcome::Comparison {
                    reference,
                    prediction,
                    rmsd: fit.rmsd,
                    matched: fit.matched,
                })
            }),
        }
    }
}

/// The corners of the box holding every atom.
fn bounds(s: &Structure) -> ([f32; 3], [f32; 3]) {
    let mut lo = [f32::MAX; 3];
    let mut hi = [f32::MIN; 3];
    for a in &s.atoms {
        for i in 0..3 {
            lo[i] = lo[i].min(a.pos[i]);
            hi[i] = hi[i].max(a.pos[i]);
        }
    }
    if lo[0] > hi[0] {
        return ([0.0; 3], [0.0; 3]);
    }
    (lo, hi)
}

/// Lay predictions out in a row as one structure, each becoming a chain under
/// the name it was folded from. SimpleFold sees one sequence at a time and
/// knows nothing of the assembly, so the chains stand side by side rather than
/// pretending to be docked.
fn assemble(parts: Vec<(String, Structure)>, id: String, title: String) -> Structure {
    let mut out = Structure::new(id, title, Source::SimpleFold);
    out.has_confidence = !parts.is_empty() && parts.iter().all(|(_, s)| s.has_confidence);

    // Enough air between them to read as separate molecules.
    const GAP: f32 = 8.0;
    let mut cursor = 0.0f32;

    for (name, part) in parts {
        let (lo, hi) = bounds(&part);
        let shift = [
            cursor - lo[0],
            -(lo[1] + hi[1]) / 2.0,
            -(lo[2] + hi[2]) / 2.0,
        ];
        cursor += (hi[0] - lo[0]) + GAP;

        let atom0 = out.atoms.len() as u32;
        let residue0 = out.residues.len() as u32;
        let chain = out.chains.len() as u16;

        out.bonds
            .extend(part.bonds.iter().map(|b| [b[0] + atom0, b[1] + atom0]));
        out.atoms.extend(part.atoms.iter().map(|a| Atom {
            pos: [
                a.pos[0] + shift[0],
                a.pos[1] + shift[1],
                a.pos[2] + shift[2],
            ],
            residue: a.residue + residue0,
            ..*a
        }));

        // The part's own chain order is the order to keep; a prediction has
        // one chain, but the fallback covers anything else that turns up.
        let order: Vec<u32> = if part.chains.is_empty() {
            (0..part.residues.len() as u32).collect()
        } else {
            part.chains
                .iter()
                .flat_map(|c| c.residues.iter().copied())
                .collect()
        };
        out.chains.push(Chain {
            id: name,
            residues: order.iter().map(|i| i + residue0).collect(),
        });
        out.residues.extend(part.residues.iter().map(|r| {
            let mut r = r.clone();
            r.chain = chain;
            r.ca = r.ca.map(|i| i + atom0);
            r.n = r.n.map(|i| i + atom0);
            r.c = r.c.map(|i| i + atom0);
            r.o = r.o.map(|i| i + atom0);
            r
        }));
    }
    out
}

// ---------------------------------------------------------------------------
// View
// ---------------------------------------------------------------------------

/// The list under the MODEL pill. Each row says whether its checkpoint is
/// already downloaded, since the larger ones are several gigabytes and are
/// fetched on first use.
fn model_list(fold: &SimpleFold, cx: &mut Context<Viewer>) -> impl IntoElement {
    let current = fold.model;
    let engine = fold.engine.clone();
    div()
        .absolute()
        .top(px(25.))
        .left_0()
        .w(px(186.))
        .p(px(3.))
        .rounded(px(8.))
        .bg(theme::surface())
        .border_1()
        .border_color(theme::border())
        .shadow_md()
        .flex()
        .flex_col()
        .gap(px(1.))
        // Keep the press from reaching what lies under the list.
        .on_mouse_down(gpui::MouseButton::Left, |_, _, cx| cx.stop_propagation())
        .on_mouse_down_out(cx.listener(|this, _, _, cx| {
            this.simplefold.menu_open = false;
            cx.notify();
        }))
        .children(MODELS.iter().enumerate().map(move |(i, name)| {
            let ready = engine.has_checkpoint(name);
            div()
                .id(("model-option", i))
                .flex()
                .items_baseline()
                .justify_between()
                .gap(px(8.))
                .px(px(8.))
                .py(px(4.))
                .rounded(px(5.))
                .cursor_pointer()
                .when(i == current, |d| d.bg(theme::sunken()))
                .when(i != current, |d| d.hover(|s| s.bg(theme::sunken())))
                .child(
                    div()
                        .text_size(px(11.5))
                        .text_color(if i == current { theme::ink() } else { theme::text() })
                        .child(model_label(name).to_string()),
                )
                .child(
                    div()
                        .text_size(px(9.))
                        .text_color(theme::faint())
                        .child(if ready { "on disk" } else { "downloads" }),
                )
                .on_click(cx.listener(move |this, _, _, cx| {
                    this.simplefold.model = i;
                    this.simplefold.menu_open = false;
                    cx.notify();
                }))
        }))
}

/// Model, steps and pLDDT: the settings for Fold and Compare.
pub fn options(fold: &SimpleFold, cx: &mut Context<Viewer>) -> impl IntoElement {
    let model: SharedString = model_label(MODELS[fold.model]).to_string().into();
    let steps = fold.steps;
    let open = fold.menu_open;
    let plddt = fold.plddt;
    div()
        .flex()
        .items_center()
        .gap(px(14.))
        .px(px(20.))
        .child(label("MODEL"))
        .child(
            div()
                .relative()
                .child(
                    pill("model", format!("{model}  ▾"), open)
                        .debug_selector(|| "model".into())
                        .when(!open, |d| d.border_1().border_color(theme::border()))
                        .cursor_pointer()
                        .on_click(cx.listener(|this, _, _, cx| {
                            this.simplefold.menu_open = !this.simplefold.menu_open;
                            cx.notify();
                        })),
                )
                // Painted last so it sits over whatever follows it.
                .when(open, |d| d.child(deferred(model_list(fold, cx)).with_priority(1))),
        )
        .child(label("STEPS"))
        .child(
            pill("steps", steps.to_string(), false)
                .border_1()
                .border_color(theme::border())
                .cursor_pointer()
                .on_click(cx.listener(|this, _, _, cx| {
                    let f = &mut this.simplefold;
                    f.steps = match f.steps {
                        100 => 200,
                        200 => 500,
                        _ => 100,
                    };
                    cx.notify();
                })),
        )
        .child(
            pill("plddt", "pLDDT", plddt)
                .when(!plddt, |d| d.border_1().border_color(theme::border()))
                .cursor_pointer()
                .on_click(cx.listener(|this, _, _, cx| {
                    this.simplefold.plddt = !this.simplefold.plddt;
                    cx.notify();
                })),
        )
}
