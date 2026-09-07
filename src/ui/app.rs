//! The root view: one input on the left, the molecule on the right.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::mpsc::{self, Sender};
use std::sync::Arc;
use std::time::Duration;

use gpui::{
    div, prelude::*, px, relative, Animation, AnimationExt, App, ClipboardItem, Context, Entity,
    FocusHandle, Focusable, MouseButton, MouseDownEvent, MouseMoveEvent, MouseUpEvent, Pixels,
    Point, ScrollDelta, ScrollWheelEvent, SharedString, Task, Window,
};

use crate::align;
use crate::fetch::{self, Input};
use crate::fold::{self, Engine, MODELS};
use crate::job::{self, Event, Outcome};
use crate::protein::{parse_any, Source, Structure};
use crate::render::{self, Camera, Coloring, Palette, Repr, Scene, Style};
use crate::theme;
use crate::ui::input::{InputEvent, TextInput};
use crate::ui::viewer;

const PANEL_WIDTH: f32 = 372.0;
const TITLEBAR: f32 = 42.0;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Drag {
    Orbit,
    Pan,
}

pub struct Profold {
    focus: FocusHandle,
    input: Entity<TextInput>,
    detected: Input,

    structure: Option<Arc<Structure>>,
    /// The prediction laid over `structure`, when comparing.
    compare: Option<Arc<Structure>>,
    fit: Option<(f32, usize)>,
    show_reference: bool,
    show_prediction: bool,
    scene: Arc<Scene>,
    camera: Camera,
    repr: Repr,
    coloring: Coloring,
    style: Style,
    palette: Palette,
    cache: Rc<RefCell<viewer::Cache>>,

    drag: Option<Drag>,
    last_mouse: Point<Pixels>,

    engine: Arc<Engine>,
    model: usize,
    steps: u32,
    plddt: bool,

    copied: bool,
    copy_reset: Option<Task<()>>,
    busy: Option<SharedString>,
    progress: Option<job::Progress>,
    log: Vec<SharedString>,
    error: Option<SharedString>,
    task: Option<Task<()>>,
}

impl Profold {
    pub fn new(window: &mut Window, cx: &mut Context<Self>) -> Self {
        let input = cx.new(|cx| {
            TextInput::new("Paste a sequence, PDB ID, or UniProt accession", cx)
        });
        cx.subscribe_in(&input, window, |this, _, event, _window, cx| match event {
            InputEvent::Changed => this.reclassify(cx),
            InputEvent::Submit => this.submit(cx),
        })
        .detach();

        let focus = cx.focus_handle();
        window.focus(&input.read(cx).focus_handle(cx));

        Self {
            focus,
            input,
            detected: Input::Empty,
            structure: None,
            compare: None,
            fit: None,
            show_reference: true,
            show_prediction: true,
            scene: Arc::new(Scene::default()),
            camera: Camera {
                yaw: 0.0,
                pitch: 0.0,
                zoom: 1.0,
                pan: [0.0, 0.0],
                center: [0.0; 3],
                radius: 1.0,
            },
            repr: Repr::Cartoon,
            coloring: Coloring::Chain,
            style: Style::Soft,
            palette: Palette::Spectrum,
            cache: Rc::new(RefCell::new(viewer::Cache::default())),
            drag: None,
            last_mouse: Point::default(),
            engine: Arc::new(Engine::new()),
            model: 0,
            steps: 200,
            plddt: false,
            copied: false,
            copy_reset: None,
            busy: None,
            progress: None,
            log: Vec::new(),
            error: None,
            task: None,
        }
    }

    // -- state ------------------------------------------------------------

    fn reclassify(&mut self, cx: &mut Context<Self>) {
        self.detected = fetch::classify(self.input.read(cx).text());
        self.error = None;
        cx.notify();
    }

    fn note(&mut self, line: impl Into<SharedString>) {
        self.log.push(line.into());
        if self.log.len() > 200 {
            self.log.drain(..self.log.len() - 200);
        }
    }

    fn show(&mut self, structure: Structure, cx: &mut Context<Self>) {
        self.coloring = if structure.has_confidence {
            Coloring::Confidence
        } else {
            Coloring::Chain
        };
        self.camera = Camera::frame(&structure);
        self.structure = Some(Arc::new(structure));
        self.compare = None;
        self.fit = None;
        self.rebuild(cx);
    }

    fn show_comparison(
        &mut self,
        mut reference: Structure,
        mut prediction: Structure,
        rmsd: f32,
        matched: usize,
        cx: &mut Context<Self>,
    ) {
        reference.tint = Some(align::REFERENCE_TINT);
        prediction.tint = Some(align::PREDICTION_TINT);
        self.coloring = Coloring::Pair;
        self.camera = Camera::frame_all(&[&reference, &prediction]);
        self.structure = Some(Arc::new(reference));
        self.compare = Some(Arc::new(prediction));
        self.fit = Some((rmsd, matched));
        self.show_reference = true;
        self.show_prediction = true;
        self.rebuild(cx);
    }

    /// The structures currently drawn, honouring the comparison toggles.
    fn visible(&self) -> Vec<&Arc<Structure>> {
        let mut out = Vec::new();
        if self.show_reference || self.compare.is_none() {
            out.extend(self.structure.iter());
        }
        if self.show_prediction {
            out.extend(self.compare.iter());
        }
        out
    }

    fn rebuild(&mut self, cx: &mut Context<Self>) {
        let residues: usize = self.visible().iter().map(|s| s.residues.len()).sum();
        // Long chains get a coarser tube so that dragging stays fluid.
        let quality = match residues {
            0..=800 => 1.0,
            801..=2200 => 0.6,
            _ => 0.3,
        };
        let mut scene = Scene::default();
        for s in self.visible() {
            scene.merge(render::build_scene(
                s,
                self.repr,
                self.coloring,
                quality,
                self.style,
                self.palette,
            ));
        }
        self.scene = Arc::new(scene);
        cx.notify();
    }

    fn finish(&mut self, result: anyhow::Result<Structure>, cx: &mut Context<Self>) {
        self.busy = None;
        self.progress = None;
        match result {
            Ok(s) => {
                self.note(format!(
                    "{} · {} residues · {} atoms",
                    s.id,
                    s.residues.len(),
                    s.atoms.len()
                ));
                self.show(s, cx);
            }
            Err(e) => {
                self.error = Some(format!("{e:#}").into());
                cx.notify();
            }
        }
    }

    // -- actions ----------------------------------------------------------

    fn submit(&mut self, cx: &mut Context<Self>) {
        if self.busy.is_some() {
            return;
        }
        self.error = None;
        match self.detected.clone() {
            Input::Empty | Input::Unknown => {}
            Input::Pdb(id) => self.retrieve(Retrieval::Pdb(id), cx),
            Input::UniProt(acc) => self.retrieve(Retrieval::AlphaFold(acc), cx),
            Input::Sequence(seq) => self.fold(seq, cx),
        }
    }

    fn submit_comparison(&mut self, cx: &mut Context<Self>) {
        if self.busy.is_some() {
            return;
        }
        self.error = None;
        match self.detected.clone() {
            Input::Pdb(id) => self.compare(Retrieval::Pdb(id), cx),
            Input::UniProt(acc) => self.compare(Retrieval::AlphaFold(acc), cx),
            _ => {}
        }
    }

    /// Run one background job, streaming its stages and progress into the UI.
    fn run_job<F>(&mut self, stage: &str, work: F, cx: &mut Context<Self>)
    where
        F: FnOnce(&Sender<Event>) -> anyhow::Result<Outcome> + Send + 'static,
    {
        self.busy = Some(stage.to_string().into());
        self.progress = None;
        self.error = None;
        cx.notify();

        let (tx, rx) = mpsc::channel::<Event>();
        std::thread::spawn(move || {
            let result = work(&tx);
            let _ = tx.send(Event::Finished(result));
        });

        self.task = Some(cx.spawn(async move |this, cx| {
            loop {
                let mut dirty = false;
                let mut done = None;
                loop {
                    match rx.try_recv() {
                        Ok(Event::Stage(s)) => {
                            let _ = this.update(cx, |this, _| {
                                this.note(s.clone());
                                this.busy = Some(s.into());
                                this.progress = None;
                            });
                            dirty = true;
                        }
                        Ok(Event::Log(line)) => {
                            let _ = this.update(cx, |this, _| this.note(line));
                            dirty = true;
                        }
                        Ok(Event::Progress(p)) => {
                            let _ = this.update(cx, |this, _| this.progress = Some(p));
                            dirty = true;
                        }
                        Ok(Event::Finished(result)) => {
                            done = Some(result);
                            break;
                        }
                        Err(mpsc::TryRecvError::Empty) => break,
                        Err(mpsc::TryRecvError::Disconnected) => {
                            done = Some(Err(anyhow::anyhow!("the job stopped unexpectedly")));
                            break;
                        }
                    }
                }
                if dirty {
                    let _ = this.update(cx, |_, cx| cx.notify());
                }
                if let Some(result) = done {
                    let _ = this.update(cx, |this, cx| this.settle(result, cx));
                    return;
                }
                cx.background_executor()
                    .timer(Duration::from_millis(70))
                    .await;
            }
        }));
    }

    fn settle(&mut self, result: anyhow::Result<Outcome>, cx: &mut Context<Self>) {
        match result {
            Ok(Outcome::Structure(s)) => self.finish(Ok(s), cx),
            Ok(Outcome::Sequence(seq)) => {
                self.note("No AlphaFold model; folding the UniProt sequence instead.");
                self.fold(seq, cx);
            }
            Ok(Outcome::Comparison {
                reference,
                prediction,
                rmsd,
                matched,
            }) => {
                self.busy = None;
                self.progress = None;
                self.note(format!("RMSD {rmsd:.2} A over {matched} Ca"));
                self.show_comparison(reference, prediction, rmsd, matched, cx);
            }
            Err(e) => self.finish(Err(e), cx),
        }
    }

    fn retrieve(&mut self, what: Retrieval, cx: &mut Context<Self>) {
        self.run_job("Fetching", move |tx| what.run(tx), cx);
    }

    /// Fetch the deposited structure, fold its own sequence, and lay one over
    /// the other.
    fn compare(&mut self, what: Retrieval, cx: &mut Context<Self>) {
        let engine = self.engine.clone();
        let model = MODELS[self.model].to_string();
        let label = fold::model_label(MODELS[self.model]).to_string();
        let steps = self.steps;
        let plddt = self.plddt;
        let installing = !engine.is_installed();
        if installing {
            self.note("First run: installing SimpleFold. This takes a while.");
        }

        self.run_job(
            "Fetching",
            move |tx| {
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
                        "Comparing the longest chain only ({} of {} chains).",
                        1,
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
            },
            cx,
        );
    }

    fn fold(&mut self, sequence: String, cx: &mut Context<Self>) {
        let engine = self.engine.clone();
        let model = MODELS[self.model].to_string();
        let label = fold::model_label(MODELS[self.model]).to_string();
        let steps = self.steps;
        let plddt = self.plddt;
        let installing = !engine.is_installed();

        if installing {
            self.note("First run: installing SimpleFold. This takes a while.");
        }
        if plddt && !engine.has_plddt_weights() {
            self.note("pLDDT scoring needs SimpleFold's 1.6B checkpoint (~6.5 GB), once.");
        }

        let opening = if installing {
            "Setting up the folding engine"
        } else {
            "Folding"
        };
        self.run_job(
            opening,
            move |tx| {
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
            },
            cx,
        );
    }

    fn copy_sequence(&mut self, sequence: String, cx: &mut Context<Self>) {
        cx.write_to_clipboard(ClipboardItem::new_string(sequence));
        self.copied = true;
        cx.notify();
        // A separate task: reusing the job task would cancel a running fold.
        self.copy_reset = Some(cx.spawn(async move |this, cx| {
            cx.background_executor()
                .timer(Duration::from_millis(1400))
                .await;
            let _ = this.update(cx, |this, cx| {
                this.copied = false;
                cx.notify();
            });
        }));
    }

    /// Load an identifier or sequence as though it had been typed.
    pub fn open(&mut self, text: &str, cx: &mut Context<Self>) {
        let text = text.to_string();
        self.input.update(cx, |input, cx| input.set_text(text, cx));
        self.reclassify(cx);
        self.submit(cx);
    }

    /// Load an identifier and immediately compare it against a prediction.
    pub fn open_comparison(&mut self, text: &str, cx: &mut Context<Self>) {
        let text = text.to_string();
        self.input.update(cx, |input, cx| input.set_text(text, cx));
        self.reclassify(cx);
        self.submit_comparison(cx);
    }

    // -- camera -----------------------------------------------------------

    fn on_viewport_down(&mut self, e: &MouseDownEvent, _: &mut Window, cx: &mut Context<Self>) {
        if e.click_count >= 2 {
            if let Some(s) = &self.structure {
                self.camera = Camera::frame(s);
            }
        } else {
            self.drag = Some(if e.modifiers.shift || e.button == MouseButton::Right {
                Drag::Pan
            } else {
                Drag::Orbit
            });
            self.last_mouse = e.position;
        }
        cx.notify();
    }

    fn on_mouse_move(&mut self, e: &MouseMoveEvent, _: &mut Window, cx: &mut Context<Self>) {
        let Some(drag) = self.drag else { return };
        let dx = f32::from(e.position.x - self.last_mouse.x);
        let dy = f32::from(e.position.y - self.last_mouse.y);
        self.last_mouse = e.position;
        match drag {
            Drag::Orbit => {
                // Negated so the molecule follows the pointer: dragging right
                // brings the near face right, rather than spinning away.
                self.camera.yaw -= dx * 0.008;
                self.camera.pitch = (self.camera.pitch - dy * 0.008).clamp(-1.55, 1.55);
            }
            Drag::Pan => {
                self.camera.pan[0] += dx;
                self.camera.pan[1] += dy;
            }
        }
        cx.notify();
    }

    fn on_mouse_up(&mut self, _: &MouseUpEvent, _: &mut Window, cx: &mut Context<Self>) {
        if self.drag.take().is_some() {
            cx.notify();
        }
    }

    fn on_scroll(&mut self, e: &ScrollWheelEvent, _: &mut Window, cx: &mut Context<Self>) {
        let dy = match e.delta {
            ScrollDelta::Pixels(p) => f32::from(p.y),
            ScrollDelta::Lines(p) => p.y * 18.0,
        };
        self.camera.zoom = (self.camera.zoom * (1.0 + dy * 0.0022)).clamp(0.12, 24.0);
        cx.notify();
    }
}

enum Retrieval {
    Pdb(String),
    AlphaFold(String),
}

impl Retrieval {
    fn run(self, tx: &Sender<Event>) -> anyhow::Result<Outcome> {
        match self {
            Retrieval::Pdb(id) => fetch::fetch_pdb(&id, tx).map(Outcome::Structure),
            Retrieval::AlphaFold(acc) => match fetch::fetch_alphafold(&acc, tx) {
                Ok(s) => Ok(Outcome::Structure(s)),
                // An accession without a model can still be folded.
                Err(e) => match fetch::fetch_uniprot_sequence(&acc) {
                    Ok(seq) => Ok(Outcome::Sequence(seq)),
                    Err(_) => Err(e),
                },
            },
        }
    }
}

impl Focusable for Profold {
    fn focus_handle(&self, _: &App) -> FocusHandle {
        self.focus.clone()
    }
}

// ---------------------------------------------------------------------------
// View
// ---------------------------------------------------------------------------

const EXAMPLES: [(&str, &str); 3] = [
    ("1CRN", "crambin, 46 aa"),
    ("P69905", "human haemoglobin α"),
    (
        "GSHMKAEELKKKAEELLKKAEELAKKLAEELKKLAEELKKKAEELLKKLA",
        "a designed helix, folded here",
    ),
];

fn pill(id: &'static str, label: impl Into<SharedString>, active: bool) -> gpui::Stateful<gpui::Div> {
    div()
        .id(id)
        .px(px(9.))
        .py(px(3.))
        .rounded(px(5.))
        .text_size(px(11.))
        .when(active, |d| d.bg(theme::INK).text_color(theme::ON_INK))
        .when(!active, |d| {
            d.text_color(theme::MUTED)
                .hover(|s| s.bg(theme::SUNKEN).text_color(theme::TEXT))
        })
        .child(label.into())
}

fn label(text: &'static str) -> impl IntoElement {
    div()
        .text_size(px(9.5))
        .text_color(theme::FAINT)
        .child(text)
}


fn swatch(c: [f32; 3]) -> gpui::Rgba {
    gpui::Rgba {
        r: c[0],
        g: c[1],
        b: c[2],
        a: 1.0,
    }
}

/// A pill with the structure's own colour beside it, used to show and hide
/// each half of a comparison.
fn dot_pill(
    id: &'static str,
    color: [f32; 3],
    label: String,
    on: bool,
) -> gpui::Stateful<gpui::Div> {
    div()
        .id(id)
        .flex()
        .items_center()
        .gap(px(6.))
        .px(px(9.))
        .py(px(3.))
        .rounded(px(5.))
        .cursor_pointer()
        .when(on, |d| d.bg(theme::SUNKEN))
        .child(
            div()
                .size(px(7.))
                .rounded_full()
                .bg(swatch(color))
                .when(!on, |d| d.opacity(0.28)),
        )
        .child(
            div()
                .text_size(px(11.))
                .text_color(if on { theme::TEXT } else { theme::FAINT })
                .child(label),
        )
}

/// A one-line label on a dark ground, for hover hints.
struct Tip(SharedString);

impl Render for Tip {
    fn render(&mut self, _: &mut Window, _: &mut Context<Self>) -> impl IntoElement {
        div()
            .px(px(7.))
            .py(px(3.))
            .rounded(px(5.))
            .bg(theme::INK)
            .text_size(px(10.5))
            .text_color(theme::ON_INK)
            .child(self.0.clone())
    }
}

/// A palette is shown as its own ramp; the name is on hover.
fn swatch_strip(p: Palette) -> impl IntoElement {
    div()
        .flex()
        .rounded(px(3.))
        .overflow_hidden()
        .children((0..7).map(move |i| {
            div()
                .w(px(5.))
                .h(px(12.))
                .bg(swatch(p.ramp(i as f32 / 6.0)))
        }))
}


fn rule() -> impl IntoElement {
    div().h(px(1.)).w_full().bg(theme::BORDER)
}

impl Profold {
    fn header(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let _ = cx;
        div()
            .id("titlebar")
            .flex()
            .flex_col()
            .gap(px(2.))
            .px(px(20.))
            .pb(px(16.))
            .on_mouse_down(MouseButton::Left, |_, window, _| window.start_window_move())
            .child(
                div()
                    .text_size(px(15.))
                    .text_color(theme::INK)
                    .child("profold"),
            )
            .child(
                div()
                    .text_size(px(10.5))
                    .text_color(theme::FAINT)
                    .child("SimpleFold · structure viewer"),
            )
    }

    fn search(&self, focused: bool, cx: &mut Context<Self>) -> impl IntoElement {
        let detected = self.detected.clone();
        let ready = !matches!(detected, Input::Empty | Input::Unknown);
        let can_compare = matches!(detected, Input::Pdb(_) | Input::UniProt(_));
        let busy = self.busy.is_some();

        div()
            .flex()
            .flex_col()
            .gap(px(10.))
            .px(px(20.))
            .child(
                div()
                    .flex()
                    .items_center()
                    .h(px(38.))
                    .px(px(11.))
                    .rounded(px(8.))
                    .bg(if focused { theme::SURFACE } else { theme::SUNKEN })
                    .border_1()
                    .border_color(if focused {
                        theme::BORDER_FOCUS
                    } else {
                        theme::BORDER
                    })
                    .text_size(px(12.5))
                    .text_color(theme::TEXT)
                    .font_family(theme::MONO)
                    .child(self.input.clone()),
            )
            .child(
                div()
                    .flex()
                    .items_center()
                    .justify_between()
                    .gap(px(8.))
                    .h(px(22.))
                    .child(match &detected {
                        Input::Empty => div().child(""),
                        Input::Unknown => div()
                            .text_size(px(10.5))
                            .text_color(theme::FAINT)
                            .child(detected.detail()),
                        other => div()
                            .flex()
                            .items_center()
                            .gap(px(6.))
                            .child(
                                div()
                                    .text_size(px(9.))
                                    .text_color(theme::MUTED)
                                    .child(other.kind()),
                            )
                            .child(
                                div()
                                    .font_family(theme::MONO)
                                    .text_size(px(11.))
                                    .text_color(theme::TEXT)
                                    .child(other.detail()),
                            ),
                    })
                    .child(
                        div()
                            .flex()
                            .items_center()
                            .gap(px(6.))
                            .when(can_compare && !busy, |d| {
                                d.child(
                                    div()
                                        .id("compare")
                                        .px(px(11.))
                                        .py(px(4.))
                                        .rounded(px(6.))
                                        .text_size(px(11.5))
                                        .border_1()
                                        .border_color(theme::BORDER)
                                        .text_color(theme::MUTED)
                                        .cursor_pointer()
                                        .hover(|s| s.bg(theme::SUNKEN).text_color(theme::TEXT))
                                        .child("Compare")
                                        .on_click(cx.listener(|this, _, _, cx| {
                                            this.submit_comparison(cx)
                                        })),
                                )
                            })
                            .child(
                                div()
                                    .id("run")
                                    .px(px(13.))
                                    .py(px(5.))
                                    .rounded(px(6.))
                                    .text_size(px(11.5))
                                    .when(ready && !busy, |d| {
                                        d.bg(theme::INK)
                                            .text_color(theme::ON_INK)
                                            .hover(|s| s.opacity(0.86))
                                            .cursor_pointer()
                                    })
                                    .when(!ready || busy, |d| {
                                        d.bg(theme::SUNKEN).text_color(theme::FAINT)
                                    })
                                    .child(detected.action())
                                    .on_click(cx.listener(|this, _, _, cx| this.submit(cx))),
                            ),
                    ),
            )
    }

    fn options(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let model = fold::model_label(MODELS[self.model]).to_string();
        let steps = self.steps;
        div()
            .flex()
            .items_center()
            .gap(px(14.))
            .px(px(20.))
            .child(label("MODEL"))
            .child(
                pill("model", model, false)
                    .border_1()
                    .border_color(theme::BORDER)
                    .cursor_pointer()
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.model = (this.model + 1) % MODELS.len();
                        cx.notify();
                    })),
            )
            .child(label("STEPS"))
            .child(
                pill("steps", steps.to_string(), false)
                    .border_1()
                    .border_color(theme::BORDER)
                    .cursor_pointer()
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.steps = match this.steps {
                            100 => 200,
                            200 => 500,
                            _ => 100,
                        };
                        cx.notify();
                    })),
            )
            .child(
                pill("plddt", "pLDDT", self.plddt)
                    .when(!self.plddt, |d| d.border_1().border_color(theme::BORDER))
                    .cursor_pointer()
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.plddt = !this.plddt;
                        cx.notify();
                    })),
            )
    }

    fn recent(&self, lines: usize) -> Vec<SharedString> {
        self.log.iter().rev().take(lines).rev().cloned().collect()
    }

    /// Left panel: only what went wrong. The running log belongs to the
    /// working panel, which is where you look while a job is going.
    fn status(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let _ = cx;
        div()
            .px(px(20.))
            .when_some(self.error.clone(), |d, err| {
                d.child(div().text_size(px(11.)).text_color(theme::DANGER).child(err))
            })
    }

    /// A determinate bar where the job knows how far along it is, and a
    /// travelling one where it does not, so the window never looks stalled.
    fn bar(&self) -> impl IntoElement {
        let fraction = self.progress.as_ref().and_then(|p| p.fraction);
        div()
            .relative()
            .h(px(4.))
            .w_full()
            .rounded_full()
            .bg(theme::SUNKEN)
            .overflow_hidden()
            .child(match fraction {
                Some(f) => div()
                    .h_full()
                    .w(relative(f.clamp(0.005, 1.0)))
                    .rounded_full()
                    .bg(theme::INK)
                    .into_any_element(),
                None => div()
                    .id("indeterminate")
                    .absolute()
                    .top_0()
                    .bottom_0()
                    .w(relative(0.28))
                    .rounded_full()
                    .bg(theme::INK)
                    .with_animation(
                        "travel",
                        Animation::new(Duration::from_millis(1600)).repeat(),
                        |el, delta| el.left(relative(-0.28 + delta * 1.28)),
                    )
                    .into_any_element(),
            })
    }

    /// Shown in place of the viewer while a job runs.
    fn working(&self, busy: SharedString) -> impl IntoElement {
        let percent = self
            .progress
            .as_ref()
            .and_then(|p| p.fraction)
            .map(|f| format!("{:.0}%", f * 100.0));
        let detail = self
            .progress
            .as_ref()
            .map(|p| p.text.clone())
            .filter(|t| !t.is_empty());

        div()
            .absolute()
            .inset_0()
            .bg(theme::PAPER)
            .flex()
            .items_center()
            .justify_center()
            .child(
                div()
                    .w(px(430.))
                    .flex()
                    .flex_col()
                    .gap(px(11.))
                    .child(
                        div()
                            .flex()
                            .items_baseline()
                            .justify_between()
                            .gap(px(10.))
                            .child(
                                div()
                                    .text_size(px(13.))
                                    .text_color(theme::INK)
                                    .child(busy),
                            )
                            .when_some(percent, |d, p| {
                                d.child(
                                    div()
                                        .font_family(theme::MONO)
                                        .text_size(px(11.5))
                                        .text_color(theme::MUTED)
                                        .child(p),
                                )
                            }),
                    )
                    .child(self.bar())
                    .child(
                        div()
                            .h(px(13.))
                            .overflow_hidden()
                            .whitespace_nowrap()
                            .font_family(theme::MONO)
                            .text_size(px(10.5))
                            .text_color(theme::MUTED)
                            .children(detail),
                    )
                    .child(
                        div()
                            .pt(px(10.))
                            .flex()
                            .flex_col()
                            .gap(px(1.))
                            .font_family(theme::MONO)
                            .text_size(px(10.))
                            .text_color(theme::FAINT)
                            .children(self.recent(6).into_iter().map(|line| {
                                div().overflow_hidden().whitespace_nowrap().child(line)
                            })),
                    ),
            )
    }

    fn examples(&self, cx: &mut Context<Self>) -> impl IntoElement {
        div()
            .flex()
            .flex_col()
            .gap(px(7.))
            .px(px(20.))
            .child(label("TRY"))
            .children(EXAMPLES.iter().enumerate().map(|(i, (value, note))| {
                let value = value.to_string();
                let shown = if value.len() > 16 {
                    format!("{}…", &value[..14])
                } else {
                    value.clone()
                };
                div()
                    .id(("example", i))
                    .flex()
                    .items_baseline()
                    .gap(px(8.))
                    .cursor_pointer()
                    .child(
                        div()
                            .font_family(theme::MONO)
                            .text_size(px(11.))
                            .text_color(theme::TEXT)
                            .child(shown),
                    )
                    .child(
                        div()
                            .text_size(px(10.5))
                            .text_color(theme::FAINT)
                            .child(*note),
                    )
                    .on_click(cx.listener(move |this, _, _, cx| {
                        this.open(&value, cx)
                    }))
            }))
    }

    /// The comparison sits over the viewport rather than in the sidebar: it
    /// describes what is on screen, and both halves can be switched off there.
    fn comparison(&self, cx: &mut Context<Self>) -> Option<impl IntoElement> {
        let (rmsd, matched) = self.fit?;
        let reference = self.structure.as_ref()?;
        let prediction = self.compare.as_ref()?;

        Some(
            div()
                .absolute()
                .top(px(TITLEBAR + 8.))
                .left(px(20.))
                .flex()
                .items_center()
                .gap(px(10.))
                .child(
                    div()
                        .flex()
                        .items_center()
                        .gap(px(2.))
                        .p(px(3.))
                        .rounded(px(8.))
                        .bg(theme::SURFACE)
                        .border_1()
                        .border_color(theme::BORDER)
                        .shadow_sm()
                        .child(
                            dot_pill(
                                "show-ref",
                                align::REFERENCE_TINT,
                                reference.id.to_uppercase(),
                                self.show_reference,
                            )
                            .on_click(cx.listener(|this, _, _, cx| {
                                this.show_reference = !this.show_reference;
                                this.rebuild(cx);
                            })),
                        )
                        .child(
                            dot_pill(
                                "show-pred",
                                align::PREDICTION_TINT,
                                format!("SimpleFold {}", prediction.id),
                                self.show_prediction,
                            )
                            .on_click(cx.listener(|this, _, _, cx| {
                                this.show_prediction = !this.show_prediction;
                                this.rebuild(cx);
                            })),
                        ),
                )
                .child(
                    div()
                        .flex()
                        .items_baseline()
                        .gap(px(6.))
                        .child(
                            div()
                                .font_family(theme::MONO)
                                .text_size(px(13.))
                                .text_color(theme::INK)
                                .child(format!("{rmsd:.2} Å")),
                        )
                        .child(
                            div()
                                .text_size(px(10.))
                                .text_color(theme::MUTED)
                                .child(format!("Cα RMSD · {matched} residues")),
                        ),
                ),
        )
    }

    fn info(&self, structure: &Structure, cx: &mut Context<Self>) -> impl IntoElement {
        let sequence = structure.sequence();
        // Enough to recognise it by: both termini, and the length. The whole
        // thing goes to the clipboard.
        let preview = if sequence.len() <= 34 {
            sequence.clone()
        } else {
            format!(
                "{} … {}",
                &sequence[..20],
                &sequence[sequence.len() - 10..]
            )
        };
        let title = if structure.title.is_empty() {
            structure.id.clone()
        } else {
            structure.title.clone()
        };
        let mut meta = vec![
            structure.source.label().to_string(),
            format!(
                "chain{} {}",
                if structure.chains.len() == 1 { "" } else { "s" },
                structure
                    .chains
                    .iter()
                    .map(|c| c.id.clone())
                    .collect::<Vec<_>>()
                    .join(" ")
            ),
            match (structure.residues.first(), structure.residues.last()) {
                (Some(a), Some(b)) => {
                    format!("{} residues ({}–{})", structure.residues.len(), a.seq, b.seq)
                }
                _ => format!("{} residues", structure.residues.len()),
            },
            format!("{} atoms", structure.atoms.len()),
        ];
        if let Some(c) = structure.mean_confidence() {
            meta.push(format!("pLDDT {c:.0}"));
        }
        let meta = meta.join("  ·  ");

        div()
            .flex()
            .flex_col()
            .gap(px(10.))
            .px(px(20.))
            .pb(px(20.))
            .child(rule())
            .child(
                div()
                    .pt(px(4.))
                    .text_size(px(12.))
                    .text_color(theme::TEXT)
                    .child(title),
            )
            .child(
                div()
                    .text_size(px(10.5))
                    .text_color(theme::MUTED)
                    .child(meta),
            )
            .child(
                div()
                    .flex()
                    .items_center()
                    .justify_between()
                    .gap(px(8.))
                    .child(
                        div()
                            .flex()
                            .items_baseline()
                            .gap(px(7.))
                            .child(label("SEQUENCE"))
                            .child(
                                div()
                                    .font_family(theme::MONO)
                                    .text_size(px(10.))
                                    .text_color(theme::MUTED)
                                    .child(format!("{} aa", sequence.len())),
                            ),
                    )
                    .child(
                        div()
                            .id("copy-sequence")
                            .px(px(7.))
                            .py(px(2.))
                            .rounded(px(4.))
                            .text_size(px(10.))
                            .cursor_pointer()
                            .when(self.copied, |d| {
                                d.text_color(theme::TEXT).bg(theme::SUNKEN)
                            })
                            .when(!self.copied, |d| {
                                d.text_color(theme::MUTED)
                                    .hover(|s| s.bg(theme::SUNKEN).text_color(theme::TEXT))
                            })
                            .child(if self.copied { "Copied" } else { "Copy" })
                            .on_click(cx.listener(move |this, _, _, cx| {
                                this.copy_sequence(sequence.clone(), cx)
                            })),
                    ),
            )
            .child(
                div()
                    .overflow_hidden()
                    .whitespace_nowrap()
                    .font_family(theme::MONO)
                    .text_size(px(10.))
                    .text_color(theme::FAINT)
                    .child(preview),
            )
    }

    fn controls(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let repr = self.repr;
        let coloring = self.coloring;
        let style = self.style;
        let palette = self.palette;
        let has_confidence = self
            .structure
            .as_ref()
            .is_some_and(|s| s.has_confidence);
        let options: &[Coloring] = if self.compare.is_some() {
            &Coloring::COMPARE
        } else {
            &Coloring::ALL
        };

        let group = |children: Vec<gpui::AnyElement>| {
            div()
                .flex()
                .items_center()
                .gap(px(2.))
                .p(px(3.))
                .rounded(px(8.))
                .bg(theme::SURFACE)
                .border_1()
                .border_color(theme::BORDER)
                .shadow_sm()
                .children(children)
        };

        div()
            .absolute()
            .bottom(px(20.))
            .left_0()
            .right_0()
            .flex()
            .justify_center()
            .gap(px(8.))
            .child(group(
                Palette::ALL
                    .iter()
                    .map(|p| {
                        let p = *p;
                        div()
                            .id(match p {
                                Palette::Spectrum => "pal-spectrum",
                                Palette::Viridis => "pal-viridis",
                                Palette::Ember => "pal-ember",
                                Palette::Ice => "pal-ice",
                            })
                            .p(px(3.))
                            .rounded(px(5.))
                            .cursor_pointer()
                            .when(p == palette, |d| d.bg(theme::INK))
                            .when(p != palette, |d| d.hover(|s| s.bg(theme::SUNKEN)))
                            .child(swatch_strip(p))
                            .tooltip(move |_, cx| cx.new(|_| Tip(p.label().into())).into())
                            .on_click(cx.listener(move |this, _, _, cx| {
                                this.palette = p;
                                this.rebuild(cx);
                            }))
                            .into_any_element()
                    })
                    .collect(),
            ))
            .child(group(
                Style::ALL
                    .iter()
                    .map(|s| {
                        let s = *s;
                        pill(
                            match s {
                                Style::Soft => "style-soft",
                                Style::Illustrative => "style-ink",
                            },
                            s.label(),
                            s == style,
                        )
                        .cursor_pointer()
                        .on_click(cx.listener(move |this, _, _, cx| {
                            this.style = s;
                            this.rebuild(cx);
                        }))
                        .into_any_element()
                    })
                    .collect(),
            ))
            .child(group(
                Repr::ALL
                    .iter()
                    .map(|r| {
                        let r = *r;
                        pill(
                            match r {
                                Repr::Cartoon => "repr-cartoon",
                                Repr::Tube => "repr-tube",
                                Repr::Spacefill => "repr-fill",
                            },
                            r.label(),
                            r == repr,
                        )
                        .cursor_pointer()
                        .on_click(cx.listener(move |this, _, _, cx| {
                            this.repr = r;
                            this.rebuild(cx);
                        }))
                        .into_any_element()
                    })
                    .collect(),
            ))
            .child(group(
                options
                    .iter()
                    .map(|c| {
                        let c = *c;
                        let enabled = c != Coloring::Confidence || has_confidence;
                        pill(
                            match c {
                                Coloring::Chain => "col-chain",
                                Coloring::Rainbow => "col-rainbow",
                                Coloring::Confidence => "col-plddt",
                                Coloring::Element => "col-element",
                                Coloring::Pair => "col-pair",
                                Coloring::Deviation => "col-deviation",
                            },
                            c.label(),
                            c == coloring,
                        )
                        .when(!enabled, |d| d.text_color(theme::FAINT))
                        .when(enabled, |d| d.cursor_pointer())
                        .on_click(cx.listener(move |this, _, _, cx| {
                            if enabled {
                                this.coloring = c;
                                this.rebuild(cx);
                            }
                        }))
                        .into_any_element()
                    })
                    .collect(),
            ))
    }

    fn stage(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let has_structure = self.structure.is_some();
        let busy: SharedString = self.busy.clone().unwrap_or_default();
        // While comparing, the strip along the top says all of this already.
        let badge = self
            .structure
            .as_ref()
            .filter(|_| self.compare.is_none())
            .map(|s| format!("{} · {}", s.id.to_uppercase(), s.source.label()));

        div()
            .id("stage")
            .relative()
            .flex_1()
            .h_full()
            .overflow_hidden()
            .on_mouse_down(MouseButton::Left, cx.listener(Self::on_viewport_down))
            .on_mouse_down(MouseButton::Right, cx.listener(Self::on_viewport_down))
            .on_scroll_wheel(cx.listener(Self::on_scroll))
            .child(viewer::viewport(
                self.scene.clone(),
                self.camera,
                self.cache.clone(),
                self.drag.is_some(),
                self.style,
            ))
            .children(self.comparison(cx))
            .when_some(badge, |d, badge| {
                d.child(
                    div()
                        .absolute()
                        .top(px(TITLEBAR + 6.))
                        .right(px(20.))
                        .font_family(theme::MONO)
                        .text_size(px(10.))
                        .text_color(theme::FAINT)
                        .child(badge),
                )
            })
            .when(!has_structure && self.busy.is_none(), |d| {
                d.child(
                    div()
                        .absolute()
                        .inset_0()
                        .flex()
                        .flex_col()
                        .items_center()
                        .justify_center()
                        .gap(px(6.))
                        .child(
                            div()
                                .text_size(px(12.))
                                .text_color(theme::FAINT)
                                .child("Nothing loaded"),
                        )
                        .child(
                            div()
                                .text_size(px(10.5))
                                .text_color(theme::FAINT)
                                .child("Paste a sequence or an identifier to begin"),
                        ),
                )
            })
            .when(self.busy.is_some(), |d| d.child(self.working(busy.clone())))
            .when(
                has_structure && self.busy.is_none() && self.visible().is_empty(),
                |d| {
                    d.child(
                        div()
                            .absolute()
                            .inset_0()
                            .flex()
                            .items_center()
                            .justify_center()
                            .text_size(px(11.))
                            .text_color(theme::FAINT)
                            .child("Both structures are hidden"),
                    )
                },
            )
            .when(has_structure && self.busy.is_none(), |d| {
                d.child(self.controls(cx)).child(
                    div()
                        .absolute()
                        .bottom(px(4.))
                        .left_0()
                        .right_0()
                        .flex()
                        .justify_center()
                        .text_size(px(9.5))
                        .text_color(theme::FAINT)
                        .child("drag to orbit · shift-drag to pan · scroll to zoom · double-click to reset"),
                )
            })
    }
}

impl Render for Profold {
    fn render(&mut self, window: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let focused = self.input.read(cx).focus_handle(cx).is_focused(window);
        let sequence_mode = matches!(self.detected, Input::Sequence(_));
        let structure = self.structure.clone();

        div()
            .flex()
            .size_full()
            .bg(theme::PAPER)
            .text_color(theme::TEXT)
            .track_focus(&self.focus)
            .on_mouse_move(cx.listener(Self::on_mouse_move))
            .on_mouse_up(MouseButton::Left, cx.listener(Self::on_mouse_up))
            .on_mouse_up_out(MouseButton::Left, cx.listener(Self::on_mouse_up))
            .on_mouse_up(MouseButton::Right, cx.listener(Self::on_mouse_up))
            .child(
                div()
                    .w(px(PANEL_WIDTH))
                    .h_full()
                    .flex()
                    .flex_col()
                    .gap(px(16.))
                    .pt(px(TITLEBAR))
                    .bg(theme::SURFACE)
                    .border_r_1()
                    .border_color(theme::BORDER)
                    .child(self.header(cx))
                    .child(self.search(focused, cx))
                    .when(sequence_mode, |d| d.child(self.options(cx)))
                    .child(self.status(cx))
                    .when(structure.is_none(), |d| {
                        d.child(self.examples(cx)).child(div().flex_1())
                    })
                    .when_some(structure, |d, s| d.child(self.info(&s, cx)))
                    .child(div().flex_1()),
            )
            .child(self.stage(cx))
    }
}
