//! The root view: one input on the left, the molecule on the right.

use std::cell::{Cell as Shared, RefCell};
use std::collections::HashMap;
use std::path::PathBuf;
use std::rc::Rc;
use std::sync::mpsc::{self, Sender};
use std::sync::Arc;
use std::time::Duration;

use gpui::{
    actions, canvas, deferred, div, prelude::*, px, relative, Animation, AnimationExt, App,
    Bounds, ClipboardItem, Context, Entity, ExternalPaths, FocusHandle, Focusable, MouseButton,
    MouseDownEvent, MouseMoveEvent, MouseUpEvent, PathPromptOptions, Pixels, Point, ScrollDelta,
    ScrollWheelEvent, SharedString, Task, Window,
};

use crate::align;
use crate::fetch::{self, Input, Retrieval};
use crate::formats;
use crate::job::{self, Event, Outcome};
use crate::protein::Structure;
use crate::render::{self, Appearance, Camera, Coloring, Look, Palette, Repr, Scene, Species, Style};
use crate::samples::SAMPLES;
use crate::settings::{self, Settings};
use crate::theme;
use crate::ui::input::{InputEvent, TextInput};
use crate::ui::viewer;

#[cfg(feature = "simplefold")]
use crate::plugins::simplefold::{self, SimpleFold};

actions!(profold, [OpenFile]);

const PANEL_WIDTH: f32 = 372.0;
const TITLEBAR: f32 = 42.0;
/// Trajectory playback rate. Big systems play slower, at whatever pace the
/// rasterizer keeps up with.
const FRAME_INTERVAL: Duration = Duration::from_millis(60);

/// The two settings each species has a slider for.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Knob {
    Size,
    Opacity,
}

impl Knob {
    /// Largest value the slider reaches.
    fn max(self) -> f32 {
        match self {
            Knob::Size => 3.0,
            Knob::Opacity => 1.0,
        }
    }

    fn of(self, look: Look) -> f32 {
        match self {
            Knob::Size => look.size,
            Knob::Opacity => look.alpha,
        }
    }

    fn with(self, mut look: Look, v: f32) -> Look {
        match self {
            Knob::Size => look.size = v,
            Knob::Opacity => look.alpha = v,
        }
        look
    }

    fn describe(self, v: f32) -> String {
        match self {
            Knob::Size => format!("Size {v:.2}×"),
            Knob::Opacity => format!("Opacity {:.0}%", v * 100.0),
        }
    }
}

/// A value that can be typed in place of dragging for it.
#[derive(Clone, PartialEq)]
enum Field {
    Look(Species, Knob),
    /// Rotation about the screen's x (0) or y (1) axis.
    Angle(usize),
    Zoom,
    Hex(Species),
    Rgb(Species),
}

/// The colours offered in the picker, after the species' own.
const SWATCHES: [[f32; 3]; 15] = [
    [0.86, 0.24, 0.22],
    [0.95, 0.50, 0.15],
    [0.98, 0.78, 0.20],
    [0.55, 0.78, 0.25],
    [0.20, 0.62, 0.35],
    [0.15, 0.65, 0.65],
    [0.30, 0.70, 0.90],
    [0.20, 0.40, 0.85],
    [0.40, 0.30, 0.75],
    [0.65, 0.35, 0.80],
    [0.90, 0.40, 0.65],
    [0.55, 0.38, 0.25],
    [0.93, 0.93, 0.93],
    [0.55, 0.56, 0.58],
    [0.15, 0.15, 0.17],
];

#[derive(Clone, Copy, PartialEq, Eq)]
enum Drag {
    Orbit,
    Pan,
}

/// The menus in the top bar. One at a time.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Dropdown {
    Samples,
    Palette,
    Camera,
    Image,
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
    /// The last orbit step, which becomes the throw when the mouse comes up.
    last_step: [f32; 2],
    /// Yaw and pitch per frame while the view is still coasting.
    spin: Option<[f32; 2]>,
    spin_task: Option<Task<()>>,
    momentum: bool,

    open: Option<Dropdown>,
    export: Option<Task<()>>,
    /// The open dialog, while it is up.
    prompt: Option<Task<()>>,
    /// Whether the periodic cell is drawn.
    show_cell: bool,
    /// Size and opacity per element or atom type.
    appearance: Appearance,
    /// Elements in their standard colours rather than the palette's.
    fixed_colors: bool,
    /// The slider being dragged, if any.
    sliding: Option<(Species, Knob)>,
    /// Where each slider was last laid out.
    sliders: Rc<RefCell<HashMap<(Species, Knob), Bounds<Pixels>>>>,
    /// One text field, moved to whichever value is being typed.
    editor: Entity<TextInput>,
    editing: Option<Field>,
    /// The species whose colour picker is open.
    picking: Option<Species>,
    /// The settings as last written, so they are written only on change.
    saved: String,

    /// Trajectory playback.
    playing: bool,
    play_task: Option<Task<()>>,
    /// Where the scrubber was last laid out, for turning a click into a frame.
    scrubber: Rc<Shared<Bounds<Pixels>>>,
    scrubbing: bool,

    #[cfg(feature = "simplefold")]
    pub(crate) simplefold: SimpleFold,

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
        let placeholder = if cfg!(feature = "simplefold") {
            "PDB ID, UniProt accession, file path, or a sequence"
        } else {
            "PDB ID, UniProt accession, or file path"
        };
        let input = cx.new(|cx| TextInput::new(placeholder, cx));
        cx.subscribe_in(&input, window, |this, _, event, _window, cx| match event {
            InputEvent::Changed => this.reclassify(cx),
            InputEvent::Submit => this.submit(cx),
            InputEvent::Cancel => {}
        })
        .detach();

        let editor = cx.new(|cx| TextInput::new("", cx));
        cx.subscribe_in(&editor, window, |this, _, event, window, cx| match event {
            InputEvent::Submit => {
                this.commit_edit(cx);
                window.focus(&this.focus);
            }
            InputEvent::Cancel => {
                this.editing = None;
                window.focus(&this.focus);
                cx.notify();
            }
            InputEvent::Changed => {}
        })
        .detach();
        // Clicking away from a field takes what was typed, as Return does.
        let editor_focus = editor.read(cx).focus_handle(cx);
        cx.on_blur(&editor_focus, window, |this, _, cx| this.commit_edit(cx))
            .detach();

        let focus = cx.focus_handle();
        window.focus(&input.read(cx).focus_handle(cx));
        let settings = Settings::load();

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
            style: settings.style,
            palette: settings.palette,
            cache: Rc::new(RefCell::new(viewer::Cache::default())),
            drag: None,
            last_mouse: Point::default(),
            last_step: [0.0; 2],
            spin: None,
            spin_task: None,
            momentum: settings.momentum,
            open: None,
            export: None,
            prompt: None,
            show_cell: true,
            appearance: settings.appearance.clone(),
            fixed_colors: settings.fixed_colors,
            sliding: None,
            sliders: Rc::new(RefCell::new(HashMap::new())),
            editor,
            editing: None,
            picking: None,
            saved: settings.to_json(),
            playing: false,
            play_task: None,
            scrubber: Rc::new(Shared::new(Bounds::default())),
            scrubbing: false,
            #[cfg(feature = "simplefold")]
            simplefold: SimpleFold::new(),
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
        let polymer = structure.has_polymer();
        self.coloring = if structure.has_confidence {
            Coloring::Confidence
        } else if polymer {
            Coloring::Chain
        } else {
            Coloring::Element
        };
        // One protein after another keeps whatever style was picked; anything
        // else opens in the style that suits it.
        let was_polymer = self.structure.as_ref().is_some_and(|s| s.has_polymer());
        if !(polymer && was_polymer) {
            self.repr = Repr::default_for(&structure);
        }
        // A protein's crystal cell is bookkeeping; a simulation box is not.
        self.show_cell = !polymer;
        self.stop_playing();
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
        self.repr = Repr::Cartoon;
        self.stop_playing();
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
        let residues: usize = self
            .visible()
            .iter()
            .map(|s| s.residues.iter().filter(|r| r.ca.is_some()).count())
            .sum();
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
                self.fixed_colors,
                self.show_cell,
                &self.appearance,
            ));
        }
        self.scene = Arc::new(scene.stamp());
        self.persist();
        cx.notify();
    }

    fn finish(&mut self, result: anyhow::Result<Structure>, cx: &mut Context<Self>) {
        self.busy = None;
        self.progress = None;
        match result {
            Ok(s) => {
                let frames = s.frame_count();
                self.note(format!(
                    "{} · {} atoms{}",
                    s.id,
                    thousands(s.atoms.len()),
                    if frames > 1 {
                        format!(" · {frames} frames")
                    } else {
                        String::new()
                    }
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
            Input::File(path) => self.open_path(path, cx),
            Input::Pdb(id) => self.retrieve(Retrieval::Pdb(id), cx),
            Input::UniProt(acc) => self.retrieve(Retrieval::AlphaFold(acc), cx),
            Input::Sequence(seq) => self.fold(seq, cx),
        }
    }

    #[cfg(feature = "simplefold")]
    fn start(&mut self, job: simplefold::Job, cx: &mut Context<Self>) {
        for note in job.notes {
            self.note(note);
        }
        self.run_job(job.stage, job.work, cx);
    }

    #[cfg(feature = "simplefold")]
    fn fold(&mut self, sequence: String, cx: &mut Context<Self>) {
        let job = self.simplefold.fold(sequence);
        self.start(job, cx);
    }

    #[cfg(not(feature = "simplefold"))]
    fn fold(&mut self, _sequence: String, cx: &mut Context<Self>) {
        self.busy = None;
        self.progress = None;
        self.error = Some("Folding needs the SimpleFold plugin (build with --features simplefold).".into());
        cx.notify();
    }

    /// Fold what is in the box, whatever it is. A pasted sequence goes
    /// straight to the engine; an identifier has its sequences fetched first.
    #[cfg(feature = "simplefold")]
    fn submit_fold(&mut self, cx: &mut Context<Self>) {
        if self.busy.is_some() {
            return;
        }
        self.error = None;
        let job = match self.detected.clone() {
            Input::Sequence(seq) => self.simplefold.fold(seq),
            Input::Pdb(id) => self.simplefold.fold_entry(Retrieval::Pdb(id.clone()), id),
            Input::UniProt(acc) => {
                self.simplefold.fold_entry(Retrieval::AlphaFold(acc.clone()), acc)
            }
            _ => return,
        };
        self.start(job, cx);
    }

    #[cfg(feature = "simplefold")]
    fn submit_comparison(&mut self, cx: &mut Context<Self>) {
        if self.busy.is_some() {
            return;
        }
        self.error = None;
        let job = match self.detected.clone() {
            Input::Pdb(id) => self.simplefold.compare(Retrieval::Pdb(id)),
            Input::UniProt(acc) => self.simplefold.compare(Retrieval::AlphaFold(acc)),
            _ => return,
        };
        self.start(job, cx);
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
                if cfg!(feature = "simplefold") {
                    self.note("No AlphaFold model; folding the UniProt sequence instead.");
                    self.fold(seq, cx);
                } else {
                    self.finish(Err(anyhow::anyhow!("no AlphaFold model for this accession")), cx);
                }
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

    /// Read a structure or trajectory from disk.
    pub fn open_path(&mut self, path: PathBuf, cx: &mut Context<Self>) {
        if self.busy.is_some() {
            return;
        }
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        self.run_job(
            "Opening",
            move |tx| {
                let _ = tx.send(Event::Stage(format!("Reading {name}")));
                formats::load_path(&path).map(Outcome::Structure)
            },
            cx,
        );
    }

    /// Open one of the structures built into the app.
    fn open_sample(&mut self, i: usize, cx: &mut Context<Self>) {
        if self.busy.is_some() {
            return;
        }
        let Some(sample) = SAMPLES.get(i) else {
            return;
        };
        self.error = None;
        self.run_job(
            "Opening",
            move |_| {
                let mut s = formats::load_bytes(sample.file, sample.bytes.to_vec())?;
                s.title = sample.name.to_string();
                s.source = crate::protein::Source::Sample;
                Ok(Outcome::Structure(s))
            },
            cx,
        );
    }

    /// The system's open dialog.
    fn prompt_open(&mut self, cx: &mut Context<Self>) {
        let chosen = cx.prompt_for_paths(PathPromptOptions {
            files: true,
            directories: false,
            multiple: false,
            prompt: Some("Open".into()),
        });
        self.prompt = Some(cx.spawn(async move |this, cx| {
            let Ok(Ok(Some(paths))) = chosen.await else {
                return;
            };
            if let Some(path) = paths.into_iter().next() {
                let _ = this.update(cx, |this, cx| this.open_path(path, cx));
            }
        }));
    }

    // -- trajectories -----------------------------------------------------

    fn frame_count(&self) -> usize {
        match (&self.structure, &self.compare) {
            (Some(s), None) => s.frame_count(),
            _ => 1,
        }
    }

    fn set_frame(&mut self, i: usize, cx: &mut Context<Self>) {
        let Some(s) = self.structure.as_mut() else {
            return;
        };
        let i = i.min(s.frame_count().saturating_sub(1));
        if s.frame == i {
            return;
        }
        Arc::make_mut(s).set_frame(i);
        self.rebuild(cx);
    }

    fn stop_playing(&mut self) {
        self.playing = false;
        self.play_task = None;
    }

    fn toggle_play(&mut self, cx: &mut Context<Self>) {
        if self.playing {
            self.stop_playing();
            cx.notify();
            return;
        }
        let count = self.frame_count();
        if count < 2 {
            return;
        }
        // From the end, play starts over.
        if self.structure.as_ref().is_some_and(|s| s.frame + 1 >= count) {
            self.set_frame(0, cx);
        }
        self.playing = true;
        self.play_task = Some(cx.spawn(async move |this, cx| loop {
            cx.background_executor().timer(FRAME_INTERVAL).await;
            let going = this.update(cx, |this, cx| {
                let Some(s) = &this.structure else {
                    return false;
                };
                let next = s.frame + 1;
                if next >= s.frame_count() {
                    this.playing = false;
                    cx.notify();
                    return false;
                }
                this.set_frame(next, cx);
                true
            });
            if !matches!(going, Ok(true)) {
                return;
            }
        }));
        cx.notify();
    }

    /// Write the preferences out if they have changed since last time.
    fn persist(&mut self) {
        let json = Settings {
            momentum: self.momentum,
            style: self.style,
            palette: self.palette,
            fixed_colors: self.fixed_colors,
            appearance: self.appearance.clone(),
        }
        .to_json();
        if json != self.saved {
            Settings::save(&json);
            self.saved = json;
        }
    }

    /// Put the field editor on a value, with what it shows selected so
    /// typing replaces it.
    fn start_edit(&mut self, field: Field, text: String, window: &mut Window, cx: &mut Context<Self>) {
        if self.editing.as_ref() == Some(&field) {
            return;
        }
        self.commit_edit(cx);
        self.editing = Some(field);
        self.editor.update(cx, |e, cx| {
            e.set_text(text, cx);
            e.select_everything(cx);
        });
        window.focus(&self.editor.read(cx).focus_handle(cx));
        cx.notify();
    }

    /// Take what was typed. Something that does not parse is dropped, and
    /// the field goes back to showing the value as it was.
    fn commit_edit(&mut self, cx: &mut Context<Self>) {
        let Some(field) = self.editing.take() else {
            return;
        };
        let text = self.editor.read(cx).text().to_string();
        match field {
            Field::Look(species, knob) => {
                let value = match knob {
                    Knob::Size => settings::parse_size(&text),
                    Knob::Opacity => settings::parse_opacity(&text),
                };
                if let Some(v) = value {
                    let look = knob.with(self.appearance.get(&species), v);
                    self.appearance.set(species, look);
                    self.rebuild(cx);
                }
            }
            Field::Angle(axis) => {
                if let Some(deg) = settings::parse_degrees(&text) {
                    self.halt();
                    let r = deg.to_radians();
                    if axis == 0 {
                        self.camera.pitch = (Camera::HOME_PITCH + r).clamp(-1.55, 1.55);
                    } else {
                        self.camera.yaw = Camera::HOME_YAW + r;
                    }
                }
            }
            Field::Zoom => {
                if let Some(z) = settings::parse_zoom(&text) {
                    self.halt();
                    self.camera.zoom = z.clamp(0.12, 24.0);
                }
            }
            Field::Hex(species) => {
                if let Some(c) = settings::parse_hex(&text) {
                    self.set_color(species, Some(c), cx);
                }
            }
            Field::Rgb(species) => {
                if let Some(c) = settings::parse_rgb(&text) {
                    self.set_color(species, Some(c), cx);
                }
            }
        }
        cx.notify();
    }

    /// Give a species its own colour, or `None` for its element's. A colour
    /// shows when colouring by element, so that is switched to.
    fn set_color(&mut self, species: Species, color: Option<[f32; 3]>, cx: &mut Context<Self>) {
        let mut look = self.appearance.get(&species);
        look.color = color;
        self.appearance.set(species, look);
        if self.compare.is_none() {
            self.coloring = Coloring::Element;
        }
        self.rebuild(cx);
    }

    /// A value that turns into a text field when clicked.
    fn field(
        &self,
        id: impl Into<gpui::ElementId>,
        field: Field,
        text: String,
        width: f32,
        cx: &mut Context<Self>,
    ) -> gpui::AnyElement {
        if self.editing.as_ref() == Some(&field) {
            return div()
                .w(px(width))
                .h(px(18.))
                .px(px(4.))
                .flex()
                .items_center()
                .rounded(px(4.))
                .bg(theme::SURFACE)
                .border_1()
                .border_color(theme::BORDER_FOCUS)
                .font_family(theme::MONO)
                .text_size(px(10.5))
                .text_color(theme::INK)
                .on_mouse_down(MouseButton::Left, |_, _, cx| cx.stop_propagation())
                .child(self.editor.clone())
                .into_any_element();
        }
        div()
            .id(id.into())
            .w(px(width))
            .h(px(18.))
            .px(px(4.))
            .flex()
            .items_center()
            .justify_end()
            .rounded(px(4.))
            .font_family(theme::MONO)
            .text_size(px(10.5))
            .text_color(theme::MUTED)
            .cursor(gpui::CursorStyle::IBeam)
            .hover(|s| s.bg(theme::SUNKEN).text_color(theme::INK))
            .child(text.clone())
            .on_mouse_down(
                MouseButton::Left,
                cx.listener(move |this, _, window, cx| {
                    cx.stop_propagation();
                    this.start_edit(field.clone(), text.clone(), window, cx);
                }),
            )
            .into_any_element()
    }

    /// Set a species' size or opacity from where the pointer is on its
    /// slider. Values near the default snap to it, so it can be found again.
    fn slide(&mut self, species: Species, knob: Knob, x: Pixels, cx: &mut Context<Self>) {
        let Some(bounds) = self.sliders.borrow().get(&(species.clone(), knob)).copied() else {
            return;
        };
        let width = f32::from(bounds.size.width);
        if width <= 0.0 {
            return;
        }
        let t = (f32::from(x - bounds.origin.x) / width).clamp(0.0, 1.0);
        let mut v = t * knob.max();
        if (v - 1.0).abs() < knob.max() * 0.03 {
            v = 1.0;
        }
        let look = knob.with(self.appearance.get(&species), v);
        if look != self.appearance.get(&species) {
            self.appearance.set(species, look);
            self.rebuild(cx);
        }
    }

    /// The frame under the pointer on the scrubber.
    fn scrub_to(&mut self, x: Pixels, cx: &mut Context<Self>) {
        let bounds = self.scrubber.get();
        let width = f32::from(bounds.size.width);
        if width <= 0.0 {
            return;
        }
        let t = (f32::from(x - bounds.origin.x) / width).clamp(0.0, 1.0);
        let last = self.frame_count().saturating_sub(1);
        self.set_frame((t * last as f32).round() as usize, cx);
    }

    /// Write what is on screen to a PNG, at a multiple of the viewport size.
    /// The molecule only: no panel, no controls, no pointer hints.
    fn save_image(&mut self, scale: usize, cx: &mut Context<Self>) {
        let (width, height) = self.cache.borrow().size;
        if width == 0 || height == 0 {
            return;
        }
        let (width, height) = (width * scale, height * scale);
        let scene = self.scene.clone();
        let camera = self.camera;
        let style = self.style;
        let stem = self
            .structure
            .as_ref()
            .map(|s| s.id.to_lowercase())
            .unwrap_or_else(|| "profold".into());
        let suggested = format!("{stem}.png");
        let directory = dirs::home_dir()
            .map(|home| home.join("Desktop"))
            .unwrap_or_else(|| std::path::PathBuf::from("."));

        let chosen = cx.prompt_for_new_path(&directory, Some(&suggested));
        self.export = Some(cx.spawn(async move |this, cx| {
            let Ok(Ok(Some(path))) = chosen.await else {
                return;
            };
            let drawn = cx
                .background_executor()
                .spawn(async move {
                    let frame =
                        render::render(&scene, &camera, width, height, theme::VIEW, style);
                    // The rasterizer hands back BGRA; PNG wants it the other
                    // way round.
                    let mut rgba = frame.bgra;
                    for p in rgba.chunks_exact_mut(4) {
                        p.swap(0, 2);
                    }
                    image::save_buffer(
                        &path,
                        &rgba,
                        width as u32,
                        height as u32,
                        image::ColorType::Rgba8,
                    )
                    .map(|()| path)
                })
                .await;
            let _ = this.update(cx, |this, cx| {
                match drawn {
                    Ok(path) => {
                        this.note(format!("Saved {}", path.display()));
                        cx.reveal_path(&path);
                    }
                    Err(e) => this.error = Some(format!("{e}").into()),
                }
                cx.notify();
            });
        }));
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
    #[cfg(feature = "simplefold")]
    pub fn open_comparison(&mut self, text: &str, cx: &mut Context<Self>) {
        let text = text.to_string();
        self.input.update(cx, |input, cx| input.set_text(text, cx));
        self.reclassify(cx);
        self.submit_comparison(cx);
    }

    /// Fold an identifier's own sequences and show the prediction alone.
    #[cfg(feature = "simplefold")]
    pub fn open_fold(&mut self, text: &str, cx: &mut Context<Self>) {
        let text = text.to_string();
        self.input.update(cx, |input, cx| input.set_text(text, cx));
        self.reclassify(cx);
        self.submit_fold(cx);
    }

    // -- camera -----------------------------------------------------------

    /// Let go of the view and it keeps turning, slowing to a stop. The throw
    /// is whatever the last step of the drag was.
    fn coast(&mut self, velocity: [f32; 2], cx: &mut Context<Self>) {
        self.spin = Some(velocity);
        self.spin_task = Some(cx.spawn(async move |this, cx| loop {
            cx.background_executor()
                .timer(Duration::from_millis(16))
                .await;
            let turning = this.update(cx, |this, cx| {
                let Some(v) = this.spin else {
                    return false;
                };
                this.camera.yaw += v[0];
                this.camera.pitch = (this.camera.pitch + v[1]).clamp(-1.55, 1.55);
                let slowed = [v[0] * 0.94, v[1] * 0.94];
                this.spin = (slowed[0].abs() + slowed[1].abs() > 0.0007).then_some(slowed);
                cx.notify();
                this.spin.is_some()
            });
            if !matches!(turning, Ok(true)) {
                return;
            }
        }));
    }

    fn halt(&mut self) {
        self.spin = None;
        self.spin_task = None;
        self.last_step = [0.0; 2];
    }

    fn on_viewport_down(&mut self, e: &MouseDownEvent, _: &mut Window, cx: &mut Context<Self>) {
        self.halt();
        self.open = None;
        if e.click_count >= 2 {
            let all: Vec<&Structure> = self.visible().into_iter().map(|s| s.as_ref()).collect();
            if !all.is_empty() {
                self.camera = Camera::frame_all(&all);
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
        if let Some((species, knob)) = self.sliding.clone() {
            self.slide(species, knob, e.position.x, cx);
            return;
        }
        if self.scrubbing {
            self.scrub_to(e.position.x, cx);
            return;
        }
        let Some(drag) = self.drag else { return };
        let dx = f32::from(e.position.x - self.last_mouse.x);
        let dy = f32::from(e.position.y - self.last_mouse.y);
        self.last_mouse = e.position;
        match drag {
            Drag::Orbit => {
                // Negated so the molecule follows the pointer: dragging right
                // brings the near face right, rather than spinning away.
                self.last_step = [-dx * 0.008, -dy * 0.008];
                self.camera.yaw += self.last_step[0];
                self.camera.pitch = (self.camera.pitch + self.last_step[1]).clamp(-1.55, 1.55);
            }
            Drag::Pan => {
                self.camera.pan[0] += dx;
                self.camera.pan[1] += dy;
            }
        }
        cx.notify();
    }

    fn on_mouse_up(&mut self, _: &MouseUpEvent, _: &mut Window, cx: &mut Context<Self>) {
        self.scrubbing = false;
        self.sliding = None;
        let was = self.drag.take();
        if was.is_some() {
            let throw = self.last_step;
            let thrown = throw[0].abs() + throw[1].abs() > 0.004;
            if self.momentum && was == Some(Drag::Orbit) && thrown {
                self.coast(throw, cx);
            }
            self.last_step = [0.0; 2];
            cx.notify();
        }
    }

    fn on_scroll(&mut self, e: &ScrollWheelEvent, _: &mut Window, cx: &mut Context<Self>) {
        self.halt();
        let dy = match e.delta {
            ScrollDelta::Pixels(p) => f32::from(p.y),
            ScrollDelta::Lines(p) => p.y * 18.0,
        };
        self.camera.zoom = (self.camera.zoom * (1.0 + dy * 0.0022)).clamp(0.12, 24.0);
        cx.notify();
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
    ("4HHB", "haemoglobin with its haems"),
];

fn examples() -> Vec<(&'static str, &'static str)> {
    #[cfg_attr(not(feature = "simplefold"), allow(unused_mut))]
    let mut out = EXAMPLES.to_vec();
    #[cfg(feature = "simplefold")]
    out.push(simplefold::EXAMPLE);
    out
}

pub(crate) fn pill(id: &'static str, label: impl Into<SharedString>, active: bool) -> gpui::Stateful<gpui::Div> {
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

/// A pill whose id carries an index, for rows built from a list.
fn pill_keyed(
    id: (&'static str, usize),
    label: impl Into<SharedString>,
    active: bool,
) -> gpui::Stateful<gpui::Div> {
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

/// A menu title in the top bar. Filled while its menu is open.
fn tab(id: &'static str, text: &'static str, open: bool) -> gpui::Stateful<gpui::Div> {
    div()
        .id(id)
        // The bar behind is the window's drag handle; a press here is a
        // click, not the start of moving the window.
        .on_mouse_down(MouseButton::Left, swallow)
        .px(px(9.))
        .py(px(3.))
        .rounded(px(5.))
        .text_size(px(11.5))
        .cursor_pointer()
        .when(open, |d| d.bg(theme::SUNKEN).text_color(theme::INK))
        .when(!open, |d| {
            d.text_color(theme::MUTED)
                .hover(|s| s.bg(theme::SUNKEN).text_color(theme::TEXT))
        })
        .child(text)
}

/// The card a menu is drawn on, hanging under its title.
fn menu(width: f32) -> gpui::Div {
    div()
        .absolute()
        .top(px(24.))
        .left(px(-3.))
        .w(px(width))
        .p(px(3.))
        .rounded(px(8.))
        .bg(theme::SURFACE)
        .border_1()
        .border_color(theme::BORDER)
        .shadow_md()
        .flex()
        .flex_col()
        .gap(px(1.))
        .on_mouse_down(MouseButton::Left, swallow)
}

/// Keep a press on a floating card to itself. A menu hangs out over the
/// viewport, which closes menus on mouse down; without this the menu would
/// be gone before the click on its item finished.
fn swallow(_: &MouseDownEvent, _: &mut Window, cx: &mut App) {
    cx.stop_propagation();
}

/// One line of a menu: a name, a tick when it is the one in use.
fn item(id: impl Into<gpui::ElementId>, text: String, chosen: bool) -> gpui::Stateful<gpui::Div> {
    div()
        .id(id.into())
        .flex()
        .items_center()
        .justify_between()
        .gap(px(10.))
        .px(px(8.))
        .py(px(4.))
        .rounded(px(5.))
        .cursor_pointer()
        .when(!chosen, |d| d.hover(|s| s.bg(theme::SUNKEN)))
        .child(
            div()
                .text_size(px(11.5))
                .text_color(if chosen { theme::INK } else { theme::TEXT })
                .child(text),
        )
        .when(chosen, |d| {
            d.child(div().text_size(px(10.)).text_color(theme::MUTED).child("✓"))
        })
}

/// The hairline between groups of menu items.
fn divide() -> impl IntoElement {
    div().my(px(3.)).h(px(1.)).bg(theme::BORDER)
}

/// An action beside the primary one: same size, outlined rather than filled.
#[cfg_attr(not(feature = "simplefold"), allow(dead_code))]
fn secondary(id: &'static str, text: &'static str) -> gpui::Stateful<gpui::Div> {
    div()
        .id(id)
        .px(px(11.))
        .py(px(4.))
        .rounded(px(6.))
        .text_size(px(11.5))
        .border_1()
        .border_color(theme::BORDER)
        .text_color(theme::MUTED)
        .cursor_pointer()
        .hover(|s| s.bg(theme::SUNKEN).text_color(theme::TEXT))
        .child(text)
}

pub(crate) fn label(text: &'static str) -> impl IntoElement {
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

/// A palette is shown as its own ramp, beside its name in the menu.
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

/// gpui has no letter spacing, so the small capitals are spaced by hand.
fn spaced(text: &str) -> String {
    let mut out = String::with_capacity(text.len() * 2);
    for (i, c) in text.chars().enumerate() {
        if i > 0 {
            out.push('\u{2009}');
        }
        out.push(c);
    }
    out
}

/// The field name in the left column of the specification table.
fn key(text: &str) -> impl IntoElement {
    div()
        .text_size(px(9.))
        .text_color(theme::FAINT)
        .child(spaced(&text.to_uppercase()))
}

/// One line of the specification table: the field, then the figure, with a
/// hairline under it so the column reads as a single ruled block.
fn spec(field: &str, value: String, mono: bool) -> impl IntoElement {
    div()
        .flex()
        .items_baseline()
        .justify_between()
        .gap(px(10.))
        .py(px(8.))
        .border_b_1()
        .border_color(theme::BORDER)
        .child(key(field))
        .child(
            div()
                .when(mono, |d| d.font_family(theme::MONO))
                .text_size(px(11.5))
                .text_color(theme::INK)
                .child(value),
        )
}

/// 12 as ₁₂, for chemical formulae.
fn subscript(n: usize) -> String {
    n.to_string()
        .chars()
        .map(|c| char::from_u32(0x2080 + c.to_digit(10).unwrap_or(0)).unwrap_or(c))
        .collect()
}

/// 4779 as 4,779: the counts are easier to size up at a glance.
fn thousands(n: usize) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// PDB titles are one long shouted sentence describing an experiment. The
/// name of the molecule is buried in the middle of it, so it is dug out for
/// the headline and the whole sentence set as prose beneath. When nothing
/// short enough can be found the identifier takes the headline instead.
fn headline(s: &Structure) -> (String, Option<String>) {
    let title = s.title.trim();
    let id = s.id.to_uppercase();
    if title.is_empty() {
        return (id, None);
    }
    let prose = sentence_case(title);
    if prose.chars().count() <= 24 {
        return (prose, None);
    }
    match subject(&prose) {
        Some(name) => (name, Some(prose)),
        None => (id, Some(prose)),
    }
}

/// The molecule a title is about: what is left after the description of the
/// experiment in front of it and the conditions behind it are taken away.
fn subject(prose: &str) -> Option<String> {
    const LEADING: [&str; 8] = [
        "the ",
        "crystal structure analysis of ",
        "crystal structure of ",
        "x-ray structure of ",
        "solution structure of ",
        "nmr structure of ",
        "atomic structure of ",
        "structure of ",
    ];
    const TRAILING: [&str; 7] = [
        " at ",
        " in complex",
        " complexed",
        " bound to",
        " with ",
        " from ",
        ", ",
    ];

    let mut rest = prose;
    // Prefixes stack: "The crystal structure of ..." sheds both.
    let mut trimming = true;
    while trimming {
        trimming = false;
        let lower = rest.to_ascii_lowercase();
        for prefix in LEADING {
            if lower.starts_with(prefix) {
                rest = &rest[prefix.len()..];
                trimming = true;
                break;
            }
        }
    }
    let lower = rest.to_ascii_lowercase();
    let cut = TRAILING
        .iter()
        .filter_map(|marker| lower.find(marker))
        .min()
        .unwrap_or(rest.len());
    let name = rest[..cut].trim_end_matches(|c: char| !c.is_alphanumeric());

    let length = name.chars().count();
    if !(3..=26).contains(&length) || !name.chars().any(|c| c.is_alphabetic()) {
        return None;
    }
    let mut chars = name.chars();
    Some(chars.next()?.to_uppercase().collect::<String>() + chars.as_str())
}

/// Take the shout out of an all-capitals title, leaving identifiers, units and
/// single letters (chain names, the A of an angstrom) as they were.
fn sentence_case(text: &str) -> String {
    let shouted = text
        .chars()
        .filter(|c| c.is_alphabetic())
        .all(|c| c.is_uppercase());
    if !shouted {
        return text.to_string();
    }
    let mut out = String::with_capacity(text.len());
    for (i, word) in text.split_whitespace().enumerate() {
        if i > 0 {
            out.push(' ');
        }
        let keep = word.chars().any(|c| c.is_ascii_digit())
            || word.chars().filter(|c| c.is_alphabetic()).count() <= 1;
        if keep {
            out.push_str(word);
        } else {
            out.extend(word.chars().map(|c| c.to_ascii_lowercase()));
        }
    }
    let mut chars = out.chars();
    match chars.next() {
        Some(first) => first.to_uppercase().collect::<String>() + chars.as_str(),
        None => out,
    }
}

impl Profold {
    fn toggle(&mut self, which: Dropdown, cx: &mut Context<Self>) {
        self.open = (self.open != Some(which)).then_some(which);
        cx.notify();
    }

    /// The bar across the top: the name, the menus that govern the viewport,
    /// and what is loaded. It is also the window's drag handle.
    fn topbar(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let loaded = self.structure.is_some();
        let badge = self.structure.as_ref().map(|s| match s.format {
            Some(format) if s.source.is_local() => {
                format!("{} · {format}", s.id)
            }
            _ => format!("{} · {}", s.id.to_uppercase(), s.source.label()),
        });

        div()
            .id("topbar")
            .flex()
            .items_center()
            .gap(px(2.))
            .h(px(TITLEBAR))
            .pl(px(84.))
            .pr(px(16.))
            .bg(theme::SURFACE)
            .border_b_1()
            .border_color(theme::BORDER)
            .on_mouse_down(MouseButton::Left, |_, window, _| window.start_window_move())
            .child(
                div()
                    .pr(px(14.))
                    .text_size(px(13.))
                    .text_color(theme::INK)
                    .child("profold"),
            )
            .child(
                tab("open-file", "Open…", false)
                    .on_click(cx.listener(|this, _, _, cx| this.prompt_open(cx))),
            )
            .child(self.tab_with_menu(
                Dropdown::Samples,
                "samples",
                "Samples",
                self.samples_menu(cx).into_any_element(),
                cx,
            ))
            .when(loaded, |d| {
                d                .child(self.tab_with_menu(
                    Dropdown::Camera,
                    "camera",
                    "Camera",
                    self.camera_menu(cx).into_any_element(),
                    cx,
                ))
                .child(self.tab_with_menu(
                    Dropdown::Image,
                    "image",
                    "Image",
                    self.image_menu(cx).into_any_element(),
                    cx,
                ))
            })
            .child(div().flex_1())
            .when(loaded, |d| d.child(self.angle_readout(cx)))
            .when_some(badge, |d, badge| {
                d.child(
                    div()
                        .font_family(theme::MONO)
                        .text_size(px(10.))
                        .text_color(theme::FAINT)
                        .child(badge),
                )
            })
    }

    /// How far the view has turned about each axis since it opened, and how
    /// far it is zoomed, live as the molecule is moved. Click a figure to
    /// type one.
    fn angle_readout(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let [x, y] = self.camera.angles();
        let zoom = self.camera.zoom;
        let name = |n: &'static str| div().text_size(px(9.)).text_color(theme::FAINT).child(n);
        // Rounded first, so a hair either side of zero is not "-0".
        let degrees = |d: f32| format!("{}°", d.round() as i32);
        div()
            .flex()
            .items_center()
            .gap(px(3.))
            .mr(px(18.))
            .child(name("X"))
            .child(self.field("angle-x", Field::Angle(0), degrees(x), 42., cx))
            .child(div().w(px(6.)))
            .child(name("Y"))
            .child(self.field("angle-y", Field::Angle(1), degrees(y), 42., cx))
            .child(div().w(px(6.)))
            .child(name("ZOOM"))
            .child(self.field("zoom", Field::Zoom, format!("{zoom:.2}×"), 50., cx))
    }

    /// A menu title with its card hanging under it, painted over everything
    /// below so it is not clipped by the panel or the viewport.
    fn tab_with_menu(
        &self,
        which: Dropdown,
        id: &'static str,
        name: &'static str,
        card: gpui::AnyElement,
        cx: &mut Context<Self>,
    ) -> impl IntoElement {
        let open = self.open == Some(which);
        div()
            .relative()
            .child(
                tab(id, name, open)
                    .debug_selector(move || format!("tab-{id}"))
                    .on_click(cx.listener(move |this, _, _, cx| {
                    this.toggle(which, cx)
                })),
            )
            .when(open, |d| d.child(deferred(card).with_priority(2)))
    }

    /// The structures built into the app, to open at any time.
    fn samples_menu(&self, cx: &mut Context<Self>) -> impl IntoElement {
        menu(260.)
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.open = None;
                cx.notify();
            }))
            .children(SAMPLES.iter().enumerate().map(|(i, sample)| {
                div()
                    .id(("sample-menu", i))
                    .flex()
                    .items_baseline()
                    .justify_between()
                    .gap(px(10.))
                    .px(px(8.))
                    .py(px(4.))
                    .rounded(px(5.))
                    .cursor_pointer()
                    .hover(|s| s.bg(theme::SUNKEN))
                    .child(
                        div()
                            .text_size(px(11.5))
                            .text_color(theme::TEXT)
                            .child(sample.name),
                    )
                    .child(
                        div()
                            .text_size(px(9.5))
                            .text_color(theme::FAINT)
                            .child(sample.file),
                    )
                    .on_click(cx.listener(move |this, _, _, cx| {
                        this.open = None;
                        this.open_sample(i, cx);
                    }))
            }))
    }

    /// The controls that govern the drawing, floating over the viewport:
    /// what is drawn and how it is coloured along the top; palette, style
    /// and the cell along the foot.
    fn controls(&self, cx: &mut Context<Self>) -> Vec<gpui::AnyElement> {
        let repr = self.repr;
        let coloring = self.coloring;
        let style = self.style;
        let structure = self.structure.as_ref();
        let has_confidence = structure.is_some_and(|s| s.has_confidence);
        let polymer = structure.is_some_and(|s| s.has_polymer());
        let has_cell = structure.is_some_and(|s| s.cell.is_some());
        let options: &[Coloring] = if self.compare.is_some() {
            &Coloring::COMPARE
        } else {
            &Coloring::ALL
        };
        // Clear of the trajectory bar when there is one.
        let bottom = if self.frame_count() > 1 { 62. } else { 20. };

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

        // A row of groups across the viewport, centred. Clicks on it are
        // its own: the viewport underneath would take them as a drag.
        let bar = || {
            div()
                .absolute()
                .left(px(16.))
                .right(px(16.))
                .flex()
                .flex_wrap()
                .justify_center()
                .gap(px(8.))
                .on_mouse_down(MouseButton::Left, |_, _, cx| cx.stop_propagation())
        };

        let top = bar()
            .top(px(12.))
            .child(group(
                Repr::ALL
                    .iter()
                    .map(|r| {
                        let r = *r;
                        // A backbone style means nothing without a backbone.
                        let enabled = !r.needs_polymer() || polymer;
                        pill_keyed(("repr", r as usize), r.label(), r == repr)
                            .when(!enabled, |d| d.text_color(theme::FAINT))
                            .when(enabled, |d| d.cursor_pointer())
                            .on_click(cx.listener(move |this, _, _, cx| {
                                if enabled {
                                    this.repr = r;
                                    this.rebuild(cx);
                                }
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
                        // pLDDT means nothing without a predicted model.
                        let enabled = c != Coloring::Confidence || has_confidence;
                        let name = match c {
                            Coloring::Chain if !polymer => "Molecule",
                            _ => c.label(),
                        };
                        pill_keyed(("colour", c as usize), name, c == coloring)
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
            ));

        let bottom = bar()
            .bottom(px(bottom))
            .child(self.palette_button(cx))
            .child(group(
                Style::ALL
                    .iter()
                    .map(|st| {
                        let st = *st;
                        pill_keyed(("style", st as usize), st.label(), st == style)
                            .cursor_pointer()
                            .on_click(cx.listener(move |this, _, _, cx| {
                                this.style = st;
                                this.rebuild(cx);
                            }))
                            .into_any_element()
                    })
                    .collect(),
            ))
            .when(has_cell, |d| {
                d.child(group(vec![pill("cell", "Cell", self.show_cell)
                    .cursor_pointer()
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.show_cell = !this.show_cell;
                        this.rebuild(cx);
                    }))
                    .into_any_element()]))
            });

        vec![top.into_any_element(), bottom.into_any_element()]
    }

    /// The current palette, shown as its own ramp; its menu opens upward,
    /// since the button sits at the foot of the viewport.
    fn palette_button(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let current = self.palette;
        let open = self.open == Some(Dropdown::Palette);
        let list = menu(196.)
            .top_auto()
            .bottom(px(34.))
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.open = None;
                cx.notify();
            }))
            .children(Palette::ALL.iter().map(|p| {
                let p = *p;
                let chosen = p == current;
                div()
                    .id(("palette", p as usize))
                    .debug_selector(move || format!("palette-{}", p.label()))
                    .flex()
                    .items_center()
                    .gap(px(9.))
                    .px(px(8.))
                    .py(px(4.))
                    .rounded(px(5.))
                    .cursor_pointer()
                    .when(!chosen, |d| d.hover(|s| s.bg(theme::SUNKEN)))
                    .child(swatch_strip(p))
                    .child(
                        div()
                            .flex_1()
                            .text_size(px(11.5))
                            .text_color(if chosen { theme::INK } else { theme::TEXT })
                            .child(p.label()),
                    )
                    .when(chosen, |d| {
                        d.child(div().text_size(px(10.)).text_color(theme::MUTED).child("✓"))
                    })
                    .on_click(cx.listener(move |this, _, _, cx| {
                        this.palette = p;
                        this.open = None;
                        this.rebuild(cx);
                    }))
            }));
        // Fixed colours leave the palette nothing to do, so it says so by
        // going pale and staying shut.
        let fixed = self.fixed_colors;

        div()
            .relative()
            .child(
                div()
                    .id("palette-button")
                    .debug_selector(|| "palette-button".into())
                    .flex()
                    .items_center()
                    .gap(px(7.))
                    .p(px(3.))
                    .pr(px(9.))
                    .rounded(px(8.))
                    .bg(theme::SURFACE)
                    .border_1()
                    .border_color(theme::BORDER)
                    .shadow_sm()
                    .when(fixed, |d| d.opacity(0.4))
                    .when(!fixed, |d| d.cursor_pointer().hover(|s| s.bg(theme::SUNKEN)))
                    .tooltip(move |_, cx| {
                        cx.new(|_| {
                            Tip(if fixed {
                                "Colours are fixed; turn off Fixed colours to use a palette".into()
                            } else {
                                "Palette".into()
                            })
                        })
                        .into()
                    })
                    .child(div().p(px(3.)).child(swatch_strip(current)))
                    .child(
                        div()
                            .text_size(px(11.))
                            .text_color(theme::TEXT)
                            .child(current.label()),
                    )
                    .child(
                        div()
                            .text_size(px(12.))
                            .text_color(theme::MUTED)
                            .child(if open { "▾" } else { "▴" }),
                    )
                    .on_click(cx.listener(|this, _, _, cx| {
                        if !this.fixed_colors {
                            this.toggle(Dropdown::Palette, cx)
                        }
                    })),
            )
            .when(open && !fixed, |d| d.child(deferred(list).with_priority(2)))
    }

    fn camera_menu(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let momentum = self.momentum;
        menu(184.)
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.open = None;
                cx.notify();
            }))
            .child(
                item("momentum", "Momentum".into(), momentum).on_click(cx.listener(
                    |this, _, _, cx| {
                        this.momentum = !this.momentum;
                        if !this.momentum {
                            this.halt();
                        }
                        this.persist();
                        cx.notify();
                    },
                )),
            )
            .child(divide())
            .child(
                item("recentre", "Centre".into(), false).on_click(cx.listener(|this, _, _, cx| {
                    this.camera.pan = [0.0, 0.0];
                    this.open = None;
                    cx.notify();
                })),
            )
            .child(
                item("reset-view", "Reset view".into(), false)
                    .debug_selector(|| "reset-view".into())
                    .on_click(cx.listener(
                    |this, _, _, cx| {
                        this.halt();
                        // Framed before the camera is touched, so the borrow
                        // of the structures is over by then.
                        let framed = {
                            let all: Vec<&Structure> =
                                this.visible().into_iter().map(|s| s.as_ref()).collect();
                            (!all.is_empty()).then(|| Camera::frame_all(&all))
                        };
                        if let Some(camera) = framed {
                            this.camera = camera;
                        }
                        this.open = None;
                        cx.notify();
                    },
                )),
            )
    }

    /// Write the viewport out as a picture, at a multiple of its size.
    fn image_menu(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let (w, h) = self.cache.borrow().size;
        menu(200.)
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.open = None;
                cx.notify();
            }))
            .children([1usize, 2, 3].into_iter().map(move |scale| {
                let size = format!("{} × {}", w * scale, h * scale);
                div()
                    .id(("image", scale))
                    .flex()
                    .items_baseline()
                    .justify_between()
                    .gap(px(10.))
                    .px(px(8.))
                    .py(px(4.))
                    .rounded(px(5.))
                    .cursor_pointer()
                    .hover(|s| s.bg(theme::SUNKEN))
                    .child(
                        div()
                            .text_size(px(11.5))
                            .text_color(theme::TEXT)
                            .child(format!("Save PNG · {scale}×")),
                    )
                    .child(
                        div()
                            .font_family(theme::MONO)
                            .text_size(px(9.5))
                            .text_color(theme::FAINT)
                            .child(size),
                    )
                    .on_click(cx.listener(move |this, _, _, cx| {
                        this.open = None;
                        this.save_image(scale, cx);
                    }))
            }))
    }

    fn search(&self, focused: bool, cx: &mut Context<Self>) -> impl IntoElement {
        let detected = self.detected.clone();
        let ready = !matches!(detected, Input::Empty | Input::Unknown)
            && (cfg!(feature = "simplefold") || !matches!(detected, Input::Sequence(_)));
        // An identifier can be fetched, folded from its own sequences, or both
        // at once. A pasted sequence can only be folded, and Fold is then the
        // primary action rather than one of three.
        let identifier = matches!(detected, Input::Pdb(_) | Input::UniProt(_));
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
                            .when(identifier && !busy, |d| self.plugin_actions(d, cx))
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

    /// Buttons a plugin puts beside the primary action for an identifier.
    #[cfg(feature = "simplefold")]
    fn plugin_actions(&self, d: gpui::Div, cx: &mut Context<Self>) -> gpui::Div {
        d.child(
            secondary("fold", "Fold").on_click(cx.listener(|this, _, _, cx| this.submit_fold(cx))),
        )
        .child(
            secondary("compare", "Compare")
                .on_click(cx.listener(|this, _, _, cx| this.submit_comparison(cx))),
        )
    }

    #[cfg(feature = "simplefold")]
    fn plugin_options(&self, cx: &mut Context<Self>) -> Option<gpui::AnyElement> {
        Some(simplefold::options(&self.simplefold, cx).into_any_element())
    }

    #[cfg(not(feature = "simplefold"))]
    fn plugin_options(&self, _: &mut Context<Self>) -> Option<gpui::AnyElement> {
        None
    }

    #[cfg(not(feature = "simplefold"))]
    fn plugin_actions(&self, d: gpui::Div, _: &mut Context<Self>) -> gpui::Div {
        d
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
            .child(label("SAMPLES"))
            .children(SAMPLES.iter().enumerate().map(|(i, sample)| {
                div()
                    .id(("sample", i))
                    .flex()
                    .items_baseline()
                    .gap(px(8.))
                    .cursor_pointer()
                    .child(
                        div()
                            .text_size(px(11.5))
                            .text_color(theme::TEXT)
                            .hover(|s| s.text_color(theme::INK))
                            .child(sample.name),
                    )
                    .child(
                        div()
                            .text_size(px(10.5))
                            .text_color(theme::FAINT)
                            .child(sample.note),
                    )
                    .on_click(cx.listener(move |this, _, _, cx| this.open_sample(i, cx)))
            }))
            .child(div().h(px(8.)))
            .child(label("FETCH"))
            .children(examples().into_iter().enumerate().map(|(i, (value, note))| {
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
                            .child(note),
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
                // Under the representation and colouring bar.
                .top(px(56.))
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
                                format!("{} {}", prediction.source.label(), prediction.id),
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

    /// One row per chain, in the colour that chain is painted. The dots go
    /// hollow when the model is coloured by something other than chain, so the
    /// key never claims a colour the viewport is not using.
    fn legend(&self, structure: &Structure) -> impl IntoElement {
        let by_chain = self.coloring == Coloring::Chain;
        let palette = self.palette;
        let total = structure.chains.len();
        let shown = total.min(6);
        let rest = total - shown;

        div()
            .pt(px(14.))
            .flex()
            .flex_col()
            .gap(px(7.))
            .children(
                structure
                    .chains
                    .iter()
                    .take(shown)
                    .enumerate()
                    .map(move |(i, chain)| {
                        let name = if chain.id.trim().is_empty() {
                            format!("Chain {}", i + 1)
                        } else {
                            format!("Chain {}", chain.id)
                        };
                        div()
                            .flex()
                            .items_center()
                            .gap(px(9.))
                            .child(
                                div()
                                    .size(px(7.))
                                    .rounded_full()
                                    .when(by_chain, |d| d.bg(swatch(palette.chain(i))))
                                    .when(!by_chain, |d| {
                                        d.border_1().border_color(theme::FAINT)
                                    }),
                            )
                            .child(
                                div()
                                    .flex_1()
                                    .text_size(px(11.5))
                                    .text_color(theme::TEXT)
                                    .child(name),
                            )
                            .child(
                                div()
                                    .font_family(theme::MONO)
                                    .text_size(px(10.))
                                    .text_color(theme::FAINT)
                                    .child(format!("{} aa", chain.residues.len())),
                            )
                    }),
            )
            .when(rest > 0, |d| {
                d.child(
                    div()
                        .pl(px(16.))
                        .text_size(px(10.5))
                        .text_color(theme::FAINT)
                        .child(format!("and {rest} more")),
                )
            })
    }

    /// One row per element — or per type, for atoms of no element — with
    /// its colour, its size and opacity as sliders and as figures to type
    /// over, and its count. The colour opens a picker.
    fn species_table(&self, structure: &Structure, cx: &mut Context<Self>) -> impl IntoElement {
        let by_element = self.coloring == Coloring::Element;
        let palette = self.palette;
        let colors = render::element_colors(structure, palette, self.fixed_colors);
        let mut rows: Vec<(Species, String, [f32; 3], usize)> = structure
            .composition()
            .into_iter()
            .filter(|(e, _)| e.is_known())
            .map(|(e, n)| {
                (Species::Element(e.0), e.symbol().to_string(), colors[e.0 as usize], n)
            })
            .collect();
        let mut by_kind = vec![0usize; structure.kinds.len()];
        for a in structure.atoms.iter().filter(|a| !a.element.is_known()) {
            if let Some(n) = by_kind.get_mut(a.kind as usize) {
                *n += 1;
            }
        }
        rows.extend(by_kind.iter().enumerate().filter(|(_, n)| **n > 0).map(|(k, n)| {
            let name = structure.kinds[k].clone();
            (Species::Kind(name.clone()), format!("Type {name}"), palette.chain(k), *n)
        }));
        let changed = !self.appearance.is_default();

        let header = div()
            .flex()
            .items_center()
            .gap(px(8.))
            .child(div().flex_1().child(key("atoms")))
            .child(div().w(px(96.)).child(key("size")))
            .child(div().w(px(96.)).child(key("opacity")))
            .child(
                div()
                    .id("reset-looks")
                    .w(px(36.))
                    .flex()
                    .justify_end()
                    .text_size(px(9.5))
                    .when(changed, |d| {
                        d.text_color(theme::MUTED)
                            .cursor_pointer()
                            .hover(|s| s.text_color(theme::INK))
                            .child("Reset")
                            .on_click(cx.listener(|this, _, _, cx| {
                                this.appearance.reset();
                                this.rebuild(cx);
                            }))
                    }),
            );

        div()
            .pt(px(14.))
            .flex()
            .flex_col()
            .gap(px(5.))
            .child(header)
            .children(rows.into_iter().enumerate().map(|(i, (species, name, default, n))| {
                let look = self.appearance.get(&species);
                let shown = look.color.unwrap_or(default);
                let picking = self.picking.as_ref() == Some(&species);
                let pick = species.clone();
                div()
                    .relative()
                    .flex()
                    .items_center()
                    .gap(px(8.))
                    .child(
                        div()
                            .flex_1()
                            .flex()
                            .items_center()
                            .gap(px(8.))
                            .overflow_hidden()
                            .when(look.hidden(), |d| d.opacity(0.45))
                            .child(
                                div()
                                    .id(("colour-dot", i))
                                    .debug_selector({
                                        let name = name.clone();
                                        move || format!("dot-{name}")
                                    })
                                    .flex_none()
                                    .size(px(13.))
                                    .flex()
                                    .items_center()
                                    .justify_center()
                                    .rounded_full()
                                    .cursor_pointer()
                                    .when(picking, |d| d.border_1().border_color(theme::INK))
                                    .hover(|d| d.bg(theme::SUNKEN))
                                    .child(
                                        div()
                                            .size(px(8.))
                                            .rounded_full()
                                            .when(by_element, |d| {
                                                d.bg(swatch(shown)).opacity(look.alpha.max(0.2))
                                            })
                                            .when(!by_element, |d| {
                                                d.border_1().border_color(theme::FAINT)
                                            }),
                                    )
                                    .tooltip(|_, cx| {
                                        cx.new(|_| Tip("Choose a colour".into())).into()
                                    })
                                    .on_click(cx.listener(move |this, _, _, cx| {
                                        this.picking = Some(pick.clone());
                                        cx.notify();
                                    })),
                            )
                            .child(
                                div()
                                    .text_size(px(11.5))
                                    .text_color(theme::TEXT)
                                    .whitespace_nowrap()
                                    .child(name.clone()),
                            ),
                    )
                    .child(
                        div()
                            .w(px(96.))
                            .flex()
                            .items_center()
                            .gap(px(4.))
                            .child(self.slider(i, species.clone(), Knob::Size, look, cx))
                            .child(self.field(
                                ("size-value", i),
                                Field::Look(species.clone(), Knob::Size),
                                format!("{:.2}×", look.size),
                                40.,
                                cx,
                            )),
                    )
                    .child(
                        div()
                            .w(px(96.))
                            .flex()
                            .items_center()
                            .gap(px(4.))
                            .child(self.slider(i, species.clone(), Knob::Opacity, look, cx))
                            .child(self.field(
                                ("opacity-value", i),
                                Field::Look(species.clone(), Knob::Opacity),
                                format!("{:.0}%", look.alpha * 100.0),
                                40.,
                                cx,
                            )),
                    )
                    .child(
                        div()
                            .w(px(36.))
                            .flex()
                            .justify_end()
                            .font_family(theme::MONO)
                            .text_size(px(10.))
                            .text_color(theme::FAINT)
                            .child(thousands(n)),
                    )
                    .when(picking, |d| {
                        d.child(
                            deferred(self.color_picker(species.clone(), name, default, look, cx))
                                .with_priority(3),
                        )
                    })
            }))
            .child(self.fixed_colors_toggle(cx))
    }

    /// Whether elements wear the palette's colours or the standard ones.
    fn fixed_colors_toggle(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let on = self.fixed_colors;
        div()
            .id("fixed-colours")
            .debug_selector(|| "fixed-colours".into())
            .mt(px(4.))
            .flex()
            .items_center()
            .gap(px(8.))
            .cursor_pointer()
            .child(
                // A switch: a track with the knob at whichever end.
                div()
                    .w(px(22.))
                    .h(px(12.))
                    .p(px(2.))
                    .rounded_full()
                    .flex()
                    .when(on, |d| d.bg(theme::INK).justify_end())
                    .when(!on, |d| d.bg(theme::BORDER))
                    .child(div().size(px(8.)).rounded_full().bg(theme::SURFACE)),
            )
            .child(
                div()
                    .text_size(px(11.))
                    .text_color(if on { theme::INK } else { theme::MUTED })
                    .child("Fixed colours"),
            )
            .child(
                div()
                    .text_size(px(9.5))
                    .text_color(theme::FAINT)
                    .child(if on {
                        "standard element colours"
                    } else {
                        "from the palette"
                    }),
            )
            .on_click(cx.listener(|this, _, _, cx| {
                this.fixed_colors = !this.fixed_colors;
                if this.fixed_colors && this.open == Some(Dropdown::Palette) {
                    this.open = None;
                }
                this.rebuild(cx);
            }))
    }

    /// Presets, then the colour as hex and as RGB to type over. The first
    /// preset is the species' own colour, which clears any choice.
    fn color_picker(
        &self,
        species: Species,
        name: String,
        default: [f32; 3],
        look: Look,
        cx: &mut Context<Self>,
    ) -> impl IntoElement {
        let current = look.color.unwrap_or(default);
        let swatch_button = |id: usize, c: [f32; 3], chosen: bool, set: Option<[f32; 3]>, species: Species| {
            div()
                .id(("swatch", id))
                .debug_selector(move || format!("swatch-{id}"))
                .size(px(20.))
                .rounded(px(5.))
                .p(px(2.))
                .cursor_pointer()
                .when(chosen, |d| d.border_1().border_color(theme::INK))
                .when(!chosen, |d| d.hover(|s| s.bg(theme::SUNKEN)))
                .child(div().size_full().rounded(px(3.)).bg(swatch(c)))
                .on_click(cx.listener(move |this, _, _, cx| {
                    this.set_color(species.clone(), set, cx)
                }))
        };
        let label = |t: &'static str| {
            div()
                .w(px(30.))
                .text_size(px(9.))
                .text_color(theme::FAINT)
                .child(t)
        };

        div()
            .absolute()
            .top(px(20.))
            .left(px(-4.))
            .w(px(232.))
            .p(px(10.))
            .rounded(px(8.))
            .bg(theme::SURFACE)
            .border_1()
            .border_color(theme::BORDER)
            .shadow_md()
            .flex()
            .flex_col()
            .gap(px(9.))
            .on_mouse_down(MouseButton::Left, swallow)
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.picking = None;
                cx.notify();
            }))
            .child(
                div()
                    .flex()
                    .items_center()
                    .justify_between()
                    .child(
                        div()
                            .text_size(px(11.5))
                            .text_color(theme::INK)
                            .child(format!("Colour of {name}")),
                    )
                    .child(
                        div()
                            .id("picker-done")
                            .text_size(px(10.5))
                            .text_color(theme::MUTED)
                            .cursor_pointer()
                            .hover(|s| s.text_color(theme::INK))
                            .child("Done")
                            .on_click(cx.listener(|this, _, _, cx| {
                                this.picking = None;
                                cx.notify();
                            })),
                    ),
            )
            .child(
                div()
                    .flex()
                    .flex_wrap()
                    .gap(px(5.))
                    .child(
                        swatch_button(0, default, look.color.is_none(), None, species.clone())
                            .tooltip(|_, cx| cx.new(|_| Tip("Its own colour".into())).into()),
                    )
                    .children(SWATCHES.iter().enumerate().map(|(k, c)| {
                        swatch_button(k + 1, *c, look.color == Some(*c), Some(*c), species.clone())
                    })),
            )
            .child(
                div()
                    .flex()
                    .items_center()
                    .gap(px(6.))
                    .child(label("HEX"))
                    .child(
                        div()
                            .size(px(14.))
                            .rounded(px(3.))
                            .bg(swatch(current))
                            .border_1()
                            .border_color(theme::BORDER),
                    )
                    .child(self.field(
                        "picker-hex",
                        Field::Hex(species.clone()),
                        settings::hex(current),
                        150.,
                        cx,
                    )),
            )
            .child(
                div()
                    .flex()
                    .items_center()
                    .gap(px(6.))
                    .child(label("RGB"))
                    .child(div().size(px(14.)))
                    .child(self.field(
                        "picker-rgb",
                        Field::Rgb(species),
                        settings::rgb(current),
                        150.,
                        cx,
                    )),
            )
    }

    /// A small horizontal slider. A tick marks the default; double-click
    /// goes back to it.
    fn slider(
        &self,
        row: usize,
        species: Species,
        knob: Knob,
        look: Look,
        cx: &mut Context<Self>,
    ) -> impl IntoElement {
        let value = knob.of(look);
        let t = (value / knob.max()).clamp(0.0, 1.0);
        let home = 1.0 / knob.max();
        let store = self.sliders.clone();
        let slot = (species.clone(), knob);
        let tip: SharedString = knob.describe(value).into();
        div()
            .id((if knob == Knob::Size { "size" } else { "opacity" }, row))
            .relative()
            .w(px(52.))
            .h(px(14.))
            .flex()
            .items_center()
            .cursor_pointer()
            .tooltip(move |_, cx| cx.new(|_| Tip(tip.clone())).into())
            .child(
                canvas(
                    move |bounds, _, _| {
                        store.borrow_mut().insert(slot.clone(), bounds);
                    },
                    |_, _, _, _| {},
                )
                .absolute()
                .size_full(),
            )
            .child(
                div()
                    .relative()
                    .w_full()
                    .h(px(3.))
                    .rounded_full()
                    .bg(theme::SUNKEN)
                    .child(div().h_full().w(relative(t)).rounded_full().bg(theme::MUTED))
                    // The default, so it can be found by eye.
                    .when(knob == Knob::Size, |d| {
                        d.child(
                            div()
                                .absolute()
                                .top(px(-2.))
                                .left(relative(home))
                                .w(px(1.))
                                .h(px(7.))
                                .bg(theme::FAINT),
                        )
                    }),
            )
            .child(
                div()
                    .absolute()
                    .top(px(3.))
                    .left(relative(t))
                    .ml(px(-4.))
                    .size(px(8.))
                    .rounded_full()
                    .bg(theme::INK),
            )
            .on_mouse_down(
                MouseButton::Left,
                cx.listener(move |this, e: &MouseDownEvent, _, cx| {
                    cx.stop_propagation();
                    if e.click_count >= 2 {
                        let look = knob.with(this.appearance.get(&species), 1.0);
                        this.appearance.set(species.clone(), look);
                        this.rebuild(cx);
                        return;
                    }
                    this.sliding = Some((species.clone(), knob));
                    this.slide(species.clone(), knob, e.position.x, cx);
                }),
            )
    }

    /// The info panel for anything that is not a protein: a molecule, a
    /// crystal, a simulation box. Atoms, bonds, the cell and the frames are
    /// what there is to say.
    fn info_atoms(&self, structure: &Structure, cx: &mut Context<Self>) -> impl IntoElement {
        let name = if structure.title.trim().is_empty() {
            structure.id.clone()
        } else {
            structure.title.trim().to_string()
        };
        let (name, blurb) = if name.chars().count() > 28 {
            (structure.id.clone(), Some(name))
        } else {
            (name, None)
        };
        let mut composition = structure.composition();
        // Hill order: carbon, then hydrogen, then the rest alphabetically;
        // without carbon, everything alphabetically.
        let has_carbon = composition.iter().any(|(e, _)| e.symbol() == "C");
        composition.sort_by_key(|(e, _)| {
            let rank = match e.symbol() {
                "C" if has_carbon => 0,
                "H" if has_carbon => 1,
                _ => 2,
            };
            (rank, e.symbol())
        });
        let formula = (composition.iter().all(|(e, _)| e.is_known()) && composition.len() <= 8)
            .then(|| {
                composition
                    .iter()
                    .map(|(e, n)| {
                        if *n == 1 {
                            e.symbol().to_string()
                        } else {
                            format!("{}{}", e.symbol(), subscript(*n))
                        }
                    })
                    .collect::<Vec<_>>()
                    .join("")
            })
            .filter(|f| f.chars().count() <= 30);
        let mass = structure.atomic_mass();
        let category = match structure.format {
            Some(format) => format!("{} / {}", structure.source.label(), format),
            None => structure.source.label().to_string(),
        };

        div()
            .flex()
            .flex_col()
            .px(px(20.))
            .pb(px(20.))
            .child(rule())
            .child(
                div()
                    .pt(px(15.))
                    .pb(px(7.))
                    .text_size(px(9.))
                    .text_color(theme::FAINT)
                    .child(spaced(&category.to_uppercase())),
            )
            .child(
                div()
                    .font_family(theme::SERIF)
                    .text_size(px(27.))
                    .line_height(px(31.))
                    .text_color(theme::INK)
                    .overflow_hidden()
                    .child(name),
            )
            .when_some(blurb, |d, blurb| {
                d.child(
                    div()
                        .pt(px(10.))
                        .text_size(px(11.5))
                        .line_height(px(18.))
                        .text_color(theme::TEXT)
                        .child(blurb),
                )
            })
            .child(self.species_table(structure, cx))
            .child(
                div()
                    .pt(px(16.))
                    .flex()
                    .flex_col()
                    .when_some(formula, |d, f| d.child(spec("formula", f, false)))
                    .child(spec("atoms", thousands(structure.atoms.len()), false))
                    .child(spec("bonds", thousands(structure.bonds.len()), false))
                    .when(structure.residues.len() > 1, |d| {
                        d.child(spec(
                            "molecules",
                            thousands(structure.residues.len()),
                            false,
                        ))
                    })
                    .when(mass > 0.0 && mass < 1.0e6, |d| {
                        d.child(spec("mass", format!("{mass:.2} g/mol"), false))
                    })
                    .when_some(structure.cell, |d, cell| {
                        let [a, b, c] = cell.lengths();
                        d.child(spec("cell", format!("{a:.2} × {b:.2} × {c:.2} Å"), false))
                    })
                    .when(structure.frame_count() > 1, |d| {
                        d.child(spec("frames", thousands(structure.frame_count()), false))
                    }),
            )
    }

    /// What is on screen, set like the opening of an article: the category in
    /// small capitals, the name large, the sentence beneath it, then the
    /// figures in a ruled table that reads straight down the left edge.
    fn info(&self, structure: &Structure, cx: &mut Context<Self>) -> gpui::AnyElement {
        if !structure.has_polymer() {
            return self.info_atoms(structure, cx).into_any_element();
        }
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
        let (name, blurb) = headline(structure);
        // The numbering is worth showing only when the file does not start at
        // one, which is the sign of a fragment or a trimmed construct.
        let residues = match (structure.residues.first(), structure.residues.last()) {
            (Some(a), Some(b)) if a.seq != 1 => format!(
                "{} aa  ·  {}–{}",
                thousands(structure.residues.len()),
                a.seq,
                b.seq
            ),
            _ => format!("{} aa", thousands(structure.residues.len())),
        };
        let chains = if structure.chains.len() <= 4 {
            let ids = structure
                .chains
                .iter()
                .map(|c| c.id.clone())
                .collect::<Vec<_>>()
                .join(" ");
            format!("{}  ·  {}", structure.chains.len(), ids)
        } else {
            structure.chains.len().to_string()
        };

        div()
            .flex()
            .flex_col()
            .px(px(20.))
            .pb(px(20.))
            .child(rule())
            .child(
                div()
                    .pt(px(15.))
                    .pb(px(7.))
                    .text_size(px(9.))
                    .text_color(theme::FAINT)
                    .child(spaced(&format!(
                        "{} / {}",
                        structure.source.label().to_uppercase(),
                        // A file says what format it was; a download, whether
                        // it was measured or predicted.
                        match structure.format {
                            Some(format) if structure.source.is_local() => format,
                            _ => structure.source.kind(),
                        }
                        .to_uppercase()
                    ))),
            )
            .child(
                div()
                    .font_family(theme::SERIF)
                    .text_size(px(27.))
                    .line_height(px(31.))
                    .text_color(theme::INK)
                    .child(name),
            )
            .when_some(blurb, |d, blurb| {
                d.child(
                    div()
                        .pt(px(10.))
                        .text_size(px(11.5))
                        .line_height(px(18.))
                        .text_color(theme::TEXT)
                        .child(blurb),
                )
            })
            .child(self.legend(structure))
            .child(self.species_table(structure, cx))
            .child(
                div()
                    .pt(px(16.))
                    .flex()
                    .flex_col()
                    .child(spec(
                        "structure",
                        structure.id.to_uppercase(),
                        true,
                    ))
                    .child(spec("residues", residues, false))
                    .child(spec("chains", chains, false))
                    .child(spec(
                        "atoms",
                        thousands(structure.atoms.len()),
                        false,
                    ))
                    .child(spec(
                        "mass",
                        format!("{:.1} kDa", structure.mass_kda()),
                        false,
                    ))
                    .when_some(structure.mean_confidence(), |d, c| {
                        d.child(spec("pLDDT", format!("{c:.0}"), false))
                    })
                    .child(
                        div()
                            .flex()
                            .items_baseline()
                            .justify_between()
                            .gap(px(10.))
                            .py(px(8.))
                            .border_b_1()
                            .border_color(theme::BORDER)
                            .child(key("sequence"))
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
                    ),
            )
            .child(
                div()
                    .pt(px(9.))
                    .overflow_hidden()
                    .whitespace_nowrap()
                    .font_family(theme::MONO)
                    .text_size(px(10.))
                    .text_color(theme::FAINT)
                    .child(preview),
            )
            .when(structure.frame_count() > 1, |d| {
                d.child(
                    div()
                        .pt(px(8.))
                        .text_size(px(10.5))
                        .text_color(theme::FAINT)
                        .child(format!("{} frames", thousands(structure.frame_count()))),
                )
            })
            .into_any_element()
    }

    fn stage(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let has_structure = self.structure.is_some();
        let busy: SharedString = self.busy.clone().unwrap_or_default();

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
                                .child("Drop a structure file here, open one with ⌘O,"),
                        )
                        .child(
                            div()
                                .text_size(px(10.5))
                                .text_color(theme::FAINT)
                                .child("or paste an identifier to fetch it"),
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
            .when(has_structure && self.busy.is_none() && !self.visible().is_empty(), |d| {
                d.children(self.controls(cx))
            })
            .when(has_structure && self.busy.is_none() && self.frame_count() > 1, |d| {
                d.child(self.transport(cx))
            })
            .when(has_structure && self.busy.is_none() && self.frame_count() <= 1, |d| {
                d.child(
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

impl Profold {
    /// Play, pause and a scrubber for a trajectory, along the foot of the
    /// viewport.
    fn transport(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let count = self.frame_count();
        let frame = self.structure.as_ref().map(|s| s.frame).unwrap_or(0);
        let t = frame as f32 / (count - 1).max(1) as f32;
        let store = self.scrubber.clone();

        div()
            .absolute()
            .bottom(px(12.))
            .left(px(20.))
            .right(px(20.))
            .flex()
            .items_center()
            .gap(px(10.))
            .px(px(10.))
            .py(px(6.))
            .rounded(px(8.))
            .bg(theme::SURFACE)
            .border_1()
            .border_color(theme::BORDER)
            .shadow_sm()
            // The viewport underneath would take this as the start of a drag.
            .on_mouse_down(MouseButton::Left, |_, _, cx| cx.stop_propagation())
            .child(
                pill("step-back", "‹", false)
                    .cursor_pointer()
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.stop_playing();
                        let f = this.structure.as_ref().map(|s| s.frame).unwrap_or(0);
                        this.set_frame(f.saturating_sub(1), cx);
                    })),
            )
            .child(
                pill("play", if self.playing { "Pause" } else { "Play" }, self.playing)
                    .cursor_pointer()
                    .on_click(cx.listener(|this, _, _, cx| this.toggle_play(cx))),
            )
            .child(
                pill("step-forward", "›", false)
                    .cursor_pointer()
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.stop_playing();
                        let f = this.structure.as_ref().map(|s| s.frame).unwrap_or(0);
                        this.set_frame(f + 1, cx);
                    })),
            )
            .child(
                div()
                    .id("scrubber")
                    .relative()
                    .flex_1()
                    .h(px(18.))
                    .flex()
                    .items_center()
                    .cursor_pointer()
                    .child(
                        canvas(move |bounds, _, _| store.set(bounds), |_, _, _, _| {})
                            .absolute()
                            .size_full(),
                    )
                    .child(
                        div()
                            .relative()
                            .w_full()
                            .h(px(4.))
                            .rounded_full()
                            .bg(theme::SUNKEN)
                            .child(
                                div()
                                    .h_full()
                                    .w(relative(t))
                                    .rounded_full()
                                    .bg(theme::INK),
                            ),
                    )
                    .on_mouse_down(
                        MouseButton::Left,
                        cx.listener(|this, e: &MouseDownEvent, _, cx| {
                            this.stop_playing();
                            this.scrubbing = true;
                            this.scrub_to(e.position.x, cx);
                            cx.stop_propagation();
                        }),
                    ),
            )
            .child(
                div()
                    .font_family(theme::MONO)
                    .text_size(px(10.5))
                    .text_color(theme::MUTED)
                    .child(format!("{} / {}", frame + 1, count)),
            )
    }
}

impl Render for Profold {
    fn render(&mut self, window: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        if viewer::spinning() {
            self.camera.yaw += 0.02;
            window.request_animation_frame();
        }
        let _timer = viewer::Phase::Tree.timer();
        let focused = self.input.read(cx).focus_handle(cx).is_focused(window);
        // Model, steps and pLDDT are the settings for Fold and Compare, so they
        // appear exactly when those can be run.
        let can_fold = matches!(
            self.detected,
            Input::Sequence(_) | Input::Pdb(_) | Input::UniProt(_)
        );
        let structure = self.structure.clone();

        div()
            .flex()
            .flex_col()
            .size_full()
            .bg(theme::PAPER)
            .text_color(theme::TEXT)
            .track_focus(&self.focus)
            .on_action(cx.listener(|this, _: &OpenFile, _, cx| this.prompt_open(cx)))
            // Files dropped anywhere on the window are opened.
            .on_drop(cx.listener(|this, paths: &ExternalPaths, _, cx| {
                if let Some(path) = paths.paths().first() {
                    this.open_path(path.clone(), cx);
                }
            }))
            .on_mouse_move(cx.listener(Self::on_mouse_move))
            .on_mouse_up(MouseButton::Left, cx.listener(Self::on_mouse_up))
            .on_mouse_up_out(MouseButton::Left, cx.listener(Self::on_mouse_up))
            .on_mouse_up(MouseButton::Right, cx.listener(Self::on_mouse_up))
            .child(self.topbar(cx))
            .child(
                div()
                .flex()
                .flex_1()
                .overflow_hidden()
            .child(
                div()
                    .w(px(PANEL_WIDTH))
                    .h_full()
                    .flex()
                    .flex_col()
                    .gap(px(16.))
                    .pt(px(20.))
                    .bg(theme::SURFACE)
                    .border_r_1()
                    .border_color(theme::BORDER)
                    .child(self.search(focused, cx))
                    .when(can_fold, |d| d.children(self.plugin_options(cx)))
                    .child(self.status(cx))
                    .when(structure.is_none(), |d| {
                        d.child(self.examples(cx)).child(div().flex_1())
                    })
                    .when_some(structure, |d, s| d.child(self.info(&s, cx)))
                    .child(div().flex_1()),
            )
            .child(self.stage(cx)),
            )
    }
}

#[cfg(test)]
mod ui_tests {
    //! The app driven as a user would, headless: real layout, simulated
    //! clicks on the real buttons, and the scene checked afterwards.

    use super::*;
    use gpui::{Modifiers, TestAppContext, VisualTestContext};

    fn open(cx: &mut TestAppContext) -> (Entity<Profold>, &mut VisualTestContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Profold::new(window, cx));
        view.update(cx, |this, cx| {
            let s = formats::load_bytes(SAMPLES[0].file, SAMPLES[0].bytes.to_vec()).unwrap();
            this.show(s, cx);
        });
        cx.run_until_parked();
        (view, cx)
    }

    fn click(cx: &mut VisualTestContext, selector: &'static str) {
        let bounds = cx
            .debug_bounds(selector)
            .unwrap_or_else(|| panic!("{selector} is not on screen"));
        cx.simulate_click(bounds.center(), Modifiers::none());
        cx.run_until_parked();
    }

    fn colors(view: &Entity<Profold>, cx: &mut VisualTestContext) -> Vec<[f32; 3]> {
        view.update(cx, |this, _| this.scene.spheres.iter().map(|s| s.color).collect())
    }

    #[gpui::test]
    fn palette_menu_recolours_elements(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        let before = colors(&view, cx);

        click(cx, "palette-button");
        assert!(
            view.update(cx, |this, _| this.open == Some(Dropdown::Palette)),
            "the palette menu did not open"
        );
        click(cx, "palette-Viridis");
        assert!(view.update(cx, |this, _| this.palette == Palette::Viridis));
        assert_ne!(colors(&view, cx), before, "the palette changed nothing");
    }

    #[gpui::test]
    fn camera_menu_resets_the_view(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        view.update(cx, |this, cx| {
            this.camera.yaw += 1.0;
            this.camera.zoom = 3.0;
            cx.notify();
        });
        cx.run_until_parked();
        click(cx, "tab-camera");
        click(cx, "reset-view");
        let (angles, zoom) = view.update(cx, |this, _| (this.camera.angles(), this.camera.zoom));
        assert_eq!(angles, [0.0, 0.0]);
        assert_eq!(zoom, 1.0);
    }

    #[gpui::test]
    fn fixed_colours_hold_elements_still(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        click(cx, "fixed-colours");
        assert!(view.update(cx, |this, _| this.fixed_colors));
        let fixed = colors(&view, cx);

        // The palette is shut while colours are fixed.
        click(cx, "palette-button");
        assert!(view.update(cx, |this, _| this.open.is_none()));
        assert_eq!(colors(&view, cx), fixed);

        click(cx, "fixed-colours");
        assert_ne!(colors(&view, cx), fixed);
    }

    #[gpui::test]
    fn picking_a_colour_paints_the_element(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        click(cx, "dot-C");
        assert!(view.update(cx, |this, _| this.picking.is_some()), "the picker did not open");
        click(cx, "swatch-1");
        let carbon = view.update(cx, |this, _| {
            this.appearance.get(&Species::Element(6)).color
        });
        assert_eq!(carbon, Some(SWATCHES[0]));
        assert!(colors(&view, cx).contains(&SWATCHES[0]) || {
            // The illustrative style saturates colours; the scene holds the
            // saturated one.
            view.update(cx, |this, _| this.style == Style::Illustrative)
        });
    }
}
