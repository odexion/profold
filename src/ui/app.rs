//! One tab: a search field and what is loaded on the left, the molecule on
//! the right.
//!
//! A tab is a viewer, or, with the SimpleFold plugin, a place to fold a
//! protein. The two share everything but the left panel: a viewer searches
//! and fetches, a fold tab takes a sequence and shows the prediction.

use std::cell::{Cell as Shared, RefCell};
use std::collections::HashMap;
use std::path::PathBuf;
use std::rc::Rc;
use std::sync::mpsc::{self, Sender};
use std::sync::Arc;
use std::time::Duration;

use gpui::{
    canvas, deferred, div, prelude::*, px, relative, Animation, AnimationExt, App,
    Bounds, ClipboardItem, Context, Entity, FocusHandle, Focusable, MouseButton,
    MouseDownEvent, MouseMoveEvent, MouseUpEvent, PathPromptOptions, Pixels, Point, ScrollDelta,
    ScrollWheelEvent, SharedString, Task, Window,
};

use crate::align;
use crate::fetch::{self, Input, Retrieval};
use crate::formats;
use crate::job::{self, Event, Outcome};
use crate::protein::{Source, Structure};
use crate::render::{
    self, Appearance, Camera, Coloring, Look, Palette, PaletteGroup, Repr, Scene, Species, Style,
};
use crate::samples::SAMPLES;
use crate::search::{self, Hit, Kind};
use crate::settings::{self, Settings};
use crate::theme;
use crate::ui::input::{InputEvent, TextInput};
use crate::ui::viewer;

#[cfg(feature = "simplefold")]
use crate::plugins::simplefold::{self, SimpleFold};

const PANEL_WIDTH: f32 = 372.0;
const TOOLBAR: f32 = 36.0;
/// How long typing must pause before a search goes out.
const SEARCH_DELAY: Duration = Duration::from_millis(320);

/// What a tab is for.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Mode {
    /// Search, fetch and open structures, and look at them.
    View,
    /// Fold a sequence with the SimpleFold plugin and look at the result.
    #[cfg_attr(not(feature = "simplefold"), allow(dead_code))]
    Fold,
}
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

pub struct Viewer {
    mode: Mode,
    focus: FocusHandle,
    input: Entity<TextInput>,
    /// The text `detected` was worked out from.
    typed: String,
    /// How what is loaded was opened, as something `reopen` takes: a path,
    /// an identifier, or `sample:` and a sample's file name. The window
    /// keeps these to reopen its tabs at the next launch.
    origin: Option<String>,
    /// The same for the job under way, until it succeeds.
    opening: Option<String>,
    detected: Input,

    /// Search results for what is in the field, while it holds a name or a
    /// sequence.
    hits: Vec<Hit>,
    /// The text the hits answer.
    hits_for: String,
    /// The text being searched for, while a search runs.
    searching: Option<String>,
    search_task: Option<Task<()>>,
    search_error: Option<SharedString>,
    /// The row Return opens.
    highlight: usize,
    /// Whether the results list is showing.
    dropdown: bool,
    /// The results list's scroll, so the highlighted row stays in view.
    results_scroll: gpui::ScrollHandle,
    /// Return was pressed before the results came; the best is opened once
    /// they do.
    open_when_found: bool,

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
    /// The logo's strand, for the head of the panel.
    mark: Arc<gpui::Image>,

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

impl Viewer {
    pub fn new(mode: Mode, window: &mut Window, cx: &mut Context<Self>) -> Self {
        let placeholder = match mode {
            Mode::View => "Search by name or gene, or paste an ID",
            Mode::Fold => "Sequence, PDB ID or accession",
        };
        let input = cx.new(|cx| TextInput::new(placeholder, cx));
        cx.subscribe_in(&input, window, |this, _, event, _window, cx| match event {
            InputEvent::Changed => this.reclassify(cx),
            InputEvent::Submit => this.submit(cx),
            InputEvent::Cancel => {
                this.dropdown = false;
                cx.notify();
            }
            InputEvent::Up => this.step_highlight(-1, cx),
            InputEvent::Down => this.step_highlight(1, cx),
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
            InputEvent::Changed | InputEvent::Up | InputEvent::Down => {}
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
            mode,
            focus,
            input,
            typed: String::new(),
            origin: None,
            opening: None,
            detected: Input::Empty,
            hits: Vec::new(),
            hits_for: String::new(),
            searching: None,
            search_task: None,
            search_error: None,
            highlight: 0,
            dropdown: false,
            results_scroll: gpui::ScrollHandle::new(),
            open_when_found: false,
            structure: None,
            compare: None,
            fit: None,
            show_reference: true,
            show_prediction: true,
            scene: Arc::new(Scene::default()),
            camera: Camera {
                orientation: Camera::home(),
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
            mark: Arc::new(gpui::Image::from_bytes(
                gpui::ImageFormat::Png,
                include_bytes!("../../assets/mark.png").to_vec(),
            )),
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
        let text = self.input.read(cx).text().to_string();
        // Setting the text in code reports a change of its own afterwards.
        if text == self.typed {
            return;
        }
        self.typed = text.clone();
        self.detected = fetch::classify(&text);
        // A fold tab has nothing to search, so words that could be residues
        // are residues there.
        if self.mode == Mode::Fold {
            if let Input::Query(_) | Input::Unknown = self.detected {
                if let Some(seq) = fetch::as_sequence(&text, true) {
                    self.detected = Input::Sequence(seq);
                }
            }
        }
        self.error = None;
        self.open_when_found = false;
        match self.search_text() {
            Some(query) if query == self.hits_for => self.dropdown = true,
            Some(query) => self.schedule_search(query, SEARCH_DELAY, cx),
            None => {
                self.search_task = None;
                self.searching = None;
                self.dropdown = false;
            }
        }
        cx.notify();
    }

    /// What there is to search for in the field: words, or a sequence to
    /// find in the PDB. Only a viewer searches.
    fn search_text(&self) -> Option<String> {
        if self.mode != Mode::View {
            return None;
        }
        match &self.detected {
            Input::Query(q) if q.chars().filter(|c| c.is_alphanumeric()).count() >= 3 => {
                Some(q.clone())
            }
            Input::Sequence(seq) => Some(seq.clone()),
            _ => None,
        }
    }

    /// Search once typing has paused for `delay`. A newer search replaces
    /// this one, and results for text no longer in the field are dropped.
    fn schedule_search(&mut self, query: String, delay: Duration, cx: &mut Context<Self>) {
        self.searching = Some(query.clone());
        self.search_error = None;
        let sequence = matches!(self.detected, Input::Sequence(_));
        self.search_task = Some(cx.spawn(async move |this, cx| {
            cx.background_executor().timer(delay).await;
            let asked = query.clone();
            let found = cx
                .background_executor()
                .spawn(async move {
                    if sequence {
                        search::search_sequence(&asked)
                    } else {
                        search::search(&asked)
                    }
                })
                .await;
            let _ = this.update(cx, |this, cx| this.receive(query, found, cx));
        }));
    }

    /// Results are back for `query`. Enter may have gone ahead of them: then
    /// a hit named exactly what was typed opens outright, and anything else
    /// is left in the list to choose from.
    fn receive(&mut self, query: String, found: anyhow::Result<Vec<Hit>>, cx: &mut Context<Self>) {
        if self.searching.as_deref() != Some(query.as_str()) {
            return;
        }
        self.searching = None;
        match found {
            Ok(hits) => {
                self.hits = hits;
                self.highlight = 0;
                self.dropdown = true;
                if std::mem::take(&mut self.open_when_found) {
                    if let Some(i) = self.hits.iter().position(|h| h.is_named(&query)) {
                        self.pick(i, cx);
                    }
                }
                self.hits_for = query;
            }
            Err(e) => {
                self.hits.clear();
                self.hits_for = query;
                self.dropdown = true;
                self.open_when_found = false;
                self.search_error = Some(format!("{e:#}").into());
            }
        }
        cx.notify();
    }

    /// Move the highlighted result up or down, opening the list if it was
    /// put away.
    fn step_highlight(&mut self, by: isize, cx: &mut Context<Self>) {
        if self.hits.is_empty() || self.search_text().as_deref() != Some(self.hits_for.as_str()) {
            return;
        }
        if !self.dropdown {
            self.dropdown = true;
        } else {
            let last = self.hits.len() as isize - 1;
            self.highlight = (self.highlight as isize + by).clamp(0, last) as usize;
        }
        // The list's children are the rows with a heading before each
        // group, so the row's place counts the headings above it.
        let i = self.highlight;
        let headings = (0..=i)
            .filter(|&k| k == 0 || self.hits[k].kind != self.hits[k - 1].kind)
            .count();
        self.results_scroll.scroll_to_item(i + headings);
        cx.notify();
    }

    /// Open a search result in this tab. What was typed stays in the field,
    /// so another result is a click away.
    fn pick(&mut self, i: usize, cx: &mut Context<Self>) {
        let Some(hit) = self.hits.get(i).cloned() else {
            return;
        };
        self.highlight = i;
        self.dropdown = false;
        self.note(format!("{} · {}", hit.label(), hit.title));
        self.retrieve(hit.retrieval(), cx);
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
        self.origin = None;
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

    /// Draw again after the theme has changed: the scene's colours are
    /// settled for the page they sit on.
    pub fn retheme(&mut self, cx: &mut Context<Self>) {
        self.rebuild(cx);
    }

    /// Show the structure another way. A backbone style is refused for a
    /// structure without a backbone, as its pill is greyed out.
    pub fn set_repr(&mut self, repr: Repr, cx: &mut Context<Self>) {
        let Some(s) = self.structure.as_ref() else { return };
        if repr.needs_polymer() && !s.has_polymer() {
            return;
        }
        self.repr = repr;
        self.rebuild(cx);
    }

    pub fn set_style(&mut self, style: Style, cx: &mut Context<Self>) {
        if self.structure.is_none() {
            return;
        }
        self.style = style;
        self.rebuild(cx);
    }

    /// Bring the structure back to the middle, keeping the angle and zoom.
    pub fn recentre(&mut self, cx: &mut Context<Self>) {
        self.camera.pan = [0.0, 0.0];
        self.open = None;
        cx.notify();
    }

    /// Back to the angle and framing it opened at.
    pub fn reset_view(&mut self, cx: &mut Context<Self>) {
        self.halt();
        // Framed before the camera is touched, so the borrow of the
        // structures is over by then.
        let framed = {
            let all: Vec<&Structure> = self.visible().into_iter().map(|s| s.as_ref()).collect();
            (!all.is_empty()).then(|| Camera::frame_all(&all))
        };
        if let Some(camera) = framed {
            self.camera = camera;
        }
        self.open = None;
        cx.notify();
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
                self.origin = self.opening.take();
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
                self.opening = None;
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
        if self.mode == Mode::Fold {
            self.submit_fold(cx);
            return;
        }
        match self.detected.clone() {
            Input::Empty | Input::Unknown => {}
            Input::Missing(path) => {
                self.error = Some(format!("No file at {path}").into());
                cx.notify();
            }
            Input::File(path) => self.open_path(path, cx),
            Input::Pdb(id) => self.retrieve(Retrieval::Pdb(id), cx),
            Input::UniProt(acc) => self.retrieve(Retrieval::AlphaFold(acc), cx),
            Input::Cid(cid) => self.retrieve(Retrieval::PubChem(cid), cx),
            Input::Query(_) | Input::Sequence(_) => {
                let Some(query) = self.search_text() else { return };
                if query == self.hits_for && self.searching.is_none() {
                    if self.hits.is_empty() {
                        self.dropdown = true;
                        cx.notify();
                    } else {
                        self.pick(self.highlight.min(self.hits.len() - 1), cx);
                    }
                } else {
                    // Not back yet, or still waiting out the typing pause:
                    // go now, and open it when it comes if it is named
                    // exactly what was typed.
                    self.open_when_found = true;
                    self.schedule_search(query, Duration::ZERO, cx);
                    cx.notify();
                }
            }
        }
    }

    #[cfg(feature = "simplefold")]
    fn start(&mut self, job: simplefold::Job, cx: &mut Context<Self>) {
        for note in job.notes {
            self.note(note);
        }
        self.run_job(job.stage, job.work, cx);
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
            Input::Empty => return,
            _ => {
                self.error = Some(
                    "Folding takes a protein sequence, or a PDB ID or UniProt accession to fold the chains of."
                        .into(),
                );
                cx.notify();
                return;
            }
        };
        self.start(job, cx);
    }

    #[cfg(not(feature = "simplefold"))]
    fn submit_fold(&mut self, _: &mut Context<Self>) {}

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
                        Ok(Event::Heartbeat) => {}
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
        if self.busy.is_some() {
            return;
        }
        self.error = None;
        self.opening = Some(match &what {
            Retrieval::Pdb(id) | Retrieval::AlphaFold(id) => id.clone(),
            Retrieval::PubChem(cid) => format!("CID {cid}"),
        });
        // Without a model, an accession can still be predicted; say where.
        let hint = cfg!(feature = "simplefold") && matches!(what, Retrieval::AlphaFold(_));
        self.run_job(
            "Fetching",
            move |tx| {
                what.run(tx).map_err(|e| {
                    if hint {
                        anyhow::anyhow!("{e:#}. A Fold tab (+ in the tab bar) can predict it.")
                    } else {
                        e
                    }
                })
            },
            cx,
        );
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
        let absolute = std::fs::canonicalize(&path).unwrap_or_else(|_| path.clone());
        self.opening = Some(absolute.to_string_lossy().into_owned());
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
        self.opening = Some(format!("sample:{}", sample.file));
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
    pub fn prompt_open(&mut self, cx: &mut Context<Self>) {
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
            dark: theme::dark(),
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
                    let mut angles = self.camera.angles();
                    angles[axis] = deg;
                    self.camera.set_angles(angles);
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
                .bg(theme::surface())
                .border_1()
                .border_color(theme::border_focus())
                .font_family(theme::MONO)
                .text_size(px(10.5))
                .text_color(theme::ink())
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
            .text_color(theme::muted())
            .cursor(gpui::CursorStyle::IBeam)
            .hover(|s| s.bg(theme::sunken()).text_color(theme::ink()))
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
    pub fn save_image(&mut self, scale: usize, cx: &mut Context<Self>) {
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
                        render::render(&scene, &camera, width, height, theme::view(), style);
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

    /// How what is loaded was opened, for reopening it at the next launch.
    pub fn origin(&self) -> Option<&str> {
        // What is still being opened counts: a session saved mid-download
        // should not forget the tab.
        match &self.opening {
            Some(opening) if self.busy.is_some() => Some(opening),
            _ => self.origin.as_deref().filter(|_| self.structure.is_some()),
        }
    }

    /// Open again what `origin` described.
    pub fn reopen(&mut self, origin: &str, cx: &mut Context<Self>) {
        if let Some(file) = origin.strip_prefix("sample:") {
            if let Some(i) = SAMPLES.iter().position(|s| s.file == file) {
                self.open_sample(i, cx);
            }
            return;
        }
        let path = PathBuf::from(origin);
        if path.is_absolute() {
            if path.is_file() {
                self.open_path(path, cx);
            }
            return;
        }
        self.open(origin, cx);
    }

    /// Put text in the field as though it had been typed.
    pub fn type_text(&mut self, text: &str, cx: &mut Context<Self>) {
        let text = text.to_string();
        self.input.update(cx, |input, cx| input.set_text(text, cx));
        self.reclassify(cx);
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
                this.camera.rotate(v[0], v[1]);
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
                self.camera.rotate(self.last_step[0], self.last_step[1]);
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

impl Focusable for Viewer {
    fn focus_handle(&self, _: &App) -> FocusHandle {
        self.focus.clone()
    }
}

// ---------------------------------------------------------------------------
// View
// ---------------------------------------------------------------------------

/// Things to type, listed on an empty viewer.
const EXAMPLES: [(&str, &str); 5] = [
    ("green fluorescent protein", "search by name"),
    ("caffeine", "a molecule from PubChem"),
    ("4HHB", "haemoglobin, by PDB ID"),
    ("P69905", "an AlphaFold model, by accession"),
    ("CID 2244", "aspirin, by PubChem CID"),
];

/// Things to fold, listed on an empty fold tab.
#[cfg(feature = "simplefold")]
const FOLD_EXAMPLES: [(&str, &str); 3] = [
    simplefold::EXAMPLE,
    ("1CRN", "crambin's chain, 46 aa"),
    ("P69905", "human haemoglobin α"),
];
#[cfg(not(feature = "simplefold"))]
const FOLD_EXAMPLES: [(&str, &str); 0] = [];

pub(crate) fn pill(id: &'static str, label: impl Into<SharedString>, active: bool) -> gpui::Stateful<gpui::Div> {
    div()
        .id(id)
        .px(px(9.))
        .py(px(3.))
        .rounded(px(5.))
        .text_size(px(11.))
        .when(active, |d| d.bg(theme::ink()).text_color(theme::on_ink()))
        .when(!active, |d| {
            d.text_color(theme::muted())
                .hover(|s| s.bg(theme::sunken()).text_color(theme::text()))
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
        .when(active, |d| d.bg(theme::ink()).text_color(theme::on_ink()))
        .when(!active, |d| {
            d.text_color(theme::muted())
                .hover(|s| s.bg(theme::sunken()).text_color(theme::text()))
        })
        .child(label.into())
}

/// A one-line label on a dark ground, for hover hints.
pub(crate) struct Tip(pub(crate) SharedString);

impl Render for Tip {
    fn render(&mut self, _: &mut Window, _: &mut Context<Self>) -> impl IntoElement {
        div()
            .px(px(7.))
            .py(px(3.))
            .rounded(px(5.))
            .bg(theme::ink())
            .text_size(px(10.5))
            .text_color(theme::on_ink())
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
        .when(open, |d| d.bg(theme::sunken()).text_color(theme::ink()))
        .when(!open, |d| {
            d.text_color(theme::muted())
                .hover(|s| s.bg(theme::sunken()).text_color(theme::text()))
        })
        .child(text)
}

/// The card a menu is drawn on, hanging under its title.
pub(crate) fn menu(width: f32) -> gpui::Div {
    div()
        .absolute()
        .top(px(24.))
        .left(px(-3.))
        .w(px(width))
        .p(px(3.))
        .rounded(px(8.))
        .bg(theme::surface())
        .border_1()
        .border_color(theme::border())
        .shadow_md()
        .flex()
        .flex_col()
        .gap(px(1.))
        .on_mouse_down(MouseButton::Left, swallow)
}

/// Keep a press on a floating card to itself. A menu hangs out over the
/// viewport, which closes menus on mouse down; without this the menu would
/// be gone before the click on its item finished.
pub(crate) fn swallow(_: &MouseDownEvent, _: &mut Window, cx: &mut App) {
    cx.stop_propagation();
}

/// One line of a menu: a name, a tick when it is the one in use.
pub(crate) fn item(id: impl Into<gpui::ElementId>, text: String, chosen: bool) -> gpui::Stateful<gpui::Div> {
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
        .when(!chosen, |d| d.hover(|s| s.bg(theme::sunken())))
        .child(
            div()
                .text_size(px(11.5))
                .text_color(if chosen { theme::ink() } else { theme::text() })
                .child(text),
        )
        .when(chosen, |d| {
            d.child(div().text_size(px(10.)).text_color(theme::muted()).child("✓"))
        })
}

/// The hairline between groups of menu items.
pub(crate) fn divide() -> impl IntoElement {
    div().my(px(3.)).h(px(1.)).bg(theme::border())
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
        .border_color(theme::border())
        .text_color(theme::muted())
        .cursor_pointer()
        .hover(|s| s.bg(theme::sunken()).text_color(theme::text()))
        .child(text)
}

pub(crate) fn label(text: impl Into<gpui::SharedString>) -> impl IntoElement {
    div()
        .text_size(px(9.5))
        .text_color(theme::faint())
        .child(text.into())
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
        .when(on, |d| d.bg(theme::sunken()))
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
                .text_color(if on { theme::text() } else { theme::faint() })
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
    div().h(px(1.)).w_full().bg(theme::border())
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
        .text_color(theme::faint())
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
        .border_color(theme::border())
        .child(key(field))
        .child(
            div()
                .when(mono, |d| d.font_family(theme::MONO))
                .text_size(px(11.5))
                .text_color(theme::ink())
                .child(value),
        )
}

/// Text cut to a number of characters with an ellipsis. gpui measures a
/// line that does not wrap once, at its full width, and keeps that size, so
/// its own ellipsis never appears in a flexible row; this does it by count.
fn clip(text: &str, max: usize) -> String {
    if text.chars().count() <= max {
        return text.to_string();
    }
    let cut: String = text.chars().take(max - 1).collect();
    format!("{}…", cut.trim_end())
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

impl Viewer {
    fn toggle(&mut self, which: Dropdown, cx: &mut Context<Self>) {
        self.open = (self.open != Some(which)).then_some(which);
        cx.notify();
    }

    /// The bar across the top of the tab: the menus that govern the
    /// viewport, and what is loaded.
    fn toolbar(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let viewing = self.mode == Mode::View;
        let loaded = self.structure.is_some();
        let badge = self.structure.as_ref().map(|s| match s.format {
            Some(format) if s.source.is_local() => {
                format!("{} · {format}", s.id)
            }
            _ => format!("{} · {}", s.id.to_uppercase(), s.source.label()),
        });

        div()
            .id("toolbar")
            .flex()
            .flex_none()
            .items_center()
            .gap(px(2.))
            .h(px(TOOLBAR))
            .pl(px(12.))
            .pr(px(16.))
            .bg(theme::surface())
            .border_b_1()
            .border_color(theme::border())
            .when(viewing, |d| {
                d.child(
                    tab("open-file", "Open", false)
                        .on_click(cx.listener(|this, _, _, cx| this.prompt_open(cx))),
                )
                .child(self.tab_with_menu(
                    Dropdown::Samples,
                    "samples",
                    "Samples",
                    self.samples_menu(cx).into_any_element(),
                    cx,
                ))
            })
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
                        .text_color(theme::faint())
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
        let (fps, live) = self.cache.borrow().fps();
        let name = |n: &'static str| div().text_size(px(9.)).text_color(theme::faint()).child(n);
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
            .child(div().w(px(6.)))
            .child(name("FPS"))
            .child(
                div()
                    .w(px(30.))
                    .h(px(18.))
                    .px(px(4.))
                    .flex()
                    .items_center()
                    .justify_end()
                    .font_family(theme::MONO)
                    .text_size(px(10.5))
                    .text_color(if live { theme::muted() } else { theme::faint() })
                    .child(format!("{}", fps.round() as u32)),
            )
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
            .children(crate::samples::grouped().into_iter().enumerate().flat_map(
                |(n, (group, samples))| {
                    let heading = div()
                        .px(px(8.))
                        .pt(px(if n == 0 { 3. } else { 0. }))
                        .pb(px(2.))
                        .child(label(group.to_uppercase()))
                        .into_any_element();
                    let rows = samples.into_iter().map(|(i, sample)| {
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
                            .hover(|s| s.bg(theme::sunken()))
                            .child(
                                div()
                                    .text_size(px(11.5))
                                    .text_color(theme::text())
                                    .child(sample.name),
                            )
                            .child(
                                div()
                                    .text_size(px(9.5))
                                    .text_color(theme::faint())
                                    .child(sample.file),
                            )
                            .on_click(cx.listener(move |this, _, _, cx| {
                                this.open = None;
                                this.open_sample(i, cx);
                            }))
                            .into_any_element()
                    });
                    (n > 0)
                        .then(|| divide().into_any_element())
                        .into_iter()
                        .chain(std::iter::once(heading))
                        .chain(rows)
                        .collect::<Vec<_>>()
                },
            ))
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
                .bg(theme::surface())
                .border_1()
                .border_color(theme::border())
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
                            .when(!enabled, |d| d.text_color(theme::faint()))
                            .when(enabled, |d| d.cursor_pointer())
                            .on_click(cx.listener(move |this, _, _, cx| this.set_repr(r, cx)))
                            .into_any_element()
                    })
                    .collect(),
            ))
            .child(group(
                options
                    .iter()
                    // pLDDT means nothing without a predicted model, so it
                    // is offered only for one.
                    .filter(|c| **c != Coloring::Confidence || has_confidence)
                    .map(|c| {
                        let c = *c;
                        let name = match c {
                            Coloring::Chain if !polymer => "Molecule",
                            _ => c.label(),
                        };
                        pill_keyed(("colour", c as usize), name, c == coloring)
                            .cursor_pointer()
                            .on_click(cx.listener(move |this, _, _, cx| {
                                this.coloring = c;
                                this.rebuild(cx);
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
                            .on_click(cx.listener(move |this, _, _, cx| this.set_style(st, cx)))
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
        let item = |p: Palette| {
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
                .when(!chosen, |d| d.hover(|s| s.bg(theme::sunken())))
                .child(swatch_strip(p))
                .child(
                    div()
                        .flex_1()
                        .text_size(px(11.5))
                        .text_color(if chosen { theme::ink() } else { theme::text() })
                        .child(p.label()),
                )
                .when(chosen, |d| {
                    d.child(div().text_size(px(10.)).text_color(theme::muted()).child("✓"))
                })
                .on_click(cx.listener(move |this, _, _, cx| {
                    this.palette = p;
                    this.open = None;
                    this.rebuild(cx);
                }))
                .into_any_element()
        };
        // Each family under its name. The two rainbows lead in a row across
        // the top; below, Perceptual's seven stand beside Thematic and Muted,
        // which together come out the same height, so the menu is short
        // enough to open upward in a small window and squares off at the foot.
        // A heading after another group stands a little further off it.
        let heading = |g: PaletteGroup, first: bool| {
            div()
                .px(px(8.))
                .pt(px(if first { 6. } else { 10. }))
                .pb(px(2.))
                .text_size(px(9.))
                .text_color(theme::faint())
                .child(spaced(&g.label().to_uppercase()))
        };
        let members = |g: PaletteGroup| Palette::ALL.into_iter().filter(move |p| p.group() == g);
        let column = |groups: &[PaletteGroup]| {
            let mut col = div().flex_1().flex().flex_col().gap(px(1.));
            for (i, g) in groups.iter().enumerate() {
                col = col.child(heading(*g, i == 0));
                for p in members(*g) {
                    col = col.child(item(p));
                }
            }
            col
        };
        let top = div()
            .flex()
            .flex_col()
            .gap(px(1.))
            .child(heading(PaletteGroup::Classic, true))
            .child(
                div()
                    .flex()
                    .gap(px(4.))
                    .children(members(PaletteGroup::Classic).map(|p| div().flex_1().child(item(p)))),
            );
        let left = column(&[PaletteGroup::Perceptual]);
        let right = column(&[PaletteGroup::Thematic, PaletteGroup::Muted]);
        let list = menu(372.)
            .top_auto()
            .bottom(px(34.))
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.open = None;
                cx.notify();
            }))
            .child(top)
            .child(div().flex().items_start().gap(px(4.)).child(left).child(right));
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
                    .bg(theme::surface())
                    .border_1()
                    .border_color(theme::border())
                    .shadow_sm()
                    .when(fixed, |d| d.opacity(0.4))
                    .when(!fixed, |d| d.cursor_pointer().hover(|s| s.bg(theme::sunken())))
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
                            .text_color(theme::text())
                            .child(current.label()),
                    )
                    .child(
                        div()
                            .text_size(px(12.))
                            .text_color(theme::muted())
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
                item("recentre", "Centre".into(), false)
                    .on_click(cx.listener(|this, _, _, cx| this.recentre(cx))),
            )
            .child(
                item("reset-view", "Reset view".into(), false)
                    .debug_selector(|| "reset-view".into())
                    .on_click(cx.listener(|this, _, _, cx| this.reset_view(cx))),
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
                    .hover(|s| s.bg(theme::sunken()))
                    .child(
                        div()
                            .text_size(px(11.5))
                            .text_color(theme::text())
                            .child(format!("Save PNG · {scale}×")),
                    )
                    .child(
                        div()
                            .font_family(theme::MONO)
                            .text_size(px(9.5))
                            .text_color(theme::faint())
                            .child(size),
                    )
                    .on_click(cx.listener(move |this, _, _, cx| {
                        this.open = None;
                        this.save_image(scale, cx);
                    }))
            }))
    }

    /// The logo and what the tab is for, at the head of the panel.
    fn masthead(&self) -> impl IntoElement {
        let (name, purpose) = match self.mode {
            Mode::View => ("profold", "proteins, molecules and crystals"),
            Mode::Fold => ("Fold a protein", "structure prediction with Apple's SimpleFold"),
        };
        div()
            .flex()
            .items_center()
            .gap(px(11.))
            .px(px(20.))
            .child(gpui::img(self.mark.clone()).size(px(40.)).flex_none())
            .child(
                div()
                    .flex()
                    .flex_col()
                    .gap(px(1.))
                    .child(div().text_size(px(17.)).text_color(theme::ink()).child(name))
                    .child(
                        div()
                            .text_size(px(10.5))
                            .text_color(theme::faint())
                            .child(purpose),
                    ),
            )
    }

    fn search(&self, focused: bool, cx: &mut Context<Self>) -> impl IntoElement {
        let detected = self.detected.clone();
        let folding = self.mode == Mode::Fold;
        let ready = if folding {
            matches!(detected, Input::Sequence(_) | Input::Pdb(_) | Input::UniProt(_))
        } else {
            match &detected {
                Input::Empty | Input::Unknown | Input::Missing(_) => false,
                Input::Query(_) => self.search_text().is_some(),
                _ => true,
            }
        };
        // An identifier in a fold tab can be folded, or folded and laid over
        // the deposited structure.
        let comparable = folding && matches!(detected, Input::Pdb(_) | Input::UniProt(_));
        let busy = self.busy.is_some();
        let action = if folding { "Fold" } else { detected.action() };

        let run = div()
            .id("run")
            .debug_selector(|| "run".into())
            .flex_none()
            .h(px(38.))
            .min_w(px(68.))
            .px(px(15.))
            .flex()
            .items_center()
            .justify_center()
            .rounded(px(8.))
            .text_size(px(11.5))
            .when(ready && !busy, |d| {
                d.bg(theme::ink())
                    .text_color(theme::on_ink())
                    .hover(|s| s.opacity(0.86))
                    .cursor_pointer()
            })
            .when(!ready || busy, |d| d.bg(theme::sunken()).text_color(theme::faint()))
            .child(action)
            .on_click(cx.listener(|this, _, _, cx| this.submit(cx)));

        let status = match &detected {
            Input::Empty => None,
            Input::Unknown => Some(
                div()
                    .text_size(px(10.5))
                    .text_color(theme::faint())
                    .child(detected.detail())
                    .into_any_element(),
            ),
            Input::Query(_) if folding => Some(
                div()
                    .text_size(px(10.5))
                    .text_color(theme::faint())
                    .child("not a sequence, PDB ID or accession")
                    .into_any_element(),
            ),
            other => Some(
                div()
                    .flex()
                    .items_center()
                    .gap(px(6.))
                    .min_w_0()
                    .child(
                        div()
                            .flex_none()
                            .text_size(px(9.))
                            .text_color(theme::muted())
                            .child(match other {
                                Input::Sequence(_) if !folding => "SEQUENCE SEARCH",
                                _ => other.kind(),
                            }),
                    )
                    .child(
                        div()
                            .font_family(theme::MONO)
                            .text_size(px(11.))
                            .text_color(theme::text())
                            .truncate()
                            .child(match other {
                                Input::Sequence(s) if !folding => format!("{} aa in the PDB", s.len()),
                                _ => other.detail(),
                            }),
                    )
                    .when(self.searching.is_some(), |d| {
                        d.child(
                            div()
                                .flex_none()
                                .text_size(px(10.))
                                .text_color(theme::faint())
                                .child("searching…"),
                        )
                    })
                    .into_any_element(),
            ),
        };

        div()
            .relative()
            .flex()
            .flex_col()
            .gap(px(10.))
            .px(px(20.))
            .child(
                div()
                    .flex()
                    .items_center()
                    .gap(px(8.))
                    .child(
                        div()
                            .id("search-field")
                            .flex_1()
                            .min_w_0()
                            .flex()
                            .items_center()
                            .h(px(38.))
                            .px(px(11.))
                            .rounded(px(8.))
                            .bg(if focused { theme::surface() } else { theme::sunken() })
                            .border_1()
                            .border_color(if focused {
                                theme::border_focus()
                            } else {
                                theme::border()
                            })
                            .text_size(px(12.5))
                            .text_color(theme::text())
                            .font_family(theme::MONO)
                            // Back into the field, the list comes back with it.
                            .on_mouse_down(
                                MouseButton::Left,
                                cx.listener(|this, _, _, cx| {
                                    if !this.hits.is_empty()
                                        && this.search_text().as_deref() == Some(this.hits_for.as_str())
                                    {
                                        this.dropdown = true;
                                        cx.notify();
                                    }
                                }),
                            )
                            .child(self.input.clone()),
                    )
                    .child(run),
            )
            // What the field holds, and what else can be done with it; only
            // once there is something in it.
            .when_some(status, |d, status| {
                d.child(
                    div()
                        .flex()
                        .items_center()
                        .justify_between()
                        .gap(px(8.))
                        .h(px(22.))
                        .child(status)
                        .when(comparable && !busy, |d| self.plugin_actions(d, cx)),
                )
            })
            .when(self.dropdown && !folding && !busy, |d| {
                d.child(deferred(self.results(cx)).with_priority(2))
            })
    }

    /// The search results, hanging under the field and out over the
    /// viewport, grouped by where they come from.
    fn results(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let mut rows: Vec<gpui::AnyElement> = Vec::new();
        let mut last: Option<Kind> = None;
        for (i, hit) in self.hits.iter().enumerate() {
            if last != Some(hit.kind) {
                rows.push(
                    div()
                        .px(px(10.))
                        .pt(px(if last.is_none() { 6. } else { 10. }))
                        .pb(px(3.))
                        .child(label(spaced(&hit.kind.heading().to_uppercase())))
                        .into_any_element(),
                );
                last = Some(hit.kind);
            }
            let chosen = i == self.highlight;
            rows.push(
                div()
                    .id(("hit", i))
                    .debug_selector(move || format!("hit-{i}"))
                    .flex()
                    .flex_col()
                    .gap(px(1.))
                    .flex_none()
                    .px(px(10.))
                    .py(px(4.))
                    .rounded(px(6.))
                    .cursor_pointer()
                    .when(chosen, |d| d.bg(theme::sunken()))
                    .hover(|s| s.bg(theme::sunken()))
                    .child(
                        div()
                            .flex()
                            .items_baseline()
                            .gap(px(8.))
                            .child(
                                div()
                                    .flex_none()
                                    .w(px(96.))
                                    .whitespace_nowrap()
                                    .font_family(theme::MONO)
                                    .text_size(px(10.5))
                                    .text_color(theme::muted())
                                    .child(hit.label()),
                            )
                            .child(
                                div().flex_1().min_w_0().overflow_hidden().child(
                                    div()
                                        .w_full()
                                        .truncate()
                                        .text_size(px(12.))
                                        .text_color(theme::ink())
                                        .child(clip(&hit.title, if i == 0 { 50 } else { 62 })),
                                ),
                            )
                            .when(i == 0, |d| {
                                d.child(
                                    div()
                                        .flex_none()
                                        .text_size(px(9.5))
                                        .text_color(theme::faint())
                                        .child("best match ↵"),
                                )
                            }),
                    )
                    .child(
                        div().pl(px(104.)).w_full().overflow_hidden().child(
                            div()
                                .w_full()
                                .truncate()
                                .text_size(px(10.5))
                                .text_color(theme::faint())
                                .child(clip(&hit.detail, 72)),
                        ),
                    )
                    .on_mouse_move(cx.listener(move |this, _, _, cx| {
                        if this.highlight != i {
                            this.highlight = i;
                            cx.notify();
                        }
                    }))
                    .on_click(cx.listener(move |this, _, _, cx| this.pick(i, cx)))
                    .into_any_element(),
            );
        }

        let empty = self.hits.is_empty();
        let note = match (&self.search_error, empty, self.searching.is_some()) {
            (Some(err), _, _) => Some(err.clone()),
            (None, true, true) => Some("Searching…".into()),
            (None, true, false) => Some(format!("No matches for “{}”", self.hits_for).into()),
            _ => None,
        };

        div()
            .absolute()
            .top(px(44.))
            .left(px(20.))
            .id("results")
            .track_scroll(&self.results_scroll)
            .w(px(520.))
            .max_h(px(560.))
            .p(px(4.))
            .rounded(px(9.))
            .bg(theme::surface())
            .border_1()
            .border_color(theme::border())
            .shadow_lg()
            .flex()
            .flex_col()
            .overflow_y_scroll()
            .on_mouse_down(MouseButton::Left, swallow)
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.dropdown = false;
                cx.notify();
            }))
            .children(rows)
            .when_some(note, |d, note| {
                d.child(
                    div()
                        .px(px(10.))
                        .py(px(8.))
                        .text_size(px(11.))
                        .text_color(if self.search_error.is_some() {
                            theme::danger()
                        } else {
                            theme::faint()
                        })
                        .child(note),
                )
            })
            .when(!empty, |d| {
                d.child(
                    div()
                        .px(px(10.))
                        .pt(px(6.))
                        .pb(px(4.))
                        .text_size(px(9.5))
                        .text_color(theme::faint())
                        .child("↑↓ to choose · ↵ to open · esc to close"),
                )
            })
    }

    /// Compare beside Fold, for an identifier: fold its own sequence and lay
    /// the prediction over the deposited structure.
    #[cfg(feature = "simplefold")]
    fn plugin_actions(&self, d: gpui::Div, cx: &mut Context<Self>) -> gpui::Div {
        d.child(
            secondary("compare", "Compare")
                .debug_selector(|| "compare".into())
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
                d.child(div().text_size(px(11.)).text_color(theme::danger()).child(err))
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
            .bg(theme::sunken())
            .overflow_hidden()
            .child(match fraction {
                Some(f) => div()
                    .h_full()
                    .w(relative(f.clamp(0.005, 1.0)))
                    .rounded_full()
                    .bg(theme::ink())
                    .into_any_element(),
                None => div()
                    .id("indeterminate")
                    .absolute()
                    .top_0()
                    .bottom_0()
                    .w(relative(0.28))
                    .rounded_full()
                    .bg(theme::ink())
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
            .bg(theme::paper())
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
                                    .text_color(theme::ink())
                                    .child(busy),
                            )
                            .when_some(percent, |d, p| {
                                d.child(
                                    div()
                                        .font_family(theme::MONO)
                                        .text_size(px(11.5))
                                        .text_color(theme::muted())
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
                            .text_color(theme::muted())
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
                            .text_color(theme::faint())
                            .children(self.recent(6).into_iter().map(|line| {
                                div().overflow_hidden().whitespace_nowrap().child(line)
                            })),
                    ),
            )
    }

    fn examples(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let folding = self.mode == Mode::Fold;
        let (heading, examples) = if folding {
            ("TRY", &FOLD_EXAMPLES[..])
        } else {
            ("SEARCH OR FETCH", &EXAMPLES[..])
        };
        div()
            .flex()
            .flex_col()
            .gap(px(7.))
            .px(px(20.))
            .when(!folding, |d| d.child(self.sample_list(cx)).child(div().h(px(8.))))
            .child(label(heading))
            .children(examples.iter().enumerate().map(|(i, (value, note))| {
                let value = value.to_string();
                let shown = if value.len() > 26 {
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
                            .text_color(theme::text())
                            .hover(|s| s.text_color(theme::ink()))
                            .child(shown),
                    )
                    .child(
                        div()
                            .text_size(px(10.5))
                            .text_color(theme::faint())
                            .child(*note),
                    )
                    .on_click(cx.listener(move |this, _, _, cx| this.open(&value, cx)))
            }))
    }

    fn sample_list(&self, cx: &mut Context<Self>) -> impl IntoElement {
        div()
            .flex()
            .flex_col()
            .gap(px(7.))
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
                            .text_color(theme::text())
                            .hover(|s| s.text_color(theme::ink()))
                            .child(sample.name),
                    )
                    .child(
                        div()
                            .text_size(px(10.5))
                            .text_color(theme::faint())
                            .child(sample.note),
                    )
                    .on_click(cx.listener(move |this, _, _, cx| this.open_sample(i, cx)))
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
                        .bg(theme::surface())
                        .border_1()
                        .border_color(theme::border())
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
                                .text_color(theme::ink())
                                .child(format!("{rmsd:.2} Å")),
                        )
                        .child(
                            div()
                                .text_size(px(10.))
                                .text_color(theme::muted())
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
                                        d.border_1().border_color(theme::faint())
                                    }),
                            )
                            .child(
                                div()
                                    .flex_1()
                                    .text_size(px(11.5))
                                    .text_color(theme::text())
                                    .child(name),
                            )
                            .child(
                                div()
                                    .font_family(theme::MONO)
                                    .text_size(px(10.))
                                    .text_color(theme::faint())
                                    .child(format!("{} aa", chain.residues.len())),
                            )
                    }),
            )
            .when(rest > 0, |d| {
                d.child(
                    div()
                        .pl(px(16.))
                        .text_size(px(10.5))
                        .text_color(theme::faint())
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
                        d.text_color(theme::muted())
                            .cursor_pointer()
                            .hover(|s| s.text_color(theme::ink()))
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
                                    .when(picking, |d| d.border_1().border_color(theme::ink()))
                                    .hover(|d| d.bg(theme::sunken()))
                                    .child(
                                        div()
                                            .size(px(8.))
                                            .rounded_full()
                                            .when(by_element, |d| {
                                                d.bg(swatch(shown)).opacity(look.alpha.max(0.2))
                                            })
                                            .when(!by_element, |d| {
                                                d.border_1().border_color(theme::faint())
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
                                    .text_color(theme::text())
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
                            .text_color(theme::faint())
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
                    .when(on, |d| d.bg(theme::ink()).justify_end())
                    .when(!on, |d| d.bg(theme::border()))
                    .child(div().size(px(8.)).rounded_full().bg(theme::surface())),
            )
            .child(
                div()
                    .text_size(px(11.))
                    .text_color(if on { theme::ink() } else { theme::muted() })
                    .child("Fixed colours"),
            )
            .child(
                div()
                    .text_size(px(9.5))
                    .text_color(theme::faint())
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
                .when(chosen, |d| d.border_1().border_color(theme::ink()))
                .when(!chosen, |d| d.hover(|s| s.bg(theme::sunken())))
                .child(div().size_full().rounded(px(3.)).bg(swatch(c)))
                .on_click(cx.listener(move |this, _, _, cx| {
                    this.set_color(species.clone(), set, cx)
                }))
        };
        let label = |t: &'static str| {
            div()
                .w(px(30.))
                .text_size(px(9.))
                .text_color(theme::faint())
                .child(t)
        };

        div()
            .absolute()
            .top(px(20.))
            .left(px(-4.))
            .w(px(232.))
            .p(px(10.))
            .rounded(px(8.))
            .bg(theme::surface())
            .border_1()
            .border_color(theme::border())
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
                            .text_color(theme::ink())
                            .child(format!("Colour of {name}")),
                    )
                    .child(
                        div()
                            .id("picker-done")
                            .text_size(px(10.5))
                            .text_color(theme::muted())
                            .cursor_pointer()
                            .hover(|s| s.text_color(theme::ink()))
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
                            .border_color(theme::border()),
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
                    .bg(theme::sunken())
                    .child(div().h_full().w(relative(t)).rounded_full().bg(theme::muted()))
                    // The default, so it can be found by eye.
                    .when(knob == Knob::Size, |d| {
                        d.child(
                            div()
                                .absolute()
                                .top(px(-2.))
                                .left(relative(home))
                                .w(px(1.))
                                .h(px(7.))
                                .bg(theme::faint()),
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
                    .bg(theme::ink()),
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
    /// Where what is shown came from, as a quiet line under the table:
    /// the entry's page in its database, or the file in its folder.
    fn source_line(&self, structure: &Structure) -> Option<gpui::AnyElement> {
        let (text, go): (String, Box<dyn Fn(&mut App)>) = match structure.page() {
            Some(url) => (
                format!("View on {}", structure.source.label()),
                Box::new(move |cx| cx.open_url(&url)),
            ),
            None if structure.source == Source::File => {
                let path = PathBuf::from(self.origin.as_deref()?);
                let place = if cfg!(target_os = "macos") {
                    "Show in Finder"
                } else if cfg!(target_os = "windows") {
                    "Show in Explorer"
                } else {
                    "Show in folder"
                };
                (place.to_string(), Box::new(move |cx| cx.reveal_path(&path)))
            }
            None => return None,
        };
        Some(
            div()
                .pt(px(14.))
                .flex()
                .child(
                    div()
                        .id("source-link")
                        .debug_selector(|| "source-link".into())
                        .text_size(px(10.5))
                        .text_color(theme::faint())
                        .cursor_pointer()
                        .hover(|s| s.text_color(theme::text()))
                        .child(text)
                        .on_click(move |_, _, cx| go(cx)),
                )
                .into_any_element(),
        )
    }

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
                    .text_color(theme::faint())
                    .child(spaced(&category.to_uppercase())),
            )
            .child(
                div()
                    .font_family(theme::SERIF)
                    .text_size(px(27.))
                    .line_height(px(31.))
                    .text_color(theme::ink())
                    .overflow_hidden()
                    .child(name),
            )
            .when_some(blurb, |d, blurb| {
                d.child(
                    div()
                        .pt(px(10.))
                        .text_size(px(11.5))
                        .line_height(px(18.))
                        .text_color(theme::text())
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
            .children(self.source_line(structure))
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
                    .text_color(theme::faint())
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
                    .text_color(theme::ink())
                    .child(name),
            )
            .when_some(blurb, |d, blurb| {
                d.child(
                    div()
                        .pt(px(10.))
                        .text_size(px(11.5))
                        .line_height(px(18.))
                        .text_color(theme::text())
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
                            .border_color(theme::border())
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
                                        d.text_color(theme::text()).bg(theme::sunken())
                                    })
                                    .when(!self.copied, |d| {
                                        d.text_color(theme::muted())
                                            .hover(|s| s.bg(theme::sunken()).text_color(theme::text()))
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
                    .text_color(theme::faint())
                    .child(preview),
            )
            .when(structure.frame_count() > 1, |d| {
                d.child(
                    div()
                        .pt(px(8.))
                        .text_size(px(10.5))
                        .text_color(theme::faint())
                        .child(format!("{} frames", thousands(structure.frame_count()))),
                )
            })
            .children(self.source_line(structure))
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
                        .children(self.empty_message().into_iter().enumerate().map(|(i, line)| {
                            div()
                                .text_size(px(if i == 0 { 12. } else { 10.5 }))
                                .text_color(theme::faint())
                                .child(line)
                        })),
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
                            .text_color(theme::faint())
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
                        .text_color(theme::faint())
                        .child("drag to orbit · shift-drag to pan · scroll to zoom · double-click to reset"),
                )
            })
    }
}

impl Viewer {
    /// What an empty viewport says.
    fn empty_message(&self) -> [&'static str; 3] {
        match self.mode {
            Mode::View => [
                "Nothing loaded",
                "Search for a protein or molecule by name, paste an identifier,",
                if cfg!(target_os = "macos") {
                    "drop a structure file here, or open one with ⌘O"
                } else {
                    "drop a structure file here, or open one with Ctrl+O"
                },
            ],
            Mode::Fold => [
                "Nothing folded yet",
                "Paste a protein sequence, or a PDB ID or UniProt accession",
                "to fold its chains. The prediction appears here.",
            ],
        }
    }

    /// Stop trajectory playback, as when the tab goes behind another: a
    /// hidden tab has no one to play to.
    pub fn pause(&mut self, cx: &mut Context<Self>) {
        if self.playing {
            self.stop_playing();
            cx.notify();
        }
    }

    pub fn mode(&self) -> Mode {
        self.mode
    }

    /// Whether the tab has nothing in it and is doing nothing, so that
    /// opening something can use it rather than a new one.
    pub fn is_blank(&self) -> bool {
        self.structure.is_none() && self.busy.is_none()
    }

    /// Whether a job is running, and how far along it is when it knows.
    pub fn activity(&self) -> Option<Option<f32>> {
        self.busy
            .as_ref()
            .map(|_| self.progress.as_ref().and_then(|p| p.fraction))
    }

    /// The tab's name: what is loaded, or what the tab is for.
    pub fn title(&self) -> String {
        let loaded = self.structure.as_ref().map(|s| match s.source {
            crate::protein::Source::File | crate::protein::Source::Sample
                if !s.has_polymer() && s.title.trim().is_empty() =>
            {
                s.id.clone()
            }
            _ if s.has_polymer() => headline(s).0,
            _ if s.title.trim().is_empty() => s.id.clone(),
            _ => s.title.trim().to_string(),
        });
        match (self.mode, loaded) {
            (Mode::View, Some(name)) => name,
            (Mode::View, None) if self.busy.is_some() => "Loading…".into(),
            (Mode::View, None) => "New tab".into(),
            (Mode::Fold, Some(name)) => format!("Fold · {name}"),
            (Mode::Fold, None) => "Fold".into(),
        }
    }

    /// Put the caret in the search field.
    pub fn focus_input(&self, window: &mut Window, cx: &mut Context<Self>) {
        window.focus(&self.input.read(cx).focus_handle(cx));
    }

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
            .bg(theme::surface())
            .border_1()
            .border_color(theme::border())
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
                            .bg(theme::sunken())
                            .child(
                                div()
                                    .h_full()
                                    .w(relative(t))
                                    .rounded_full()
                                    .bg(theme::ink()),
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
                    .text_color(theme::muted())
                    .child(format!("{} / {}", frame + 1, count)),
            )
    }
}

impl Render for Viewer {
    fn render(&mut self, window: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        if viewer::spinning() {
            self.camera.rotate(0.02, 0.0);
            window.request_animation_frame();
        }
        let _timer = viewer::Phase::Tree.timer();
        let focused = self.input.read(cx).focus_handle(cx).is_focused(window);
        let folding = self.mode == Mode::Fold;
        let structure = self.structure.clone();

        div()
            .flex()
            .flex_col()
            .size_full()
            .bg(theme::paper())
            .text_color(theme::text())
            .track_focus(&self.focus)
            .on_mouse_move(cx.listener(Self::on_mouse_move))
            .on_mouse_up(MouseButton::Left, cx.listener(Self::on_mouse_up))
            .on_mouse_up_out(MouseButton::Left, cx.listener(Self::on_mouse_up))
            .on_mouse_up(MouseButton::Right, cx.listener(Self::on_mouse_up))
            // A fold tab has nothing in its toolbar until it has a prediction.
            .when(!folding || structure.is_some(), |d| d.child(self.toolbar(cx)))
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
                    .bg(theme::surface())
                    .border_r_1()
                    .border_color(theme::border())
                    .child(self.masthead())
                    .child(self.search(focused, cx))
                    // Everything under the field scrolls when the window is
                    // too short for it; the field itself stays put.
                    .child(
                        div()
                            .id("panel")
                            .flex_1()
                            .min_h_0()
                            .overflow_y_scroll()
                            .child(
                                div()
                                    .flex()
                                    .flex_col()
                                    .gap(px(16.))
                                    .pb(px(20.))
                                    // Model, steps and pLDDT: the settings for Fold and Compare.
                                    .when(folding, |d| d.children(self.plugin_options(cx)))
                                    .child(self.status(cx))
                                    .when(structure.is_none(), |d| d.child(self.examples(cx)))
                                    .when_some(structure, |d, s| d.child(self.info(&s, cx))),
                            ),
                    ),
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

    fn open(cx: &mut TestAppContext) -> (Entity<Viewer>, &mut VisualTestContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
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

    fn colors(view: &Entity<Viewer>, cx: &mut VisualTestContext) -> Vec<[f32; 3]> {
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
    fn palette_menu_columns_end_level(cx: &mut TestAppContext) {
        let (_view, cx) = open(cx);
        click(cx, "palette-button");
        let mut bounds = |selector: &'static str| {
            cx.debug_bounds(selector)
                .unwrap_or_else(|| panic!("{selector} is not in the menu"))
        };
        let (left, right) = (bounds("palette-Rocket"), bounds("palette-Greyscale"));
        let gap = (f32::from(left.bottom()) - f32::from(right.bottom())).abs();
        assert!(gap <= 1.0, "the columns end {gap}px apart");
        // The rainbows share the top row.
        let (spectrum, turbo, viridis) = (
            bounds("palette-Spectrum"),
            bounds("palette-Turbo"),
            bounds("palette-Viridis"),
        );
        assert_eq!(spectrum.top(), turbo.top());
        assert!(turbo.bottom() <= viridis.top());
    }

    #[gpui::test]
    fn camera_menu_resets_the_view(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        view.update(cx, |this, cx| {
            this.camera.rotate(1.0, 0.0);
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
            // The inked styles tone colours; the scene holds the toned one.
            view.update(cx, |this, _| this.style.inked())
        });
    }

    fn hit(kind: Kind, id: &str, title: &str) -> Hit {
        Hit::for_test(kind, id, title)
    }

    #[gpui::test]
    fn results_list_opens_the_chosen_hit(cx: &mut TestAppContext) {
        cx.update(crate::ui::input::bind_keys);
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        cx.simulate_input("green fluorescent protein");
        let (detected, searching) =
            view.update(cx, |this, _| (this.detected.clone(), this.searching.clone()));
        assert_eq!(detected, Input::Query("green fluorescent protein".into()));
        assert_eq!(searching.as_deref(), Some("green fluorescent protein"), "no search was started");

        // Stand in for the network: results arrive for what is typed.
        view.update(cx, |this, cx| {
            this.search_task = None;
            this.searching = None;
            this.hits = vec![
                hit(Kind::Pdb, "1GFL", "Structure of green fluorescent protein"),
                hit(Kind::AlphaFold, "P42212", "Green fluorescent protein"),
                hit(Kind::Molecule, "2519", "Caffeine"),
            ];
            this.hits_for = "green fluorescent protein".into();
            this.dropdown = true;
            cx.notify();
        });
        cx.run_until_parked();
        assert!(cx.debug_bounds("hit-2").is_some(), "the results list is not on screen");

        // Down twice, up once: the AlphaFold model.
        cx.simulate_keystrokes("down down up");
        assert_eq!(view.update(cx, |this, _| this.highlight), 1);
        cx.simulate_keystrokes("escape");
        assert!(!view.update(cx, |this, _| this.dropdown), "escape left the list up");
        cx.simulate_keystrokes("down");
        assert!(view.update(cx, |this, _| this.dropdown), "down did not bring the list back");

        click(cx, "hit-0");
        let busy = view.update(cx, |this, _| (this.dropdown, this.busy.is_some(), this.log.last().cloned()));
        assert!(!busy.0, "the list stayed up after a pick");
        assert!(busy.1, "picking a result started nothing");
        assert!(busy.2.is_some_and(|l| l.contains("1GFL")));
    }

    #[gpui::test]
    fn a_download_says_where_it_came_from(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        // A sample built into the app has no page to go to.
        assert!(cx.debug_bounds("source-link").is_none());
        view.update(cx, |this, cx| {
            let mut s = formats::load_bytes(SAMPLES[0].file, SAMPLES[0].bytes.to_vec()).unwrap();
            s.source = Source::PubChem;
            s.id = "CID 2519".into();
            this.show(s, cx);
        });
        cx.run_until_parked();
        assert!(cx.debug_bounds("source-link").is_some(), "no link to PubChem");
        let page = view.update(cx, |this, _| this.structure.as_ref().unwrap().page());
        assert_eq!(page.as_deref(), Some("https://pubchem.ncbi.nlm.nih.gov/compound/2519"));
    }

    #[gpui::test]
    fn enter_ahead_of_the_list_opens_only_an_exact_name(cx: &mut TestAppContext) {
        cx.update(crate::ui::input::bind_keys);
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        // Enter, with the search the field started still out: the network's
        // reply is stood in for below.
        let enter = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, cx| {
                this.submit(cx);
                this.search_task = None;
                assert!(this.open_when_found, "enter did not wait for the list");
            })
        };
        cx.simulate_input("caffeine");
        enter(&view, cx);
        // Nothing is called that: the list comes up and nothing opens.
        view.update(cx, |this, cx| {
            let near = vec![hit(Kind::Molecule, "1", "Caffeine citrate")];
            this.receive("caffeine".into(), Ok(near), cx);
        });
        let (busy, dropdown) = view.update(cx, |this, _| (this.busy.is_some(), this.dropdown));
        assert!(!busy, "a near miss was opened");
        assert!(dropdown, "the list did not come up");

        // Typed again with the real thing among the hits: that one opens.
        cx.simulate_keystrokes("backspace");
        cx.simulate_input("e");
        enter(&view, cx);
        view.update(cx, |this, cx| {
            let hits = vec![
                hit(Kind::Molecule, "1", "Caffeine citrate"),
                hit(Kind::Molecule, "2519", "Caffeine"),
            ];
            this.receive("caffeine".into(), Ok(hits), cx);
        });
        let (busy, last) = view.update(cx, |this, _| (this.busy.is_some(), this.log.last().cloned()));
        assert!(busy, "the exact name did not open");
        assert!(last.is_some_and(|l| l.contains("CID 2519")));
    }

    #[gpui::test]
    fn cmd_or_ctrl_backspace_clears_words(cx: &mut TestAppContext) {
        cx.update(crate::ui::input::bind_keys);
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        let text = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, cx| this.input.read(cx).text().to_string())
        };
        cx.simulate_input("green fluorescent");
        cx.simulate_keystrokes("shift-backspace");
        assert_eq!(text(&view, cx), "green fluorescen");
        #[cfg(target_os = "macos")]
        {
            cx.simulate_keystrokes("alt-backspace");
            assert_eq!(text(&view, cx), "green ");
            cx.simulate_keystrokes("cmd-backspace");
            assert_eq!(text(&view, cx), "");
        }
        #[cfg(not(target_os = "macos"))]
        {
            cx.simulate_keystrokes("ctrl-backspace");
            assert_eq!(text(&view, cx), "green ");
            cx.simulate_keystrokes("ctrl-backspace");
            assert_eq!(text(&view, cx), "");
        }
    }

    #[gpui::test]
    fn a_fold_tab_takes_sequences_and_never_searches(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::Fold, window, cx));
        cx.simulate_input("green fluorescent protein");
        let (detected, searching) =
            view.update(cx, |this, _| (this.detected.clone(), this.searching.clone()));
        assert!(matches!(detected, Input::Sequence(_)));
        assert!(searching.is_none(), "a fold tab searched");
    }

    #[gpui::test]
    fn a_viewer_never_offers_to_fold(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        cx.simulate_input("1CRN");
        assert_eq!(view.update(cx, |this, _| this.detected.clone()), Input::Pdb("1CRN".into()));
        assert!(cx.debug_bounds("compare").is_none());
        assert!(cx.debug_bounds("model").is_none());
    }
}
