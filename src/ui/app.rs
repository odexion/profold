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
    canvas, deferred, div, prelude::*, px, relative, Animation, AnimationExt, App, Bounds,
    ClipboardItem, Context, Entity, FocusHandle, Focusable, MouseButton, MouseDownEvent,
    MouseMoveEvent, MouseUpEvent, PathPromptOptions, Pixels, Point, ScrollDelta, ScrollWheelEvent,
    SharedString, Task, Window,
};

use crate::align;
use crate::contacts;
use crate::fetch::{self, Input, Retrieval};
use crate::formats;
use crate::formats::write;
use crate::job::{self, Event, Outcome};
use crate::select::Selection;

#[path = "console.rs"]
mod console;
#[path = "export.rs"]
mod export;
#[path = "session.rs"]
mod session;
use crate::measure::{self, Measurement};
use crate::protein::{LigandPaint, Source, Structure};
use crate::render::{
    self, Appearance, Camera, Coloring, Look, Palette, PaletteGroup, Picked, Repr, Scene, Species,
    Style,
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
    Hex(Pick),
    Rgb(Pick),
    /// A value in the film dialog.
    Movie(export::Knob),
}

/// What a colour picker colours.
#[derive(Clone, PartialEq)]
enum Pick {
    Species(Species),
    /// A chain, by index.
    Chain(usize),
    /// A kind of ligand, by residue name.
    Ligand(String),
}

/// The colour a ligand takes when it is first given one: gold, which sits
/// apart from the standard element colours and most chains.
const LIGAND_GOLD: [f32; 3] = SWATCHES[2];

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

/// An atom or residue picked out in the viewport.
#[derive(Clone, Copy, PartialEq, Debug)]
struct Target {
    /// Which structure: 0 is the one loaded, 1 the prediction laid over it.
    which: usize,
    picked: Picked,
}

/// A ligand's binding site: the ligand, by residue, and the residues of
/// the polymer around it.
#[derive(Clone, PartialEq, Debug)]
struct Site {
    ligand: usize,
    residues: Vec<usize>,
    /// What holds the ligand there, and what bumps into it.
    contacts: Vec<contacts::Contact>,
    /// Kinds of contact whose lines are not drawn.
    hidden: Vec<contacts::Kind>,
}

/// How close a residue comes to a ligand to count as its binding site, in
/// ångström: a contact, or a water's width from one.
const SITE_CUTOFF: f32 = 5.0;

/// A selection, drawn its own way over the rest.
#[derive(Clone)]
struct Layer {
    selection: Selection,
    /// Its atoms in the structure shown, found anew whenever that changes.
    atoms: Arc<Vec<bool>>,
    style: render::LayerStyle,
    /// One colour for all of it, or its elements' own.
    color: Option<[f32; 3]>,
}

/// The colours a layer steps through as its dot is clicked: its elements'
/// own first.
const LAYER_COLORS: [Option<[f32; 3]>; 6] = [
    None,
    Some([0.98, 0.78, 0.20]),
    Some([0.30, 0.70, 0.90]),
    Some([0.90, 0.40, 0.65]),
    Some([0.55, 0.78, 0.25]),
    Some([0.93, 0.93, 0.93]),
];

/// What the plots work out once per structure and frame.
#[derive(Default)]
struct PlotCache {
    /// The structure and frame it is for, by address and number.
    key: Option<(usize, usize)>,
    torsions: Arc<Vec<Option<(f32, f32)>>>,
    pae: Option<Arc<gpui::RenderImage>>,
}

/// The side of each plot, in pixels: the panel's width inside its margins.
const PLOT: f32 = PANEL_WIDTH - 40.0;

/// A small picture from colours row by row, as GPUI takes them.
fn picture(
    width: usize,
    height: usize,
    rgb: impl Fn(usize, usize) -> [f32; 3],
) -> Arc<gpui::RenderImage> {
    let mut bgra = vec![0u8; width * height * 4];
    for y in 0..height {
        for x in 0..width {
            let c = rgb(x, y);
            let o = (y * width + x) * 4;
            bgra[o] = (c[2].clamp(0.0, 1.0) * 255.0) as u8;
            bgra[o + 1] = (c[1].clamp(0.0, 1.0) * 255.0) as u8;
            bgra[o + 2] = (c[0].clamp(0.0, 1.0) * 255.0) as u8;
            bgra[o + 3] = 255;
        }
    }
    let buffer =
        image::RgbaImage::from_raw(width as u32, height as u32, bgra).expect("sized to fit");
    Arc::new(gpui::RenderImage::new(smallvec::smallvec![
        image::Frame::new(buffer)
    ]))
}

/// An analysis plot as a picture `w` by `h`: a line through the values on
/// faint rules at their least and greatest, or secondary structure, frames
/// across and residues down.
fn analysis_picture(
    plot: &crate::analysis::Plot,
    frames: usize,
    w: usize,
    h: usize,
    dark: bool,
) -> Arc<gpui::RenderImage> {
    use crate::analysis::Plot;
    let page = theme::view();
    let (xs, ys): (Vec<f32>, Vec<f32>) = match plot {
        Plot::Secondary {
            frames: taken,
            residues,
        } => {
            let rows = residues.len().max(1);
            let columns = taken.len().max(1);
            return picture(w, h, |x, y| {
                let (_, ss) = &taken[(x * columns / w).min(columns - 1)];
                match ss.get((y * rows / h).min(rows - 1)) {
                    Some(crate::protein::Ss::Coil) | None => {
                        page.map(|c| c * 0.92 + if dark { 0.05 } else { 0.0 })
                    }
                    Some(s) => render::SS_COLORS
                        .iter()
                        .find(|c| c.0 == *s)
                        .map_or(page, |c| c.2),
                }
            });
        }
        Plot::Time(points) => points
            .iter()
            .map(|(f, v)| (*f as f32 / (frames - 1).max(1) as f32, *v))
            .unzip(),
        Plot::Along(values) => {
            let n = values.len().max(1);
            values
                .iter()
                .enumerate()
                .map(|(i, v)| ((i as f32 + 0.5) / n as f32, *v))
                .unzip()
        }
    };
    let (lo, hi) = plot_range(&ys);
    let to_y = |v: f32| {
        ((1.0 - PLOT_PAD) - (v - lo) / (hi - lo) * (1.0 - 2.0 * PLOT_PAD)) * (h - 1) as f32
    };
    // Each column's height on the line, joined to the next so steep
    // stretches stay unbroken.
    let column: Vec<f32> = (0..w)
        .map(|x| {
            let t = x as f32 / (w - 1).max(1) as f32;
            let k = xs.partition_point(|&q| q < t);
            let v = match (
                k.checked_sub(1).map(|i| (xs[i], ys[i])),
                xs.get(k).map(|&q| (q, ys[k])),
            ) {
                (Some((x0, y0)), Some((x1, y1))) if x1 > x0 => {
                    y0 + (y1 - y0) * (t - x0) / (x1 - x0)
                }
                (Some((_, y)), _) | (None, Some((_, y))) => y,
                _ => lo,
            };
            to_y(v)
        })
        .collect();
    let line = if dark {
        [0.62, 0.74, 0.95]
    } else {
        [0.16, 0.36, 0.72]
    };
    let rule = page.map(|c| c * 0.9 + if dark { 0.06 } else { 0.0 });
    picture(w, h, |x, y| {
        let y = y as f32;
        let (a, b) = (column[x], column[(x + 1).min(w - 1)]);
        if y >= a.min(b) - 1.5 && y <= a.max(b) + 1.5 {
            line
        } else if (y - to_y(lo)).abs() < 1.0 || (y - to_y(hi)).abs() < 1.0 {
            rule
        } else {
            page
        }
    })
}

/// How far in from top and bottom the rules at a plot's least and greatest
/// values are, as a share of its height.
const PLOT_PAD: f32 = 0.1;

/// The least and greatest of `ys`, spread apart when they are one.
fn plot_range(ys: &[f32]) -> (f32, f32) {
    let lo = ys
        .iter()
        .cloned()
        .filter(|v| v.is_finite())
        .fold(f32::MAX, f32::min);
    let hi = ys
        .iter()
        .cloned()
        .filter(|v| v.is_finite())
        .fold(f32::MIN, f32::max);
    if lo > hi {
        (0.0, 1.0)
    } else if hi - lo < 1e-6 {
        (lo - 0.5, hi + 0.5)
    } else {
        (lo, hi)
    }
}

/// Where protein backbones are found on the Ramachandran plot: the right-
/// and left-handed helices and the strands, as soft regions to read the
/// dots against.
fn ramachandran_ground() -> Arc<gpui::RenderImage> {
    // One for each page, light and dark, made when first wanted.
    static GROUND: [std::sync::OnceLock<Arc<gpui::RenderImage>>; 2] =
        [const { std::sync::OnceLock::new() }; 2];
    let dark = theme::dark();
    GROUND[dark as usize].get_or_init(|| ground(dark)).clone()
}

fn ground(dark: bool) -> Arc<gpui::RenderImage> {
    // Centres and spreads in degrees, as (phi, psi, sigma phi, sigma psi, weight).
    const REGIONS: [(f32, f32, f32, f32, f32); 5] = [
        (-63.0, -43.0, 22.0, 18.0, 1.0),
        (-120.0, 130.0, 35.0, 30.0, 0.9),
        (-75.0, 145.0, 15.0, 20.0, 0.7),
        (60.0, 45.0, 14.0, 16.0, 0.45),
        (-90.0, 0.0, 25.0, 25.0, 0.35),
    ];
    let page = theme::view();
    let tint = if dark {
        [0.32, 0.36, 0.44]
    } else {
        [0.78, 0.84, 0.92]
    };
    let n = 120;
    picture(n, n, |x, y| {
        let phi = -180.0 + (x as f32 + 0.5) / n as f32 * 360.0;
        let psi = 180.0 - (y as f32 + 0.5) / n as f32 * 360.0;
        let wrap = |d: f32| (d + 180.0).rem_euclid(360.0) - 180.0;
        let density: f32 = REGIONS
            .iter()
            .map(|&(p, q, sp, sq, w)| {
                w * (-0.5 * ((wrap(phi - p) / sp).powi(2) + (wrap(psi - q) / sq).powi(2))).exp()
            })
            .sum();
        // Two bands: favoured, and allowed.
        let k = if density > 0.5 {
            1.0
        } else if density > 0.08 {
            0.45
        } else {
            0.0
        };
        [0, 1, 2].map(|i| page[i] + (tint[i] - page[i]) * k)
    })
}

/// How a structure is shown, as undo steps back through it: everything a
/// change of style, colour, layer or measurement touches, and nothing of
/// where the camera is or what the pointer is over.
#[derive(Clone, PartialEq)]
struct ShownAs {
    repr: Repr,
    coloring: Coloring,
    style: Style,
    palette: Palette,
    standard_colors: bool,
    show_cell: bool,
    appearance: Appearance,
    layers: Vec<(String, render::LayerStyle, Option<[f32; 3]>)>,
    measurements: Vec<Measurement>,
    chain_colors: Vec<Option<[f32; 3]>>,
    ligand_paint: Vec<(String, LigandPaint)>,
}

/// Changes closer together than this are one step to undo, as a slider's
/// drag is.
const UNDO_PAUSE: Duration = Duration::from_millis(500);

/// The camera turning by itself.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Turn {
    /// Round and round about the screen's upright.
    Spin,
    /// Back and forth through a quarter turn, as on a turntable.
    Rock,
}

/// Frames the camera takes to fly to something, at about sixty a second.
const FLY_FRAMES: u32 = 20;

/// Pointers that move further than this between press and release, in
/// pixels, were dragging the view rather than clicking on it.
const CLICK_SLOP: f32 = 3.0;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Drag {
    Orbit,
    Pan,
    /// Turning in the screen's plane, about the middle of the view.
    Roll,
}

/// The menus in the top bar. One at a time.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Dropdown {
    Samples,
    Palette,
    Camera,
    Image,
    Film,
    /// Open, with the ways to lay a second structure over the first.
    Open,
    Save,
    /// The viewport's own choices, folded into menus in a narrow window.
    Repr,
    Style,
    Coloring,
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
    /// Where the left button went down on the viewport, while that may yet
    /// turn out to be a click.
    press: Option<Point<Pixels>>,
    /// What the pointer is over, and where it is in the viewport.
    hover: Option<(Target, Point<Pixels>)>,
    /// What was clicked on, which glows until something else is.
    selected: Option<Target>,
    /// The residue of the loaded structure the pointer is over in the
    /// sequence, which glows in the viewport while it is.
    seq_hover: Option<usize>,
    /// The scene without anything glowing in it, which a glow is laid over
    /// without building the rest again.
    base: Arc<Scene>,
    /// Whether clicks in the viewport pick atoms to measure between.
    measuring: bool,
    measure_kind: measure::Kind,
    /// The atoms picked towards the next measurement, all of structure
    /// `pending_which`.
    pending: Vec<usize>,
    pending_which: usize,
    measurements: Vec<Measurement>,
    /// The binding site being looked at, of a ligand in the loaded structure.
    site: Option<Site>,
    /// The camera on its way to something, while it is.
    fly_task: Option<Task<()>>,
    /// The structure as the file deposits it, while one of its assemblies
    /// is shown in its place; and which.
    deposited: Option<Arc<Structure>>,
    assembly: Option<usize>,
    /// The structure as opened, while it is shown wrapped, unwrapped or
    /// replicated across its periodic boundaries.
    unwrapped_from: Option<Arc<Structure>>,
    /// Which of those changes are on.
    boundary: Boundary,
    /// Text pinned to atoms, drawn over the view, images and films.
    pub(super) labels: Vec<Label>,
    /// How the structure's volume is drawn, and the surface last drawn for it.
    pub(crate) iso: Iso,
    iso_cache: Option<(Iso, usize, Arc<Scene>)>,
    /// Whether what is opened next is laid over the structure shown rather
    /// than taking its place.
    overlay_next: bool,
    /// The camera turning by itself, until the view is touched.
    turning: Option<Turn>,
    turn_task: Option<Task<()>>,
    /// Images are saved with nothing behind the molecule.
    transparent: bool,
    /// The still view, and saved images, are ray traced: shadows, ambient
    /// occlusion and smooth edges.
    pub(super) traced: bool,
    /// How far a traced image being saved has got: rows done, of how many.
    saving: Option<(Arc<std::sync::atomic::AtomicUsize>, usize)>,
    /// The keyboard shortcuts, laid over the viewport.
    shortcuts: bool,
    /// Selections drawn their own way, in the order made.
    layers: Vec<Layer>,
    /// Where a selection is typed, and what was wrong with the last one.
    select_input: Entity<TextInput>,
    select_error: Option<SharedString>,
    /// Viewpoints saved to come back to, in the order saved.
    views: Vec<Camera>,
    /// The console: whether it is open, where commands are typed, what came
    /// of them, the commands typed before, and a script's lines yet to run.
    console_open: bool,
    console_input: Entity<TextInput>,
    console_log: Vec<(bool, String)>,
    console_history: Vec<String>,
    console_back: usize,
    script: Vec<String>,
    /// What can be undone, the present last, and what was undone, to redo.
    history: Vec<ShownAs>,
    undone: Vec<ShownAs>,
    last_change: std::time::Instant,
    /// The plots: where each was last laid out, what the pointer is over in
    /// it, and what is worked out once per structure.
    rama_bounds: Rc<Shared<Bounds<Pixels>>>,
    rama_hover: Option<usize>,
    /// The trajectory's analysis: whether its row is open, what is shown,
    /// the atoms it is taken of (pills chosen, or a selection typed in the
    /// console), what has been worked out, the work under way and how far
    /// it is, the plot drawn, and where the pointer is along it.
    pub(super) analysis_open: bool,
    analysis_pick: Option<crate::analysis::Measure>,
    analysis_pills: Option<Vec<bool>>,
    pub(super) analysis_custom: Option<String>,
    analysis: HashMap<(usize, crate::analysis::Measure, String), Arc<crate::analysis::Plot>>,
    analysis_task: Option<(
        (usize, crate::analysis::Measure, String),
        Arc<std::sync::atomic::AtomicBool>,
        Task<()>,
    )>,
    analysis_done: Option<(Arc<std::sync::atomic::AtomicUsize>, usize)>,
    analysis_image: RefCell<Option<((usize, bool), Arc<gpui::RenderImage>)>>,
    analysis_bounds: Rc<Shared<Bounds<Pixels>>>,
    /// The atoms the selection takes, found once for each structure and
    /// selection rather than at every repaint.
    analysis_atoms_found: RefCell<Option<((usize, usize, usize, String), Arc<Vec<usize>>)>>,
    analysis_hover: Option<f32>,
    pae_bounds: Rc<Shared<Bounds<Pixels>>>,
    pae_hover: Option<(usize, usize)>,
    /// Whether each plot is open under its row; both start closed.
    rama_open: bool,
    pae_open: bool,
    plot_cache: RefCell<PlotCache>,
    /// The triangles and capsules of `base` each structure was built into,
    /// by `Target::which`, for painting one of its residues.
    parts: [(std::ops::Range<usize>, std::ops::Range<usize>); 2],

    open: Option<Dropdown>,
    export: Option<Task<()>>,
    /// The open dialog, while it is up.
    prompt: Option<Task<()>>,
    /// Whether the periodic cell is drawn.
    show_cell: bool,
    /// Size and opacity per element or atom type.
    appearance: Appearance,
    /// Elements in their standard colours rather than the palette's.
    standard_colors: bool,
    /// The slider being dragged, if any.
    sliding: Option<(Species, Knob)>,
    /// Where each slider was last laid out.
    sliders: Rc<RefCell<HashMap<(Species, Knob), Bounds<Pixels>>>>,
    /// One text field, moved to whichever value is being typed.
    editor: Entity<TextInput>,
    editing: Option<Field>,
    /// What the open colour picker colours.
    picking: Option<Pick>,
    /// Every chain listed, rather than the first few.
    all_chains: bool,
    /// The settings as last written, so they are written only on change.
    saved: String,

    /// The film dialog, its preview and export, and recording.
    movie: export::Movie,

    /// Trajectory playback.
    playing: bool,
    play_task: Option<Task<()>>,
    /// Where the scrubber was last laid out, for turning a click into a frame.
    scrubber: Rc<Shared<Bounds<Pixels>>>,
    scrubbing: bool,
    /// The divider between the panel and the viewport is being dragged.
    resizing: bool,
    /// The window is too narrow for the panel's sections side by side, the
    /// angle readout or the hints: as narrow as a phone held upright.
    narrow: bool,
    /// The viewport is too narrow for its rows of pills, which become menus.
    compact: bool,

    #[cfg(feature = "simplefold")]
    pub(crate) simplefold: SimpleFold,

    copied: bool,
    copy_reset: Option<Task<()>>,
    /// What kept its own colour through the last change of palette, said
    /// for a moment beside the palette button.
    kept: Option<String>,
    kept_task: Option<Task<()>>,
    busy: Option<SharedString>,
    progress: Option<job::Progress>,
    log: Vec<SharedString>,
    error: Option<SharedString>,
    task: Option<Task<()>>,
}

impl Viewer {
    pub fn new(mode: Mode, window: &mut Window, cx: &mut Context<Self>) -> Self {
        let placeholder = match mode {
            Mode::View => "Search or paste anything",
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

        let console_input = cx.new(|cx| TextInput::new("a command, or help", cx));
        cx.subscribe_in(
            &console_input,
            window,
            |this, input, event, window, cx| match event {
                InputEvent::Submit => {
                    let line = input.read(cx).text().to_string();
                    input.update(cx, |i, cx| i.set_text("", cx));
                    if !line.trim().is_empty() {
                        this.console_history.push(line.clone());
                    }
                    this.console_back = 0;
                    let result = this.run_command(&line, cx);
                    this.say(&line, result);
                    cx.notify();
                }
                InputEvent::Up | InputEvent::Down => {
                    // Back through what was typed, and forward again.
                    let n = this.console_history.len();
                    this.console_back = match event {
                        InputEvent::Up => (this.console_back + 1).min(n),
                        _ => this.console_back.saturating_sub(1),
                    };
                    let text = if this.console_back == 0 {
                        String::new()
                    } else {
                        this.console_history[n - this.console_back].clone()
                    };
                    input.update(cx, |i, cx| i.set_text(text, cx));
                }
                InputEvent::Cancel => {
                    this.console_open = false;
                    window.focus(&this.focus);
                    cx.notify();
                }
                InputEvent::Changed => {}
            },
        )
        .detach();

        let select_input = cx.new(|cx| TextInput::new("e.g. resn HIS and within 5 of ligand", cx));
        cx.subscribe_in(
            &select_input,
            window,
            |this, input, event, _window, cx| match event {
                InputEvent::Submit => {
                    let text = input.read(cx).text().to_string();
                    if this.add_layer(&text, cx) {
                        input.update(cx, |i, cx| i.set_text("", cx));
                    }
                }
                InputEvent::Changed => {
                    if this.select_error.take().is_some() {
                        cx.notify();
                    }
                }
                InputEvent::Cancel | InputEvent::Up | InputEvent::Down => {}
            },
        )
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
                perspective: false,
                slab: None,
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
            press: None,
            hover: None,
            selected: None,
            seq_hover: None,
            base: Arc::new(Scene::default()),
            measuring: false,
            measure_kind: measure::Kind::Distance,
            pending: Vec::new(),
            pending_which: 0,
            measurements: Vec::new(),
            site: None,
            fly_task: None,
            deposited: None,
            assembly: None,
            unwrapped_from: None,
            boundary: Boundary::default(),
            iso: Iso::default(),
            iso_cache: None,
            labels: Vec::new(),
            overlay_next: false,
            turning: None,
            turn_task: None,
            transparent: false,
            traced: false,
            saving: None,
            shortcuts: false,
            layers: Vec::new(),
            select_input,
            select_error: None,
            views: Vec::new(),
            console_open: false,
            console_input,
            console_log: Vec::new(),
            console_history: Vec::new(),
            console_back: 0,
            script: Vec::new(),
            history: Vec::new(),
            undone: Vec::new(),
            last_change: std::time::Instant::now(),
            rama_bounds: Rc::new(Shared::new(Bounds::default())),
            rama_hover: None,
            analysis_open: false,
            analysis_pick: None,
            analysis_pills: None,
            analysis_custom: None,
            analysis: HashMap::new(),
            analysis_task: None,
            analysis_done: None,
            analysis_image: RefCell::new(None),
            analysis_bounds: Rc::new(Shared::new(Bounds::default())),
            analysis_atoms_found: RefCell::new(None),
            analysis_hover: None,
            pae_bounds: Rc::new(Shared::new(Bounds::default())),
            pae_hover: None,
            rama_open: false,
            pae_open: false,
            plot_cache: RefCell::new(PlotCache::default()),
            parts: [(0..0, 0..0), (0..0, 0..0)],
            open: None,
            export: None,
            prompt: None,
            show_cell: true,
            appearance: settings.appearance.clone(),
            standard_colors: settings.standard_colors,
            sliding: None,
            sliders: Rc::new(RefCell::new(HashMap::new())),
            editor,
            editing: None,
            picking: None,
            all_chains: false,
            saved: settings.to_json(),
            movie: Default::default(),
            playing: false,
            play_task: None,
            scrubber: Rc::new(Shared::new(Bounds::default())),
            scrubbing: false,
            resizing: false,
            narrow: false,
            compact: false,
            #[cfg(feature = "simplefold")]
            simplefold: SimpleFold::new(),
            copied: false,
            copy_reset: None,
            kept: None,
            kept_task: None,
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
        // Nor are a crystal's waters, beside a protein; a simulation's
        // solvent is half of what it shows.
        self.appearance.hide_solvent = polymer;
        self.stop_playing();
        self.camera = Camera {
            perspective: self.camera.perspective,
            ..Camera::frame(&structure)
        };
        self.all_chains = false;
        self.iso = match &structure.volume {
            Some(v) => Iso {
                shown: true,
                level: v.default_level(),
                mesh: false,
            },
            None => Iso::default(),
        };
        self.structure = Some(Arc::new(structure));
        self.compare = None;
        self.fit = None;
        self.unwrapped_from = None;
        self.boundary = Boundary::default();
        self.labels.clear();
        self.deposited = None;
        self.assembly = None;
        self.layers.clear();
        self.history.clear();
        self.undone.clear();
        self.views.clear();
        self.selected = None;
        self.seq_hover = None;
        self.hover = None;
        self.pending.clear();
        self.measurements.clear();
        self.site = None;
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
        self.camera = Camera {
            perspective: self.camera.perspective,
            ..Camera::frame_all(&[&reference, &prediction])
        };
        self.structure = Some(Arc::new(reference));
        self.compare = Some(Arc::new(prediction));
        self.fit = Some((rmsd, matched));
        self.repr = Repr::Cartoon;
        self.stop_playing();
        self.show_reference = true;
        self.show_prediction = true;
        self.selected = None;
        self.seq_hover = None;
        self.hover = None;
        self.pending.clear();
        self.measurements.clear();
        self.site = None;
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

    /// Structure `which` of a `Target`, if it is drawn.
    fn shown(&self, which: usize) -> Option<&Arc<Structure>> {
        match which {
            0 if self.show_reference || self.compare.is_none() => self.structure.as_ref(),
            1 if self.show_prediction => self.compare.as_ref(),
            _ => None,
        }
    }

    /// The atom or residue drawn nearest the viewer under a point of the
    /// window, worked out against the frame last painted.
    fn pick_at(&self, pos: Point<Pixels>) -> Option<Target> {
        let (bounds, (width, height)) = {
            let cache = self.cache.borrow();
            (cache.bounds, cache.size)
        };
        if !bounds.contains(&pos) || f32::from(bounds.size.width) <= 0.0 {
            return None;
        }
        // The frame may be drawn at the display's scale; pick in its pixels.
        let k = width as f32 / f32::from(bounds.size.width);
        let x = f32::from(pos.x - bounds.origin.x) * k;
        let y = f32::from(pos.y - bounds.origin.y) * k;
        (0..2)
            .filter_map(|which| {
                let s = self.shown(which)?;
                let picked = render::pick(
                    s,
                    self.repr,
                    &self.appearance,
                    &self.camera,
                    width,
                    height,
                    x,
                    y,
                )?;
                Some(Target { which, picked })
            })
            .max_by(|a, b| a.picked.depth.total_cmp(&b.picked.depth))
    }

    /// Point at something, or at nothing.
    fn set_hover(&mut self, pos: Option<Point<Pixels>>, cx: &mut Context<Self>) {
        let origin = self.cache.borrow().bounds.origin;
        let hover = pos.and_then(|p| Some((self.pick_at(p)?, p - origin)));
        if hover != self.hover {
            self.hover = hover;
            cx.notify();
        }
    }

    /// The pointer is over a card floating on the viewport: the structure
    /// behind it is out of reach, so nothing there is pointed at. A drag
    /// still passes through, so the view keeps turning under the card.
    fn over_card(&mut self, _: &MouseMoveEvent, _: &mut Window, cx: &mut Context<Self>) {
        if self.drag.is_none() && self.sliding.is_none() && !self.scrubbing {
            self.set_hover(None, cx);
            cx.stop_propagation();
        }
    }

    /// Pick out an atom or residue, or with `None` let go of the one picked.
    fn select(&mut self, target: Option<Target>, cx: &mut Context<Self>) {
        // A second click on the same thing lets go of it.
        let target = target.filter(|t| {
            self.selected
                .is_none_or(|s| (s.which, s.picked.atom) != (t.which, t.picked.atom))
        });
        if target != self.selected {
            self.selected = target;
            self.compose(cx);
        }
    }

    /// Point at a residue in the sequence, or with `None` at none.
    fn set_seq_hover(&mut self, residue: Option<usize>, cx: &mut Context<Self>) {
        if residue != self.seq_hover {
            self.seq_hover = residue;
            self.compose(cx);
        }
    }

    /// The residue of the loaded structure that is picked, if one is.
    fn selected_residue(&self) -> Option<usize> {
        let t = self.selected.filter(|t| t.which == 0)?;
        let s = self.structure.as_ref()?;
        Some(s.atoms.get(t.picked.atom)?.residue as usize)
    }

    /// The residue of the loaded structure under the pointer in the viewport.
    fn hovered_residue(&self) -> Option<usize> {
        let (t, _) = self.hover.filter(|(t, _)| t.which == 0)?;
        let s = self.structure.as_ref()?;
        Some(s.atoms.get(t.picked.atom)?.residue as usize)
    }

    /// Draw again after the theme has changed: the scene's colours are
    /// settled for the page they sit on.
    pub fn retheme(&mut self, cx: &mut Context<Self>) {
        self.rebuild(cx);
    }

    /// The panel has moved: nothing to redraw but the window, and the
    /// layout to remember.
    pub fn relayout(&mut self, cx: &mut Context<Self>) {
        self.persist();
        cx.notify();
    }

    /// Show the structure another way. A backbone style is refused for a
    /// structure without a backbone, as its pill is greyed out.
    pub fn set_repr(&mut self, repr: Repr, cx: &mut Context<Self>) {
        let Some(s) = self.structure.as_ref() else {
            return;
        };
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
            (!all.is_empty()).then(|| Camera {
                perspective: self.camera.perspective,
                slab: self.camera.slab,
                ..Camera::frame_all(&all)
            })
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
        self.refresh_layers();
        let iso = self.iso_scene();
        let mut scene = Scene::default();
        let mut parts = [(0..0, 0..0), (0..0, 0..0)];
        // What the layers hide is the loaded structure's alone.
        let mut plain = self.appearance.clone();
        plain.hidden = None;
        for (which, s) in (0..2).filter_map(|w| Some((w, self.shown(w)?))) {
            let start = (scene.tris.len(), scene.capsules.len());
            let appearance = if which == 0 { &self.appearance } else { &plain };
            scene.merge(render::build_scene(
                s,
                self.repr,
                self.coloring,
                quality,
                self.style,
                self.palette,
                self.standard_colors,
                self.show_cell,
                appearance,
            ));
            if which == 0 {
                if let Some(iso) = &iso {
                    scene.merge((**iso).clone());
                }
                for layer in &self.layers {
                    scene.merge(render::build_layer(
                        s,
                        &layer.atoms,
                        layer.style,
                        layer.color,
                        self.style,
                        self.palette,
                        self.standard_colors,
                        &self.appearance,
                    ));
                }
            }
            parts[which] = (start.0..scene.tris.len(), start.1..scene.capsules.len());
        }
        self.parts = parts;
        self.base = Arc::new(scene);
        self.compose(cx);
        self.persist();
    }

    /// Find each layer's atoms in the structure shown, which a frame, an
    /// assembly or a new structure changes, and hide what is to be hidden.
    fn refresh_layers(&mut self) {
        let Some(s) = self.structure.clone() else {
            self.appearance.hidden = None;
            return;
        };
        let mut hidden = vec![false; s.atoms.len()];
        for layer in &mut self.layers {
            let atoms = layer.selection.atoms(&s);
            if layer.style == render::LayerStyle::Hidden {
                for (h, a) in hidden.iter_mut().zip(&atoms) {
                    *h |= *a;
                }
            }
            layer.atoms = Arc::new(atoms);
        }
        self.appearance.hidden = hidden.iter().any(|h| *h).then(|| Arc::new(hidden));
    }

    /// Add a selection, typed, as sticks over what is drawn. Says what was
    /// wrong with it if it could not be read, and whether it was added.
    fn add_layer(&mut self, text: &str, cx: &mut Context<Self>) -> bool {
        if text.trim().is_empty() || self.structure.is_none() {
            return false;
        }
        match Selection::parse(text) {
            Ok(selection) => {
                self.select_error = None;
                self.layers.push(Layer {
                    selection,
                    atoms: Arc::new(Vec::new()),
                    style: render::LayerStyle::Sticks,
                    color: None,
                });
                self.rebuild(cx);
                true
            }
            Err(e) => {
                self.select_error = Some(format!("{e}").into());
                cx.notify();
                false
            }
        }
    }

    /// Lay what glows over the scene: the atom or residue picked, and the
    /// residue pointed at in the sequence.
    fn compose(&mut self, cx: &mut Context<Self>) {
        let hovered = self.seq_hover.and_then(|r| {
            let s = self.structure.as_ref()?;
            Some(Target {
                which: 0,
                picked: residue_pick(s, r)?,
            })
        });
        // Pointed at, it is painted most of the way; picked, all the way.
        let glows: Vec<(Target, bool)> = hovered
            .map(|t| (t, true))
            .into_iter()
            .chain(self.selected.map(|t| (t, false)))
            .collect();
        let mut scene = Scene::clone(&self.base);
        for (t, faint) in glows {
            if let Some(s) = self.shown(t.which) {
                if t.picked.residue {
                    let (tris, capsules) = self.parts[t.which].clone();
                    let residue = s.atoms[t.picked.atom].residue;
                    scene.paint_residue(tris, capsules, residue, self.style, faint);
                }
                scene.merge(render::build_highlight(
                    s,
                    t.picked,
                    self.repr,
                    self.style,
                    &self.appearance,
                ));
            }
        }
        if let (Some(site), Some(s)) = (&self.site, self.shown(0)) {
            scene.merge(render::build_site(
                s,
                &site.residues,
                self.repr,
                self.style,
                self.palette,
                self.standard_colors,
                &self.appearance,
            ));
        }
        // The atoms picked towards a measurement glow; its lines are drawn
        // over the viewport, by `measure_lines`.
        if let Some(s) = self.shown(self.pending_which) {
            for &atom in &self.pending {
                let picked = Picked {
                    atom,
                    residue: false,
                    depth: 0.0,
                };
                scene.merge(render::build_highlight(
                    s,
                    picked,
                    self.repr,
                    self.style,
                    &self.appearance,
                ));
            }
        }
        self.scene = Arc::new(scene.stamp());
        self.remember();
        cx.notify();
    }

    /// How the structure is shown now.
    fn look_now(&self) -> Option<ShownAs> {
        let s = self.structure.as_ref()?;
        let mut appearance = self.appearance.clone();
        appearance.hidden = None;
        let mut ligand_paint: Vec<(String, LigandPaint)> = s
            .ligand_paint
            .iter()
            .map(|(k, v)| (k.clone(), *v))
            .collect();
        ligand_paint.sort_by(|a, b| a.0.cmp(&b.0));
        Some(ShownAs {
            repr: self.repr,
            coloring: self.coloring,
            style: self.style,
            palette: self.palette,
            standard_colors: self.standard_colors,
            show_cell: self.show_cell,
            appearance,
            layers: self
                .layers
                .iter()
                .map(|l| (l.selection.text().to_string(), l.style, l.color))
                .collect(),
            measurements: self.measurements.clone(),
            chain_colors: s.chains.iter().map(|c| c.color).collect(),
            ligand_paint,
        })
    }

    /// Keep how the structure is shown, if it changed, as a step to undo.
    fn remember(&mut self) {
        let Some(now) = self.look_now() else { return };
        if self.history.last() == Some(&now) {
            return;
        }
        // A burst of changes, a slider dragged, is one step.
        if self.history.len() > 1 && self.last_change.elapsed() < UNDO_PAUSE {
            self.history.pop();
        }
        self.history.push(now);
        self.last_change = std::time::Instant::now();
        self.undone.clear();
        if self.history.len() > 200 {
            self.history.remove(0);
        }
    }

    /// Show the structure as it was shown a step ago, or with `redo` a step
    /// on again.
    pub fn undo(&mut self, redo: bool, cx: &mut Context<Self>) {
        let target = if redo {
            let Some(next) = self.undone.pop() else {
                return;
            };
            self.history.push(next.clone());
            next
        } else {
            if self.history.len() < 2 {
                return;
            }
            let now = self.history.pop().unwrap();
            self.undone.push(now);
            self.history.last().unwrap().clone()
        };
        self.repr = target.repr;
        self.coloring = target.coloring;
        self.style = target.style;
        self.palette = target.palette;
        self.standard_colors = target.standard_colors;
        self.show_cell = target.show_cell;
        self.appearance = target.appearance.clone();
        self.layers = target
            .layers
            .iter()
            .filter_map(|(text, style, color)| {
                Some(Layer {
                    selection: Selection::parse(text).ok()?,
                    atoms: Arc::new(Vec::new()),
                    style: *style,
                    color: *color,
                })
            })
            .collect();
        self.measurements = target.measurements.clone();
        if let Some(s) = self.structure.as_mut() {
            let s = Arc::make_mut(s);
            for (c, color) in s.chains.iter_mut().zip(&target.chain_colors) {
                c.color = *color;
            }
            s.ligand_paint = target.ligand_paint.iter().cloned().collect();
        }
        // Not a new step: the look being shown is the one the history ends on.
        self.last_change = std::time::Instant::now() - UNDO_PAUSE;
        self.rebuild(cx);
    }

    /// Turn the camera's attention to a point: it slides across and zooms
    /// until a ball `extent` ångström across fills the view, and with `face`,
    /// turns too, until that direction in the structure points at the viewer.
    fn fly_to(
        &mut self,
        centre: [f32; 3],
        extent: f32,
        face: Option<[f32; 3]>,
        cx: &mut Context<Self>,
    ) {
        self.halt();
        let start = self.camera;
        let (from, zoom0, pan0) = (start.center, start.zoom, start.pan);
        let zoom1 = (start.radius / extent.max(2.0)).clamp(0.12, 24.0);
        let (axis, angle) = face.map_or(([0.0; 3], 0.0), |d| start.turn_to_face(d));
        self.fly_task = Some(cx.spawn(async move |this, cx| {
            for k in 1..=FLY_FRAMES {
                cx.background_executor()
                    .timer(Duration::from_millis(16))
                    .await;
                let t = k as f32 / FLY_FRAMES as f32;
                let e = t * t * (3.0 - 2.0 * t);
                let going = this.update(cx, |this, cx| {
                    let c = &mut this.camera;
                    c.orientation = start.turned(axis, angle * e);
                    c.center = [0, 1, 2].map(|i| from[i] + (centre[i] - from[i]) * e);
                    // Zoom is a ratio, so it is eased as one.
                    c.zoom = zoom0 * (zoom1 / zoom0).powf(e);
                    c.pan = [pan0[0] * (1.0 - e), pan0[1] * (1.0 - e)];
                    cx.notify();
                });
                if going.is_err() {
                    return;
                }
            }
        }));
    }

    /// Fly back to a viewpoint saved: its angle, centre, zoom and offset.
    fn fly_to_view(&mut self, view: Camera, cx: &mut Context<Self>) {
        self.halt();
        let start = self.camera;
        let (axis, angle) = start.turn_to(view.orientation);
        self.camera.perspective = view.perspective;
        self.camera.slab = view.slab;
        self.fly_task = Some(cx.spawn(async move |this, cx| {
            for k in 1..=FLY_FRAMES {
                cx.background_executor()
                    .timer(Duration::from_millis(16))
                    .await;
                let t = k as f32 / FLY_FRAMES as f32;
                let e = t * t * (3.0 - 2.0 * t);
                let going = this.update(cx, |this, cx| {
                    let c = &mut this.camera;
                    c.orientation = if k == FLY_FRAMES {
                        view.orientation
                    } else {
                        start.turned(axis, angle * e)
                    };
                    c.center =
                        [0, 1, 2].map(|i| start.center[i] + (view.center[i] - start.center[i]) * e);
                    c.zoom = start.zoom * (view.zoom / start.zoom).powf(e);
                    c.pan = [0, 1].map(|i| start.pan[i] + (view.pan[i] - start.pan[i]) * e);
                    cx.notify();
                });
                if going.is_err() {
                    return;
                }
            }
        }));
    }

    /// Fly to what is picked: an atom up close, a residue whole.
    pub fn focus_selected(&mut self, cx: &mut Context<Self>) {
        let Some(t) = self.selected else { return };
        let Some(s) = self.shown(t.which) else { return };
        let atoms = if t.picked.residue {
            s.residue_atoms(s.atoms[t.picked.atom].residue as usize)
        } else {
            vec![t.picked.atom]
        };
        let (centre, reach) = spread(atoms.iter().map(|&i| s.atoms[i].pos));
        self.fly_to(centre, reach * 2.0 + 10.0, None, cx);
    }

    /// Look at a ligand's binding site: its side chains drawn, each named,
    /// and the camera brought round to it.
    fn show_site(&mut self, ligand: usize, cx: &mut Context<Self>) {
        let Some(s) = self.structure.clone() else {
            return;
        };
        let residues = s.site(ligand, SITE_CUTOFF);
        let mut points: Vec<[f32; 3]> = s
            .residue_atoms(ligand)
            .iter()
            .map(|&i| s.atoms[i].pos)
            .collect();
        let (centre, _) = spread(points.iter().copied());
        points.extend(residues.iter().filter_map(|&r| s.side_chain_centre(r)));
        let found = contacts::between(&s, ligand, &residues);
        self.site = Some(Site {
            ligand,
            residues,
            contacts: found,
            hidden: Vec::new(),
        });
        self.compose(cx);
        // Round the ligand, so that it and its neighbours are all in view.
        let reach = points
            .iter()
            .map(|p| {
                ((p[0] - centre[0]).powi(2)
                    + (p[1] - centre[1]).powi(2)
                    + (p[2] - centre[2]).powi(2))
                .sqrt()
            })
            .fold(0.0, f32::max);
        // Brought round to the near side, out from the middle of the
        // structure, so the protein it sits in is not in the way.
        let (middle, _) = spread(s.atoms.iter().map(|a| a.pos));
        let out = [
            centre[0] - middle[0],
            centre[1] - middle[1],
            centre[2] - middle[2],
        ];
        let face = (out.iter().map(|v| v * v).sum::<f32>() > 1.0).then_some(out);
        self.fly_to(centre, reach * 2.0 + 4.0, face, cx);
    }

    /// The binding site of the next copy of a ligand, by name, from the one
    /// being looked at; the first, if it is another ligand's.
    fn next_site(&mut self, name: &str, cx: &mut Context<Self>) {
        let Some(s) = self.structure.as_ref() else {
            return;
        };
        let copies = s.copies(name);
        if copies.is_empty() {
            return;
        }
        let at = self
            .site
            .as_ref()
            .and_then(|site| copies.iter().position(|&r| r == site.ligand))
            .map_or(0, |i| (i + 1) % copies.len());
        self.show_site(copies[at], cx);
    }

    fn close_site(&mut self, cx: &mut Context<Self>) {
        if self.site.take().is_some() {
            self.compose(cx);
        }
    }

    /// Show one of the structure's assemblies in its place, or with `None`
    /// the structure as deposited. Style and colours stay; what was picked,
    /// measured or looked into goes, as its atoms are numbered anew.
    fn show_assembly(&mut self, k: Option<usize>, cx: &mut Context<Self>) {
        let Some(base) = self.deposited.clone().or_else(|| self.structure.clone()) else {
            return;
        };
        if k == self.assembly {
            return;
        }
        let shown = match k.and_then(|k| base.assemblies.get(k)) {
            Some(a) => {
                let mut built = crate::assembly::build(&base, a);
                if built.has_polymer() {
                    crate::dssp::assign(&mut built);
                }
                self.deposited = Some(base.clone());
                Arc::new(built)
            }
            None => {
                self.deposited = None;
                base
            }
        };
        self.assembly = k;
        self.camera = Camera {
            perspective: self.camera.perspective,
            ..Camera::frame(&shown)
        };
        self.structure = Some(shown);
        self.selected = None;
        self.seq_hover = None;
        self.hover = None;
        self.pending.clear();
        self.measurements.clear();
        self.site = None;
        self.rebuild(cx);
    }

    /// Start measuring, or stop. What was being picked is let go of either way.
    fn toggle_measuring(&mut self, cx: &mut Context<Self>) {
        self.measuring = !self.measuring;
        self.pending.clear();
        self.selected = None;
        self.compose(cx);
    }

    fn set_measure_kind(&mut self, kind: measure::Kind, cx: &mut Context<Self>) {
        self.measure_kind = kind;
        self.pending.truncate(kind.atoms() - 1);
        self.compose(cx);
    }

    /// A click while measuring: one more atom towards the measurement, which
    /// is made once there are enough. A click on nothing starts over.
    fn measure_click(&mut self, target: Option<Target>, cx: &mut Context<Self>) {
        match target {
            None => self.pending.clear(),
            Some(t) => {
                let atom = t.picked.atom;
                // One measurement is within one structure.
                if t.which != self.pending_which {
                    self.pending.clear();
                    self.pending_which = t.which;
                }
                if self.pending.last() == Some(&atom) {
                    return;
                }
                self.pending.push(atom);
                if self.pending.len() == self.measure_kind.atoms() {
                    self.measurements.push(Measurement {
                        kind: self.measure_kind,
                        which: t.which,
                        atoms: std::mem::take(&mut self.pending),
                    });
                }
            }
        }
        self.compose(cx);
    }

    fn remove_measurement(&mut self, i: usize, cx: &mut Context<Self>) {
        if i < self.measurements.len() {
            self.measurements.remove(i);
            self.compose(cx);
        }
    }

    /// Where a point of structure space is drawn in the viewport, from its
    /// top left corner, as the frame last painted puts it.
    fn project(&self, p: [f32; 3]) -> Option<Point<Pixels>> {
        let cache = self.cache.borrow();
        let (w, h) = cache.size;
        let width = f32::from(cache.bounds.size.width);
        if w == 0 || width <= 0.0 {
            return None;
        }
        let cam = &self.camera;
        let m = cam.basis();
        let q = [
            p[0] - cam.center[0],
            p[1] - cam.center[1],
            p[2] - cam.center[2],
        ];
        let d = |r: [f32; 3]| r[0] * q[0] + r[1] * q[1] + r[2] * q[2];
        let (v, _) = cam.warp([d(m[0]), d(m[1]), d(m[2])]);
        let scale = w.min(h) as f32 / (2.2 * cam.radius) * cam.zoom;
        let k = w as f32 / width;
        let x = (w as f32 * 0.5 + cam.pan[0] + v[0] * scale) / k;
        let y = (h as f32 * 0.5 + cam.pan[1] - v[1] * scale) / k;
        Some(gpui::point(px(x), px(y)))
    }

    fn finish(&mut self, result: anyhow::Result<Structure>, cx: &mut Context<Self>) {
        self.busy = None;
        self.progress = None;
        match result {
            // A trajectory with as many atoms as the structure shown moves
            // that structure's atoms, rather than standing in for it.
            Ok(s)
                if s.format == Some("trajectory")
                    && self.compare.is_none()
                    && self
                        .structure
                        .as_ref()
                        .is_some_and(|shown| shown.atoms.len() == s.atoms.len()) =>
            {
                self.opening = None;
                let frames: Vec<crate::protein::Snapshot> = if s.frames.is_empty() {
                    vec![crate::protein::Snapshot {
                        step: None,
                        pos: s.atoms.iter().map(|a| a.pos).collect(),
                        cell: s.cell,
                    }]
                } else {
                    s.frames.to_vec()
                };
                let n = frames.len();
                if let Some(moved) = self
                    .structure
                    .as_ref()
                    .and_then(|shown| shown.with_frames(frames))
                {
                    self.note(format!(
                        "{} frames from {} superimposed onto {}",
                        thousands(n),
                        s.id,
                        moved.id
                    ));
                    self.stop_playing();
                    self.structure = Some(Arc::new(moved));
                    self.rebuild(cx);
                }
            }
            // A map, which has no atoms, is laid over the structure shown.
            Ok(s)
                if s.format == Some("map")
                    && self.structure.as_ref().is_some_and(|s| !s.atoms.is_empty()) =>
            {
                self.opening = None;
                if let (Some(shown), Some(volume)) = (self.structure.as_mut(), s.volume) {
                    let [a, b, c] = volume.dims;
                    let said = format!("map {} · {a} × {b} × {c} laid over {}", s.id, shown.id);
                    self.iso = Iso {
                        shown: true,
                        level: volume.default_level(),
                        mesh: false,
                    };
                    Arc::make_mut(shown).volume = Some(volume);
                    self.note(said);
                    self.rebuild(cx);
                }
            }
            Ok(s) if std::mem::take(&mut self.overlay_next) => {
                // The tab still reopens as what it was opened as.
                self.opening = None;
                self.overlay_with(s, cx);
            }
            Ok(s) => {
                self.origin = self.opening.take();
                if let Some(o) = &self.origin {
                    crate::settings::remember(o);
                }
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
                let rest = std::mem::take(&mut self.script);
                if !rest.is_empty() {
                    self.run_script(rest, cx);
                }
            }
            Err(e) => {
                self.opening = None;
                self.overlay_next = false;
                self.script.clear();
                self.error = Some(format!("{e:#}").into());
                cx.notify();
            }
        }
    }

    /// Lay a structure over the one shown, paired residue by residue where
    /// their sequences align, and show the two as a comparison.
    fn overlay_with(&mut self, other: Structure, cx: &mut Context<Self>) {
        let Some(shown) = self.structure.as_ref() else {
            self.show(other, cx);
            return;
        };
        let mut reference = Structure::clone(shown);
        let mut mobile = other;
        match align::superpose(&mut reference, &mut mobile) {
            Ok(fit) => {
                self.note(format!(
                    "{} on {} · {:.2} Å over {} of {} paired, {}% identical",
                    mobile.id,
                    reference.id,
                    fit.rmsd,
                    fit.matched,
                    fit.aligned,
                    fit.identical * 100 / fit.aligned.max(1)
                ));
                let origin = self.origin.clone();
                self.show_comparison(reference, mobile, fit.rmsd, fit.matched, cx);
                self.origin = origin;
            }
            Err(e) => {
                self.error = Some(format!("{e:#}").into());
                cx.notify();
            }
        }
    }

    /// Take the overlaid structure away, and show the one under it alone.
    fn remove_overlay(&mut self, cx: &mut Context<Self>) {
        if self.compare.take().is_none() {
            return;
        }
        self.fit = None;
        if let Some(s) = self.structure.as_mut() {
            let s = Arc::make_mut(s);
            s.tint = None;
            s.deviation.clear();
        }
        self.show_reference = true;
        self.show_prediction = true;
        self.coloring = if self.structure.as_ref().is_some_and(|s| s.has_polymer()) {
            Coloring::Chain
        } else {
            Coloring::Element
        };
        self.selected = None;
        self.measurements.clear();
        self.pending.clear();
        self.rebuild(cx);
    }

    /// Pick a file to lay over the structure shown.
    pub fn prompt_overlay(&mut self, cx: &mut Context<Self>) {
        self.open = None;
        if self.structure.is_none() {
            return;
        }
        let chosen = cx.prompt_for_paths(PathPromptOptions {
            files: true,
            directories: false,
            multiple: false,
            prompt: Some("Overlay".into()),
        });
        self.prompt = Some(cx.spawn(async move |this, cx| {
            let Ok(Ok(Some(paths))) = chosen.await else {
                return;
            };
            if let Some(path) = paths.into_iter().next() {
                let _ = this.update(cx, |this, cx| {
                    this.overlay_next = true;
                    this.open_path(path, cx);
                });
            }
        }));
    }

    /// The open menu: a file in place of what is shown, or one laid over
    /// it, from a file or whatever is typed next into the field.
    fn open_menu(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let comparing = self.compare.is_some();
        menu(240.)
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.open = None;
                cx.notify();
            }))
            .child(
                item("open-file", "Open a file…".into(), false)
                    .debug_selector(|| "open-file".into())
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.open = None;
                        this.prompt_open(cx);
                    })),
            )
            .child(divide())
            .child(
                item("overlay-file", "Overlay a file…".into(), false)
                    .debug_selector(|| "overlay-file".into())
                    .on_click(cx.listener(|this, _, _, cx| this.prompt_overlay(cx))),
            )
            .child(
                item(
                    "overlay-entry",
                    "Overlay an entry by ID or name…".into(),
                    self.overlay_next,
                )
                .debug_selector(|| "overlay-entry".into())
                .on_click(cx.listener(|this, _, window, cx| {
                    this.open = None;
                    this.overlay_next = !this.overlay_next;
                    if this.overlay_next {
                        window.focus(&this.input.read(cx).focus_handle(cx));
                    }
                    cx.notify();
                })),
            )
            .when(comparing, |d| {
                d.child(divide()).child(
                    item("overlay-remove", "Remove the overlay".into(), false)
                        .debug_selector(|| "overlay-remove".into())
                        .on_click(cx.listener(|this, _, _, cx| {
                            this.open = None;
                            this.remove_overlay(cx);
                        })),
                )
            })
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
                let Some(query) = self.search_text() else {
                    return;
                };
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
            Input::UniProt(acc) => self
                .simplefold
                .fold_entry(Retrieval::AlphaFold(acc.clone()), acc),
            Input::Empty => return,
            _ => {
                self.error = Some(
                    "Folding takes a protein sequence, or a PDB ID or UniProt accession whose chains will be predicted."
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

        self.task = Some(cx.spawn(async move |this, cx| loop {
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
        if path.extension().is_some_and(|e| e == session::EXTENSION) {
            let absolute = std::fs::canonicalize(&path).unwrap_or_else(|_| path.clone());
            let opened = std::fs::read_to_string(&path)
                .map_err(anyhow::Error::from)
                .and_then(|text| self.open_session(&text, cx));
            match opened {
                Ok(()) => {
                    self.origin = Some(absolute.to_string_lossy().into_owned());
                    crate::settings::remember(&absolute.to_string_lossy());
                    self.note(format!("Opened the session {}", path.display()));
                }
                Err(e) => self.error = Some(format!("{e:#}").into()),
            }
            cx.notify();
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
        if self
            .structure
            .as_ref()
            .is_some_and(|s| s.frame + 1 >= count)
        {
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
            standard_colors: self.standard_colors,
            appearance: self.appearance.clone(),
            dark: theme::dark(),
            vertical: theme::vertical(),
            side_width: theme::side_width(),
            band_height: theme::band_height(),
        }
        .to_json();
        if json != self.saved {
            Settings::save(&json);
            self.saved = json;
        }
    }

    /// Put the field editor on a value, with what it shows selected so
    /// typing replaces it.
    fn start_edit(
        &mut self,
        field: Field,
        text: String,
        window: &mut Window,
        cx: &mut Context<Self>,
    ) {
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
            Field::Hex(target) => {
                if let Some(c) = settings::parse_hex(&text) {
                    self.paint(target, Some(c), cx);
                }
            }
            Field::Rgb(target) => {
                if let Some(c) = settings::parse_rgb(&text) {
                    self.paint(target, Some(c), cx);
                }
            }
            Field::Movie(knob) => self.commit_movie(knob, &text, cx),
        }
        cx.notify();
    }

    /// Colour what a picker is on, or `None` to give it back its own colour.
    fn paint(&mut self, target: Pick, color: Option<[f32; 3]>, cx: &mut Context<Self>) {
        match target {
            Pick::Species(species) => self.set_color(species, color, cx),
            Pick::Chain(i) => self.set_chain_color(i, color, cx),
            Pick::Ligand(name) => {
                let carbons = self.structure.as_ref().is_some_and(|s| {
                    s.atoms
                        .iter()
                        .any(|a| a.element.0 == 6 && s.residues[a.residue as usize].name == name)
                });
                let paint = color.map(|c| match self.ligand_paint(&name) {
                    Some(p) => p.with_color(c),
                    None if carbons => LigandPaint::Carbons(c),
                    None => LigandPaint::Solid(c),
                });
                self.set_ligand_paint(name, paint, cx)
            }
        }
    }

    fn ligand_paint(&self, name: &str) -> Option<LigandPaint> {
        self.structure.as_ref()?.ligand_paint.get(name).copied()
    }

    /// Paint a kind of ligand, or `None` to give it back its elements'
    /// colours. It shows under any scheme, so the scheme is left alone.
    /// Turn one of the periodic-boundary changes on, or off if it is on,
    /// and show the structure as opened with those now on, in one order:
    /// how atoms are placed, then frames fitted, then the cell repeated.
    /// `Periodic::Original` turns them all off. What it did, to say so.
    pub(super) fn periodic(
        &mut self,
        op: Periodic,
        cx: &mut Context<Self>,
    ) -> anyhow::Result<String> {
        let Some(shown) = self.structure.clone() else {
            anyhow::bail!("nothing is open")
        };
        let was = self.boundary;
        let mut now = was;
        let toggle = |place: Place| {
            if was.place == place {
                Place::AsRead
            } else {
                place
            }
        };
        match op {
            Periodic::Original => {
                if was == Boundary::default() {
                    anyhow::bail!("nothing to undo")
                }
                now = Boundary::default();
            }
            Periodic::Wrap => now.place = toggle(Place::Wrap),
            Periodic::Whole => now.place = toggle(Place::Whole),
            Periodic::InTime => now.place = toggle(Place::InTime),
            Periodic::Fit => now.fitted = !was.fitted,
            Periodic::Replicate(n) => {
                now.replicate = (n != [1, 1, 1] && was.replicate != Some(n)).then_some(n)
            }
        }
        let base = self.unwrapped_from.clone().unwrap_or(shown);
        let next = if now == Boundary::default() {
            base.clone()
        } else {
            let mut s = (*base).clone();
            match now.place {
                Place::AsRead => {}
                Place::Wrap => crate::pbc::wrap(&mut s)?,
                Place::Whole => crate::pbc::make_whole(&mut s)?,
                Place::InTime => crate::pbc::unwrap_in_time(&mut s)?,
            }
            if now.fitted {
                crate::analysis::fit_frames(&mut s)?;
            }
            if let Some(n) = now.replicate {
                s = crate::pbc::replicate(&s, n)?;
            }
            Arc::new(s)
        };
        let said = match op {
            Periodic::Original => "as opened".to_string(),
            Periodic::Replicate(n) => match now.replicate {
                Some(_) => format!("replicated {} × {} × {}", n[0], n[1], n[2]),
                None => "one cell".to_string(),
            },
            Periodic::Fit => format!(
                "frames {}",
                if now.fitted {
                    "fitted onto the first"
                } else {
                    "as read"
                }
            ),
            _ => match now.place {
                Place::AsRead => "atoms where the file puts them",
                Place::Wrap => "wrapped into the cell",
                Place::Whole => "molecules made whole",
                Place::InTime => "unwrapped through time",
            }
            .to_string(),
        };
        self.unwrapped_from = (now != Boundary::default()).then_some(base);
        self.boundary = now;
        if now.replicate != was.replicate {
            self.camera = Camera {
                perspective: self.camera.perspective,
                ..Camera::frame(&next)
            };
        }
        self.structure = Some(next);
        self.selected = None;
        self.hover = None;
        self.measurements.clear();
        self.rebuild(cx);
        // The frames have changed, so what was worked out from them has.
        if let (true, Some(m)) = (self.analysis_open, self.analysis_pick) {
            self.start_measure(m, cx);
        }
        Ok(said)
    }

    /// The volume's surface as drawn now, built when its level or look last
    /// changed: blue where the values pass the level, and for data with a
    /// sign, red where they pass its negative.
    pub(crate) fn iso_scene(&mut self) -> Option<Arc<Scene>> {
        let volume = self.structure.as_ref()?.volume.clone()?;
        if !self.iso.shown {
            return None;
        }
        let key = Arc::as_ptr(&volume) as usize;
        if let Some((iso, k, scene)) = &self.iso_cache {
            if *iso == self.iso && *k == key {
                return Some(scene.clone());
            }
        }
        let level = self.iso.level;
        let mut lobes = vec![(volume.isosurface(level, false), [0.22, 0.48, 0.86])];
        if volume.signed() {
            lobes.push((volume.isosurface(-level, true), [0.86, 0.30, 0.22]));
        }
        let step = volume
            .axes
            .iter()
            .map(|a| (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt())
            .fold(f32::MAX, f32::min);
        let mut scene = Scene::default();
        for (facets, color) in lobes {
            if self.iso.mesh {
                // Each edge once, as a thin rod.
                let key = |p: [f32; 3]| p.map(|x| (x * 1000.0).round() as i64);
                let mut seen = std::collections::HashSet::new();
                for f in &facets {
                    for (a, b) in [(0, 1), (1, 2), (2, 0)] {
                        let (ka, kb) = (key(f.p[a]), key(f.p[b]));
                        if seen.insert(if ka < kb { (ka, kb) } else { (kb, ka) }) {
                            scene.capsules.push(render::Capsule {
                                a: f.p[a],
                                b: f.p[b],
                                r: (step * 0.05).max(0.01),
                                ca: color,
                                cb: color,
                                alpha: 1.0,
                            });
                        }
                    }
                }
            } else {
                scene.tris.extend(facets.iter().map(|f| render::Tri {
                    p: f.p,
                    n: f.n,
                    c: color,
                }));
            }
        }
        let scene = Arc::new(scene);
        self.iso_cache = Some((self.iso, key, scene.clone()));
        Some(scene)
    }

    /// The atoms an analysis can be taken of, as pills, each with the
    /// selection it stands for: alpha carbons and the backbone of a protein,
    /// then each atom type of a simulation, or else each element.
    fn analysis_choices(s: &Structure) -> Vec<(String, String)> {
        let mut out = Vec::new();
        if s.has_polymer() {
            out.push(("Cα".to_string(), "protein and name CA".to_string()));
            out.push(("Backbone".to_string(), "backbone".to_string()));
        }
        let by_type = !s.kinds.is_empty()
            && !s
                .kinds
                .iter()
                .all(|k| crate::protein::Element::from_symbol(k).is_known());
        if by_type {
            for k in &s.kinds {
                out.push((format!("Type {k}"), format!("type {k}")));
            }
        } else {
            for (e, _) in s.composition() {
                if e.is_known() {
                    out.push((e.symbol().to_string(), format!("elem {}", e.symbol())));
                }
            }
        }
        out
    }

    /// Which pills are on: those chosen, or by default a protein's alpha
    /// carbons and otherwise every atom.
    fn analysis_chosen(&self, s: &Structure) -> Vec<bool> {
        let n = Self::analysis_choices(s).len();
        match &self.analysis_pills {
            Some(p) if p.len() == n => p.clone(),
            _ if s.has_polymer() => (0..n).map(|i| i == 0).collect(),
            _ => vec![true; n],
        }
    }

    /// The selection the analysis is taken of now.
    fn analysis_selection(&self, s: &Structure) -> String {
        if let Some(c) = &self.analysis_custom {
            return c.clone();
        }
        let parts: Vec<String> = Self::analysis_choices(s)
            .into_iter()
            .zip(self.analysis_chosen(s))
            .filter(|(_, on)| *on)
            .map(|((_, sel), _)| format!("({sel})"))
            .collect();
        if parts.is_empty() {
            "none".into()
        } else {
            parts.join(" or ")
        }
    }

    fn analysis_atoms(&self, s: &Structure) -> Arc<Vec<usize>> {
        let selection = self.analysis_selection(s);
        // The structure and its frames, by where they live, and the words.
        let key = (
            s as *const Structure as usize,
            Arc::as_ptr(&s.frames) as usize,
            s.atoms.len(),
            selection,
        );
        if let Some((k, atoms)) = self.analysis_atoms_found.borrow().as_ref() {
            if *k == key {
                return atoms.clone();
            }
        }
        let atoms: Arc<Vec<usize>> = Arc::new(
            Selection::parse(&key.3)
                .map(|sel| {
                    sel.atoms(s)
                        .iter()
                        .enumerate()
                        .filter(|(_, t)| **t)
                        .map(|(i, _)| i)
                        .collect()
                })
                .unwrap_or_default(),
        );
        *self.analysis_atoms_found.borrow_mut() = Some((key, atoms.clone()));
        atoms
    }

    fn analysis_key(
        &self,
        s: &Structure,
        m: crate::analysis::Measure,
    ) -> (usize, crate::analysis::Measure, String) {
        (
            Arc::as_ptr(&s.frames) as usize,
            m,
            self.analysis_selection(s),
        )
    }

    /// Show measure `m`, working it out in the background unless it has
    /// been for these frames and atoms already. Whatever else was being
    /// worked out is given up.
    pub(super) fn start_measure(&mut self, m: crate::analysis::Measure, cx: &mut Context<Self>) {
        self.analysis_pick = Some(m);
        self.analysis_hover = None;
        let Some(s) = self.structure.clone() else {
            return;
        };
        if s.frame_count() < 2 {
            return;
        }
        let key = self.analysis_key(&s, m);
        if self
            .analysis_task
            .as_ref()
            .is_some_and(|(k, _, _)| *k == key)
        {
            return;
        }
        if let Some((_, cancel, _)) = self.analysis_task.take() {
            cancel.store(true, std::sync::atomic::Ordering::Relaxed);
        }
        self.analysis_done = None;
        let atoms = self.analysis_atoms(&s);
        if self.analysis.contains_key(&key) || !m.applies(&s, &atoms) {
            cx.notify();
            return;
        }
        let done = Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let cancel = Arc::new(std::sync::atomic::AtomicBool::new(false));
        self.analysis_done = Some((done.clone(), crate::analysis::total(&s, m)));
        let finished = Arc::new(std::sync::atomic::AtomicBool::new(false));
        let work = cx.background_executor().spawn({
            let (finished, cancel) = (finished.clone(), cancel.clone());
            async move {
                let plot = crate::analysis::measure(&s, m, &atoms, &cancel, &done);
                finished.store(true, std::sync::atomic::Ordering::Relaxed);
                plot
            }
        });
        let task = cx.spawn({
            let (key, cancel) = (key.clone(), cancel.clone());
            async move |this, cx| {
                while !finished.load(std::sync::atomic::Ordering::Relaxed) {
                    cx.background_executor()
                        .timer(Duration::from_millis(200))
                        .await;
                    if this.update(cx, |_, cx| cx.notify()).is_err() {
                        return;
                    }
                }
                let plot = work.await;
                let _ = this.update(cx, |this, cx| {
                    if this
                        .analysis_task
                        .as_ref()
                        .is_some_and(|(_, c, _)| Arc::ptr_eq(c, &cancel))
                    {
                        this.analysis_task = None;
                        this.analysis_done = None;
                    }
                    if let Some(plot) = plot {
                        this.analysis.insert(key, Arc::new(plot));
                    }
                    cx.notify();
                });
            }
        });
        self.analysis_task = Some((key, cancel, task));
        cx.notify();
    }

    /// Turn an atom pill on or off, and measure again what is shown.
    fn toggle_analysis_atoms(&mut self, i: usize, cx: &mut Context<Self>) {
        let Some(s) = self.structure.clone() else {
            return;
        };
        let mut chosen = self.analysis_chosen(&s);
        if let Some(on) = chosen.get_mut(i) {
            *on = !*on;
        }
        self.analysis_pills = Some(chosen);
        self.analysis_custom = None;
        if let Some(m) = self.analysis_pick {
            self.start_measure(m, cx);
        }
        cx.notify();
    }

    /// Go to the frame a fraction `t` of the way along the plot.
    fn go_to_plotted(&mut self, t: f32, cx: &mut Context<Self>) {
        let n = self.frame_count();
        self.stop_playing();
        self.set_frame(
            (t.clamp(0.0, 1.0) * (n.max(1) - 1) as f32).round() as usize,
            cx,
        );
    }

    /// Words for atom `i` along an RMSF plot: its residue where it stands
    /// for one, or its name in its residue, or its number.
    fn atom_words(s: &Structure, i: usize) -> String {
        let r = s.atoms[i].residue as usize;
        if s.has_polymer() && !s.residues[r].outside_polymer() {
            if s.residues[r].ca == Some(i as u32) {
                Self::residue_label(s, r)
            } else {
                format!("{} {}", s.atoms[i].name(), Self::residue_label(s, r))
            }
        } else {
            format!("atom {}", i + 1)
        }
    }

    /// The trajectory's plots, each worked out only when picked: a value
    /// through time, how much each chosen atom moves, or helices and
    /// strands frame by frame, of the atoms picked above them. The least
    /// and greatest values are marked up the side and the span along the
    /// foot; a click along time goes to that frame.
    fn trajectory_analysis(
        &self,
        structure: &Structure,
        cx: &mut Context<Self>,
    ) -> Option<gpui::AnyElement> {
        use crate::analysis::{Measure, Plot};
        if structure.frame_count() < 2 {
            return None;
        }
        let row = |caption: String, open: bool, cx: &mut Context<Self>| {
            div()
                .flex()
                .items_center()
                .justify_between()
                .gap(px(10.))
                .py(px(8.))
                .border_b_1()
                .border_color(theme::border())
                .child(key("trajectory analysis"))
                .child(
                    div()
                        .flex()
                        .items_center()
                        .gap(px(8.))
                        .child(
                            div()
                                .text_size(px(9.5))
                                .text_color(theme::faint())
                                .child(caption),
                        )
                        .child(
                            div()
                                .id("open-analysis")
                                .debug_selector(|| "open-analysis".into())
                                .px(px(7.))
                                .py(px(2.))
                                .rounded(px(4.))
                                .text_size(px(10.))
                                .cursor_pointer()
                                .text_color(theme::muted())
                                .hover(|s| s.bg(theme::sunken()).text_color(theme::text()))
                                .child(if open { "Hide" } else { "Show" })
                                .on_click(cx.listener(|this, _, _, cx| {
                                    this.analysis_open = !this.analysis_open;
                                    cx.notify();
                                })),
                        ),
                )
        };
        if !self.analysis_open {
            return Some(row(String::new(), false, cx).into_any_element());
        }
        let frames = structure.frame_count();
        let atoms = self.analysis_atoms(structure);
        let small = |label: String| {
            div()
                .font_family(theme::MONO)
                .text_size(px(9.))
                .text_color(theme::faint())
                .child(label)
        };

        // Which atoms: pills, or the selection typed in the console.
        let choices = Self::analysis_choices(structure);
        let chosen = self.analysis_chosen(structure);
        let custom = self.analysis_custom.clone();
        // Each row's name over it, so the pills have the panel's width.
        let atom_pills = div()
            .flex()
            .flex_wrap()
            .items_center()
            .gap(px(4.))
            .children(choices.iter().enumerate().map(|(i, (label, _))| {
                pill_keyed(
                    ("analysis-atoms", i),
                    label.clone(),
                    custom.is_none() && chosen[i],
                )
                .cursor_pointer()
                .on_click(cx.listener(move |this, _, _, cx| this.toggle_analysis_atoms(i, cx)))
            }))
            .children(custom.map(|text| {
                pill("analysis-custom", text, true)
                    .cursor_pointer()
                    .tooltip(|_, cx| {
                        cx.new(
                            |_| Tip("Typed in the console; click to go back to the pills".into()),
                        )
                        .into()
                    })
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.analysis_custom = None;
                        if let Some(m) = this.analysis_pick {
                            this.start_measure(m, cx);
                        }
                    }))
            }))
            .child(small(format!("{} atoms", thousands(atoms.len()))));
        let atom_pills = div()
            .flex()
            .flex_col()
            .gap(px(5.))
            .child(key("atoms"))
            .child(atom_pills);

        // What to measure: each worked out only when picked.
        let pick = self.analysis_pick;
        let measure_pills =
            div()
                .flex()
                .items_center()
                .gap(px(2.))
                .children(Measure::ALL.iter().map(|&m| {
                    let applies = m.applies(structure, &atoms);
                    pill_keyed(("analysis-measure", m as usize), m.short(), pick == Some(m))
                        .flex_none()
                        .tooltip(move |_, cx| cx.new(|_| Tip(m.about().into())).into())
                        .when(!applies, |d| d.text_color(theme::faint()))
                        .when(applies, |d| {
                            d.cursor_pointer().on_click(
                                cx.listener(move |this, _, _, cx| this.start_measure(m, cx)),
                            )
                        })
                }));
        let measure_pills = div()
            .flex()
            .flex_col()
            .gap(px(5.))
            .child(key("measure"))
            .child(measure_pills);

        let plot = pick.and_then(|m| self.analysis.get(&self.analysis_key(structure, m)).cloned());
        let mut caption = String::new();
        let body: gpui::AnyElement = match (pick, plot) {
            (None, _) => small("Pick a measure; each is worked out only when picked.".into())
                .into_any_element(),
            (Some(m), None) => {
                let note = match &self.analysis_done {
                    Some((d, n)) => format!(
                        "Working out {} · {}%",
                        m.label(),
                        (d.load(std::sync::atomic::Ordering::Relaxed) * 100 / (*n).max(1)).min(100)
                    ),
                    None if !m.applies(structure, &atoms) => {
                        format!("{} does not apply to these atoms", m.label())
                    }
                    None => String::new(),
                };
                small(note).into_any_element()
            }
            (Some(m), Some(plot)) => {
                const GUTTER: f32 = 52.0;
                const HEIGHT: f32 = 110.0;
                let width = PLOT - GUTTER;
                let dark = theme::dark();
                let image = {
                    let key = (Arc::as_ptr(&plot) as usize, dark);
                    let mut cache = self.analysis_image.borrow_mut();
                    match cache.as_ref().filter(|(k, _)| *k == key) {
                        Some((_, image)) => image.clone(),
                        None => {
                            let image = analysis_picture(
                                &plot,
                                frames,
                                width as usize * 2,
                                HEIGHT as usize * 2,
                                dark,
                            );
                            *cache = Some((key, image.clone()));
                            image
                        }
                    }
                };
                let unit = |v: f32| {
                    if m.unit().is_empty() {
                        figure(v)
                    } else {
                        format!("{} {}", figure(v), m.unit())
                    }
                };
                // Up the side: the greatest and least values, at their rules;
                // along the foot, where the plot starts and ends.
                let ((top, bottom), (left, right)) = match &*plot {
                    Plot::Time(points) => {
                        let (lo, hi) = plot_range(&points.iter().map(|p| p.1).collect::<Vec<_>>());
                        let last = points.last().map_or(frames, |p| p.0 + 1);
                        (
                            (unit(hi), unit(lo)),
                            ("frame 1".to_string(), format!("frame {last}")),
                        )
                    }
                    Plot::Along(values) => {
                        let (lo, hi) = plot_range(values);
                        let words = |k: Option<&usize>| {
                            k.map(|&i| Self::atom_words(structure, i))
                                .unwrap_or_default()
                        };
                        (
                            (unit(hi), unit(lo)),
                            (words(atoms.first()), words(atoms.last())),
                        )
                    }
                    Plot::Secondary {
                        frames: taken,
                        residues,
                    } => {
                        let words = |r: Option<&usize>| {
                            r.map(|&r| Self::residue_label(structure, r))
                                .unwrap_or_default()
                        };
                        let last = taken.last().map_or(frames, |t| t.0 + 1);
                        (
                            (words(residues.first()), words(residues.last())),
                            ("frame 1".to_string(), format!("frame {last}")),
                        )
                    }
                };
                let secondary = matches!(&*plot, Plot::Secondary { .. });
                let timeline = !matches!(&*plot, Plot::Along(_));
                let place = |t: f32| -> Option<String> {
                    match &*plot {
                        Plot::Time(points) => {
                            let frame = (t * (frames - 1) as f32).round() as usize;
                            let (f, v) = points.iter().min_by_key(|(f, _)| f.abs_diff(frame))?;
                            Some(format!("frame {} · {}", f + 1, unit(*v)))
                        }
                        Plot::Along(values) => {
                            let k = ((t * values.len() as f32) as usize)
                                .min(values.len().checked_sub(1)?);
                            Some(format!(
                                "{} · {}",
                                Self::atom_words(structure, atoms[k]),
                                unit(values[k])
                            ))
                        }
                        Plot::Secondary { frames: taken, .. } => {
                            let (f, _) = taken.get(
                                ((t * taken.len() as f32) as usize)
                                    .min(taken.len().checked_sub(1)?),
                            )?;
                            Some(format!("frame {}", f + 1))
                        }
                    }
                };
                caption = self
                    .analysis_hover
                    .and_then(place)
                    .unwrap_or_else(|| m.label().to_string());
                let at = structure.frame as f32 / (frames - 1).max(1) as f32;
                let store = self.analysis_bounds.clone();
                let side = |label: String, y: f32| {
                    div()
                        .absolute()
                        .right(px(6.))
                        .top(px(y - 6.))
                        .whitespace_nowrap()
                        .font_family(theme::MONO)
                        .text_size(px(9.))
                        .text_color(theme::muted())
                        .child(label)
                };
                let (top_y, bottom_y) = if secondary {
                    (6.0, HEIGHT - 6.0)
                } else {
                    (HEIGHT * PLOT_PAD, HEIGHT * (1.0 - PLOT_PAD))
                };
                let canvas_plot = div()
                    .id("analysis-plot")
                    .debug_selector(|| "analysis-plot".into())
                    .relative()
                    .w(px(width))
                    .h(px(HEIGHT))
                    .rounded(px(4.))
                    .overflow_hidden()
                    .border_1()
                    .border_color(theme::border())
                    .child(
                        gpui::img(gpui::ImageSource::Render(image))
                            .absolute()
                            .inset_0()
                            .size_full(),
                    )
                    .child(
                        canvas(move |bounds, _, _| store.set(bounds), |_, _, _, _| {})
                            .absolute()
                            .inset_0()
                            .size_full(),
                    )
                    .when(timeline, |d| {
                        d.child(
                            div()
                                .absolute()
                                .top_0()
                                .bottom_0()
                                .left(px((width - 2.0) * at))
                                .w(px(1.5))
                                .bg(gpui::Hsla {
                                    a: 0.6,
                                    ..theme::ink()
                                }),
                        )
                    })
                    .on_mouse_move(cx.listener(|this, e: &MouseMoveEvent, _, cx| {
                        let b = this.analysis_bounds.get();
                        let w = f32::from(b.size.width).max(1.0);
                        this.analysis_hover =
                            Some((f32::from(e.position.x - b.origin.x) / w).clamp(0.0, 1.0));
                        cx.notify();
                    }))
                    .on_hover(cx.listener(|this, hovered: &bool, _, cx| {
                        if !*hovered {
                            this.analysis_hover = None;
                            cx.notify();
                        }
                    }))
                    .on_click(cx.listener(move |this, _, _, cx| {
                        if let (Some(t), true) = (this.analysis_hover, timeline) {
                            this.go_to_plotted(t, cx);
                        }
                    }));
                div()
                    .flex()
                    .flex_col()
                    .gap(px(3.))
                    .child(
                        div()
                            .flex()
                            .child(
                                div()
                                    .relative()
                                    .w(px(GUTTER))
                                    .h(px(HEIGHT))
                                    .flex_none()
                                    .child(side(top, top_y))
                                    .child(side(bottom, bottom_y)),
                            )
                            .child(canvas_plot),
                    )
                    .child(
                        div()
                            .flex()
                            .justify_between()
                            .pl(px(GUTTER))
                            .font_family(theme::MONO)
                            .text_size(px(9.))
                            .text_color(theme::muted())
                            .child(left)
                            .child(right),
                    )
                    .into_any_element()
            }
        };
        let fit = div()
            .id("fit-frames")
            .debug_selector(|| "fit-frames".into())
            .px(px(7.))
            .py(px(2.))
            .rounded(px(4.))
            .text_size(px(10.))
            .cursor_pointer()
            .when(self.boundary.fitted, |d| {
                d.bg(theme::sunken()).text_color(theme::ink())
            })
            .when(!self.boundary.fitted, |d| {
                d.text_color(theme::muted())
                    .hover(|s| s.bg(theme::sunken()).text_color(theme::text()))
            })
            .child("Fit frames onto the first")
            .on_click(cx.listener(|this, _, _, cx| {
                if let Err(e) = this.periodic(Periodic::Fit, cx) {
                    this.error = Some(format!("{e}").into());
                }
            }));
        Some(
            div()
                .flex()
                .flex_col()
                .child(row(caption, true, cx))
                .child(
                    div()
                        .pt(px(10.))
                        .flex()
                        .flex_col()
                        .gap(px(8.))
                        .child(atom_pills)
                        .child(measure_pills)
                        .child(body)
                        .child(div().flex().child(fit)),
                )
                .into_any_element(),
        )
    }

    /// The volume's surface as last built, if it is shown.
    pub(crate) fn iso_drawn(&self) -> Option<Arc<Scene>> {
        let volume = self.structure.as_ref()?.volume.as_ref()?;
        let (iso, key, scene) = self.iso_cache.as_ref()?;
        (self.iso.shown && *iso == self.iso && *key == Arc::as_ptr(volume) as usize)
            .then(|| scene.clone())
    }

    /// The volume's surface: shown or not, solid or a mesh, and its level.
    pub(crate) fn set_iso(&mut self, iso: Iso, cx: &mut Context<Self>) {
        self.iso = iso;
        self.rebuild(cx);
    }

    /// Controls for the volume, under the structure's name.
    fn volume_tools(
        &self,
        structure: &Structure,
        cx: &mut Context<Self>,
    ) -> Option<gpui::AnyElement> {
        let volume = structure.volume.as_ref()?;
        let iso = self.iso;
        let chip = |id: &'static str, label: String, on: bool, next: Iso| {
            div()
                .id(id)
                .px(px(7.))
                .py(px(2.))
                .rounded(px(4.))
                .text_size(px(10.))
                .cursor_pointer()
                .when(on, |d| d.bg(theme::sunken()).text_color(theme::text()))
                .when(!on, |d| {
                    d.text_color(theme::muted())
                        .hover(|s| s.bg(theme::sunken()).text_color(theme::text()))
                })
                .child(label)
                .on_click(cx.listener(move |this, _, _, cx| this.set_iso(next, cx)))
        };
        let [a, b, c] = volume.dims;
        let lobes = if volume.signed() { " · ± lobes" } else { "" };
        Some(
            div()
                .pt(px(14.))
                .flex()
                .flex_col()
                .gap(px(6.))
                .child(key(&format!("volume · {a} × {b} × {c}{lobes}")))
                .child(
                    div()
                        .flex()
                        .flex_wrap()
                        .items_center()
                        .gap(px(2.))
                        .child(chip(
                            "iso-solid",
                            "Surface".into(),
                            iso.shown && !iso.mesh,
                            Iso {
                                shown: true,
                                mesh: false,
                                ..iso
                            },
                        ))
                        .child(chip(
                            "iso-mesh",
                            "Mesh".into(),
                            iso.shown && iso.mesh,
                            Iso {
                                shown: true,
                                mesh: true,
                                ..iso
                            },
                        ))
                        .child(chip(
                            "iso-hide",
                            "Hide".into(),
                            !iso.shown,
                            Iso {
                                shown: false,
                                ..iso
                            },
                        ))
                        .child(div().w(px(8.)))
                        .child(chip(
                            "iso-down",
                            "−".into(),
                            false,
                            Iso {
                                level: iso.level / 1.25,
                                ..iso
                            },
                        ))
                        .child(
                            div()
                                .px(px(4.))
                                .font_family(theme::MONO)
                                .text_size(px(10.))
                                .text_color(theme::text())
                                .child(format!("level {}", figure(iso.level))),
                        )
                        .child(chip(
                            "iso-up",
                            "+".into(),
                            false,
                            Iso {
                                level: iso.level * 1.25,
                                ..iso
                            },
                        )),
                )
                .into_any_element(),
        )
    }

    /// Buttons under the cell's size for what can be done across its walls.
    fn periodic_tools(
        &self,
        structure: &Structure,
        cx: &mut Context<Self>,
    ) -> Option<gpui::AnyElement> {
        if !periodic(structure) {
            return None;
        }
        let b = self.boundary;
        let chip = |id: &'static str, label: &'static str, tip: &'static str, op: Periodic| {
            // Lit while it is in effect; a click on a lit one turns it off.
            let on = match op {
                Periodic::Wrap => b.place == Place::Wrap,
                Periodic::Whole => b.place == Place::Whole,
                Periodic::InTime => b.place == Place::InTime,
                Periodic::Replicate(n) => b.replicate == Some(n),
                Periodic::Fit => b.fitted,
                Periodic::Original => false,
            };
            div()
                .id(id)
                .debug_selector(move || id.into())
                .flex_none()
                .whitespace_nowrap()
                .tooltip(move |_, cx| cx.new(|_| Tip(tip.into())).into())
                .px(px(7.))
                .py(px(2.))
                .rounded(px(4.))
                .text_size(px(10.))
                .cursor_pointer()
                .when(on, |d| d.bg(theme::sunken()).text_color(theme::ink()))
                .when(!on, |d| {
                    d.text_color(theme::muted())
                        .hover(|s| s.bg(theme::sunken()).text_color(theme::text()))
                })
                .child(label)
                .on_click(cx.listener(move |this, _, _, cx| {
                    if let Err(e) = this.periodic(op, cx) {
                        this.error = Some(format!("{e}").into());
                    }
                    cx.notify();
                }))
        };
        Some(
            // One line, never wrapped: short names, each said in full on hover.
            div()
                .flex()
                .items_center()
                .gap(px(2.))
                .py(px(6.))
                .border_b_1()
                .border_color(theme::border())
                .child(chip("pbc-wrap", "Wrap", "Put every atom back inside the cell", Periodic::Wrap))
                .child(chip("pbc-whole", "Whole", "Mend molecules cut by the cell's walls", Periodic::Whole))
                .when(structure.frame_count() > 1, |d| {
                    d.child(chip(
                        "pbc-time",
                        "Unwrap",
                        "Follow atoms through the walls from frame to frame, rather than letting them jump back",
                        Periodic::InTime,
                    ))
                })
                .child(chip(
                    "pbc-replicate",
                    "2×2×2",
                    "Repeat the cell twice along each edge; replicate a b c in the console for other sizes",
                    Periodic::Replicate([2, 2, 2]),
                ))
                .when(b != Boundary::default(), |d| {
                    d.child(chip("pbc-original", "Original", "Undo these, back to the structure as opened", Periodic::Original))
                })
                .into_any_element(),
        )
    }

    /// Colour by `c`; by property, the one chosen before or else the first.
    pub(super) fn set_coloring(&mut self, c: Coloring, cx: &mut Context<Self>) {
        if c == Coloring::Property {
            if let Some(s) = self.structure.as_mut() {
                if s.prop_color.is_none() && !s.props.is_empty() {
                    Arc::make_mut(s).prop_color = Some(crate::protein::PropertyColor {
                        index: 0,
                        range: None,
                    });
                }
            }
        }
        self.coloring = c;
        self.rebuild(cx);
    }

    /// Colour by property `index`, over `range` or its own span.
    pub(super) fn colour_by_property(
        &mut self,
        index: usize,
        range: Option<(f32, f32)>,
        cx: &mut Context<Self>,
    ) {
        let Some(s) = self.structure.as_mut() else {
            return;
        };
        Arc::make_mut(s).prop_color = Some(crate::protein::PropertyColor { index, range });
        self.coloring = Coloring::Property;
        self.rebuild(cx);
    }

    fn set_ligand_paint(
        &mut self,
        name: String,
        paint: Option<LigandPaint>,
        cx: &mut Context<Self>,
    ) {
        let Some(s) = self.structure.as_mut() else {
            return;
        };
        let s = Arc::make_mut(s);
        match paint {
            Some(p) => s.ligand_paint.insert(name, p),
            None => s.ligand_paint.remove(&name),
        };
        self.rebuild(cx);
    }

    /// Give chain `i` its own colour, or `None` for the palette's. A colour
    /// shows when colouring by chain, so that is switched to.
    fn set_chain_color(&mut self, i: usize, color: Option<[f32; 3]>, cx: &mut Context<Self>) {
        let Some(s) = self.structure.as_mut() else {
            return;
        };
        if let Some(chain) = Arc::make_mut(s).chains.get_mut(i) {
            chain.color = color;
        }
        if self.compare.is_none() {
            self.coloring = Coloring::Chain;
        }
        self.rebuild(cx);
    }

    /// Give every chain back the palette's colour.
    fn reset_chain_colors(&mut self, cx: &mut Context<Self>) {
        if let Some(s) = self.structure.as_mut() {
            for chain in &mut Arc::make_mut(s).chains {
                chain.color = None;
            }
        }
        self.rebuild(cx);
    }

    /// The chains and species in view that keep a colour of their own,
    /// whatever the palette, by name.
    fn own_colors(&self) -> Vec<String> {
        let Some(s) = &self.structure else {
            return Vec::new();
        };
        let chains = s
            .chains
            .iter()
            .enumerate()
            .filter(|(_, c)| c.color.is_some())
            .map(|(i, c)| {
                if c.id.trim().is_empty() {
                    format!("Chain {}", i + 1)
                } else {
                    format!("Chain {}", c.id)
                }
            });
        let elements = s
            .composition()
            .into_iter()
            .filter(|(e, _)| e.is_known())
            .filter(|(e, _)| self.appearance.get(&Species::Element(e.0)).color.is_some())
            .map(|(e, _)| e.symbol().to_string());
        let kinds = s
            .kinds
            .iter()
            .filter(|k| {
                self.appearance
                    .get(&Species::Kind((*k).clone()))
                    .color
                    .is_some()
            })
            .map(|k| format!("Type {k}"));
        let ligands = s
            .ligands()
            .into_iter()
            .filter(|(name, _)| s.ligand_paint.contains_key(name))
            .map(|(name, _)| name);
        chains.chain(ligands).chain(elements).chain(kinds).collect()
    }

    /// Give every chain and species back the palette's colour.
    fn reset_colors(&mut self, cx: &mut Context<Self>) {
        self.kept = None;
        self.appearance.clear_colors();
        if let Some(s) = self.structure.as_mut() {
            Arc::make_mut(s).ligand_paint.clear();
        }
        self.reset_chain_colors(cx);
    }

    /// Give a species its own colour, or `None` for its element's. A colour
    /// shows when colouring by element, or on the ligands beside a polymer;
    /// without ligands, colouring by element is switched to.
    fn set_color(&mut self, species: Species, color: Option<[f32; 3]>, cx: &mut Context<Self>) {
        let mut look = self.appearance.get(&species);
        look.color = color;
        self.appearance.set(species, look);
        let ligands = self.structure.as_ref().is_some_and(|s| s.has_ligands());
        if self.compare.is_none() && !ligands {
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
        let camera = self
            .camera
            .resized((width, height), (width * scale, height * scale));
        let (width, height) = (width * scale, height * scale);
        let scene = self.scene.clone();
        let style = self.style;
        let transparent = self.transparent;
        // The labels where their atoms are, set at the picture's scale.
        let pinned: Vec<([f32; 3], String)> = self
            .shown(0)
            .map(|s| {
                self.labels
                    .iter()
                    .filter_map(|l| Some((s.atoms.get(l.atom as usize)?.pos, l.text.clone())))
                    .collect()
            })
            .unwrap_or_default();
        let label_size =
            10.5 * (width as f32 / f32::from(self.cache.borrow().bounds.size.width).max(1.0));
        let ink = theme::view().map(|c| if c > 0.5 { 0.07 } else { 0.92 });
        let page = theme::view();
        let traced = self.traced.then(|| {
            let done = Arc::new(std::sync::atomic::AtomicUsize::new(0));
            (done, crate::trace::rows(height, width))
        });
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
            let finished = Arc::new(std::sync::atomic::AtomicBool::new(false));
            let _ = this.update(cx, |this, cx| {
                this.saving = traced.clone();
                cx.notify();
            });
            let following = traced.is_some();
            let drawn = cx.background_executor().spawn({
                let finished = finished.clone();
                async move {
                    let rgba = if let Some((done, _)) = &traced {
                        let control = crate::trace::Control {
                            cancel: None,
                            done: Some(done),
                        };
                        let t =
                            crate::trace::trace(&scene, &camera, width, height, style, &control)
                                .expect("not cancelled");
                        if transparent {
                            t.transparent()
                        } else {
                            let mut rgba = t.develop(theme::view()).bgra;
                            for p in rgba.chunks_exact_mut(4) {
                                p.swap(0, 2);
                            }
                            rgba
                        }
                    } else if transparent {
                        render::render_transparent(&scene, &camera, width, height, style)
                    } else {
                        let frame =
                            render::render(&scene, &camera, width, height, theme::view(), style);
                        // The rasterizer hands back BGRA; PNG wants it the
                        // other way round.
                        let mut rgba = frame.bgra;
                        for p in rgba.chunks_exact_mut(4) {
                            p.swap(0, 2);
                        }
                        rgba
                    };
                    let mut rgba = rgba;
                    if !pinned.is_empty() {
                        if let Some(font) = crate::movie::font() {
                            let mut pic = crate::movie::Picture {
                                width,
                                height,
                                rgba,
                            };
                            let labels: Vec<([f32; 3], &str)> =
                                pinned.iter().map(|(p, t)| (*p, t.as_str())).collect();
                            crate::movie::draw_labels(
                                &mut pic, &font, &camera, &labels, label_size, ink, page,
                            );
                            rgba = pic.rgba;
                        }
                    }
                    let saved = image::save_buffer(
                        &path,
                        &rgba,
                        width as u32,
                        height as u32,
                        image::ColorType::Rgba8,
                    )
                    .map(|()| path);
                    finished.store(true, std::sync::atomic::Ordering::Relaxed);
                    saved
                }
            });
            // Shown again now and then while it is traced, so the progress moves.
            while following && !finished.load(std::sync::atomic::Ordering::Relaxed) {
                cx.background_executor()
                    .timer(std::time::Duration::from_millis(200))
                    .await;
                let _ = this.update(cx, |_, cx| cx.notify());
            }
            let drawn = drawn.await;
            let _ = this.update(cx, |this, cx| {
                this.saving = None;
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

    /// Ask where to write a file, suggesting `name` on the desktop, and write
    /// what `contents` makes there.
    fn save_text(&mut self, name: String, contents: String, cx: &mut Context<Self>) {
        self.open = None;
        let directory = dirs::home_dir()
            .map(|home| home.join("Desktop"))
            .unwrap_or_else(|| std::path::PathBuf::from("."));
        let chosen = cx.prompt_for_new_path(&directory, Some(&name));
        self.export = Some(cx.spawn(async move |this, cx| {
            let Ok(Ok(Some(path))) = chosen.await else {
                return;
            };
            let written = std::fs::write(&path, contents).map(|()| path);
            let _ = this.update(cx, |this, cx| {
                match written {
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

    /// What can be saved: the structure as PDB or mmCIF, the picked
    /// residue's chain, the binding site looked into, and the session.
    /// Write out the whole structure, as PDB or as mmCIF.
    pub fn save_structure(&mut self, cif: bool, cx: &mut Context<Self>) {
        let Some(s) = self.structure.clone() else {
            return;
        };
        let stem = s.id.to_lowercase();
        let (text, ext) = if cif {
            (write::mmcif(&s, |_| true), "cif")
        } else {
            (write::pdb(&s, |_| true), "pdb")
        };
        self.save_text(format!("{stem}.{ext}"), text, cx);
    }

    /// Write out the session: what is open and how it is shown.
    pub fn save_session(&mut self, cx: &mut Context<Self>) {
        if let Some(json) = self.session_json() {
            let stem = self
                .structure
                .as_ref()
                .map(|s| s.id.to_lowercase())
                .unwrap_or_default();
            self.save_text(format!("{stem}.{}", session::EXTENSION), json, cx);
        }
    }

    fn save_menu(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let Some(s) = self.structure.clone() else {
            return menu(240.);
        };
        let stem = s.id.to_lowercase();
        let chain = self
            .selected
            .filter(|t| t.which == 0)
            .map(|t| s.residues[s.atoms[t.picked.atom].residue as usize].chain as usize);
        let site = self.site.clone();
        let whole = |id: &'static str, label: &'static str, cif: bool| {
            item(id, label.into(), false)
                .debug_selector(move || id.into())
                .on_click(cx.listener(move |this, _, _, cx| this.save_structure(cif, cx)))
        };
        menu(240.)
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.open = None;
                cx.notify();
            }))
            .child(whole("save-pdb", "Structure as PDB…", false))
            .child(whole("save-cif", "Structure as mmCIF…", true))
            .when_some(chain, |d, c| {
                let s = s.clone();
                let id = s.chains[c].id.clone();
                let stem = format!("{stem}_{id}");
                d.child(
                    item("save-chain", format!("Chain {id} as PDB…"), false)
                        .debug_selector(|| "save-chain".into())
                        .on_click(cx.listener(move |this, _, _, cx| {
                            let text = write::pdb(&s, |i| {
                                s.residues[s.atoms[i].residue as usize].chain as usize == c
                            });
                            this.save_text(format!("{stem}.pdb"), text, cx);
                        })),
                )
            })
            .when_some(site, |d, site| {
                let s = s.clone();
                let name = s.residues[site.ligand].name.to_lowercase();
                let stem = format!("{stem}_{name}_site");
                d.child(
                    item("save-site", "Binding site as PDB…".into(), false)
                        .debug_selector(|| "save-site".into())
                        .on_click(cx.listener(move |this, _, _, cx| {
                            let text = write::pdb(&s, |i| {
                                let r = s.atoms[i].residue as usize;
                                r == site.ligand || site.residues.contains(&r)
                            });
                            this.save_text(format!("{stem}.pdb"), text, cx);
                        })),
                )
            })
            .child(divide())
            .child(
                item("save-session", "Session…".into(), false)
                    .debug_selector(|| "save-session".into())
                    .on_click(cx.listener(|this, _, _, cx| this.save_session(cx))),
            )
    }

    /// After a change of palette, say for a few seconds what kept its own
    /// colour through it, so the colours that did not move are no surprise.
    fn say_kept(&mut self, cx: &mut Context<Self>) {
        self.kept = kept_note(&self.own_colors());
        if self.kept.is_none() {
            self.kept_task = None;
            return;
        }
        self.kept_task = Some(cx.spawn(async move |this, cx| {
            cx.background_executor().timer(Duration::from_secs(5)).await;
            let _ = this.update(cx, |this, cx| {
                this.kept = None;
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
    /// Cut the structure down to a slab a fifth of it thick, through the
    /// middle, or show it whole again.
    pub fn toggle_perspective(&mut self, cx: &mut Context<Self>) {
        self.camera.perspective = !self.camera.perspective;
        cx.notify();
    }

    pub fn toggle_slab(&mut self, cx: &mut Context<Self>) {
        let r = self.camera.radius;
        self.camera.slab = match self.camera.slab {
            Some(_) => None,
            None => Some([0.2 * r, -0.2 * r]),
        };
        cx.notify();
    }

    /// Open the console, or close it, putting the cursor in it.
    pub fn toggle_console(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        self.console_open = !self.console_open;
        if self.console_open {
            window.focus(&self.console_input.read(cx).focus_handle(cx));
        } else {
            window.focus(&self.focus);
        }
        cx.notify();
    }

    /// Lines to run once what is being opened has opened, or at once if
    /// nothing is, as `--run` gives them.
    pub fn queue_script(&mut self, lines: Vec<String>, cx: &mut Context<Self>) {
        if self.busy.is_some() {
            self.script = lines;
        } else {
            self.run_script(lines, cx);
        }
    }

    /// The console, at the foot of the viewport: the last few commands and
    /// what came of them, over the line to type the next in.
    fn console_view(&self) -> Option<gpui::AnyElement> {
        if !self.console_open {
            return None;
        }
        let bottom = if self.frame_count() > 1 { 104. } else { 62. };
        let start = self.console_log.len().saturating_sub(5);
        Some(
            div()
                .id("console")
                .debug_selector(|| "console".into())
                .absolute()
                .left(px(20.))
                .right(px(20.))
                .bottom(px(bottom))
                .flex()
                .flex_col()
                .gap(px(2.))
                .p(px(8.))
                .rounded(px(8.))
                .bg(theme::surface())
                .border_1()
                .border_color(theme::border())
                .shadow_sm()
                .on_mouse_down(MouseButton::Left, |_, _, cx| cx.stop_propagation())
                .children(self.console_log[start..].iter().map(|(ok, text)| {
                    div()
                        .font_family(theme::MONO)
                        .text_size(px(10.5))
                        .text_color(if *ok { theme::muted() } else { theme::ink() })
                        .child(clip(text, 140))
                }))
                .child(
                    div()
                        .flex()
                        .items_center()
                        .gap(px(6.))
                        .font_family(theme::MONO)
                        .text_size(px(11.5))
                        .text_color(theme::text())
                        .child(div().text_color(theme::faint()).child("›"))
                        .child(div().flex_1().child(self.console_input.clone())),
                )
                .into_any_element(),
        )
    }

    /// Show the keyboard shortcuts, or put them away.
    pub fn toggle_shortcuts(&mut self, cx: &mut Context<Self>) {
        self.shortcuts = !self.shortcuts;
        cx.notify();
    }

    /// The slab's place in the structure: a gauge at the right of the
    /// viewport, the structure's depth along it and the slab a band on it,
    /// with how to move it.
    fn slab_gauge(&self) -> Option<gpui::AnyElement> {
        let [near, far] = self.camera.slab?;
        let r = self.camera.radius.max(1e-3);
        const H: f32 = 160.0;
        // Near at the top.
        let at = |z: f32| ((r - z) / (2.0 * r)).clamp(0.0, 1.0) * H;
        let (top, bottom) = (at(near), at(far));
        let hint = format!(
            "{}-scroll: slab depth · {}: slab thickness",
            theme::keys("⌥"),
            theme::keys("⌥⇧")
        );
        Some(
            div()
                .absolute()
                .right(px(14.))
                .top(px(120.))
                .flex()
                .flex_col()
                .items_end()
                .gap(px(6.))
                .child(
                    div()
                        .text_size(px(9.))
                        .text_color(theme::faint())
                        .child(spaced("SLAB")),
                )
                .child(
                    div()
                        .relative()
                        .w(px(6.))
                        .h(px(H))
                        .rounded(px(3.))
                        .bg(theme::sunken())
                        .child(
                            div()
                                .absolute()
                                .left_0()
                                .right_0()
                                .top(px(top))
                                .h(px((bottom - top).max(2.0)))
                                .rounded(px(3.))
                                .bg(theme::ink()),
                        ),
                )
                .child(
                    div()
                        .w(px(90.))
                        .text_size(px(9.))
                        .text_color(theme::faint())
                        .child(hint),
                )
                .into_any_element(),
        )
    }

    /// Every shortcut, keys and mouse, over the viewport until dismissed.
    fn shortcut_sheet(&self, cx: &mut Context<Self>) -> Option<gpui::AnyElement> {
        if !self.shortcuts {
            return None;
        }
        let groups: [(&str, &[(&str, &str)]); 4] = [
            (
                "Tabs and files",
                &[
                    ("⌘O", "Open a file"),
                    ("⌘T", "New tab"),
                    ("⇧⌘T", "New fold tab"),
                    ("⌘W", "Close the tab"),
                    ("⌘1 – ⌘9", "Go to a tab"),
                    ("⌃Tab", "Next tab"),
                    ("⇧⌘S", "Export an image"),
                    ("⇧⌘E", "Export a film"),
                ],
            ),
            (
                "Viewing",
                &[
                    ("⌥⌘1 – ⌥⌘6", "Cartoon … Surface"),
                    ("⌃⌘1 – ⌃⌘4", "Illustrative … Studio"),
                    ("⌘0", "Reset the view"),
                    ("⇧⌘D", "Toggle light or dark theme"),
                    ("⇧⌘L", "Panel beside or under the view"),
                    ("⌘K", "The console"),
                    ("⌘Z  ⇧⌘Z", "Undo, redo"),
                    ("⌘/", "Show keyboard shortcuts"),
                ],
            ),
            (
                "In the viewport",
                &[
                    ("Drag", "Turn"),
                    ("⇧ Drag", "Move"),
                    ("⌘ Drag", "Roll clockwise or back"),
                    ("Scroll", "Zoom"),
                    ("Double-click", "Fit the structure to the view"),
                    ("Click", "Pick an atom or residue"),
                    ("⌘ Click", "Fly to an atom or residue"),
                    ("⌥ Scroll", "Move the slab"),
                    ("⌥⇧ Scroll", "Increase slab thickness"),
                ],
            ),
            (
                "In the field",
                &[
                    ("Return", "Open the best match"),
                    ("↑ ↓", "Choose a result"),
                    ("Esc", "Close the results"),
                ],
            ),
        ];
        // Written the Mac way; elsewhere a chord is spelled out, and a lone
        // modifier before a mouse action named.
        let say = |k: &str| -> String {
            k.split(' ')
                .map(|w| match w {
                    "⌥" if !cfg!(target_os = "macos") => "Alt".to_string(),
                    "⇧" if !cfg!(target_os = "macos") => "Shift".to_string(),
                    "⌥⇧" if !cfg!(target_os = "macos") => "Alt+Shift".to_string(),
                    w if w.starts_with(['⌘', '⌥', '⌃', '⇧']) => theme::keys(w),
                    w => w.to_string(),
                })
                .collect::<Vec<_>>()
                .join(" ")
        };
        let row = |k: &str, what: &str| {
            let k = say(k);
            div()
                .flex()
                .justify_between()
                .gap(px(16.))
                .child(
                    div()
                        .text_size(px(11.))
                        .text_color(theme::text())
                        .child(what.to_string()),
                )
                .child(
                    div()
                        .font_family(theme::MONO)
                        .text_size(px(10.5))
                        .text_color(theme::muted())
                        .child(k),
                )
        };
        Some(
            div()
                .id("shortcuts")
                .debug_selector(|| "shortcuts".into())
                .absolute()
                .inset_0()
                .flex()
                .items_center()
                .justify_center()
                .bg(gpui::Hsla {
                    a: 0.35,
                    ..theme::ink()
                })
                .on_mouse_down(MouseButton::Left, |_, _, cx| cx.stop_propagation())
                .on_click(cx.listener(|this, _, _, cx| this.toggle_shortcuts(cx)))
                .child(
                    div()
                        .w(px(560.))
                        .p(px(20.))
                        .rounded(px(10.))
                        .bg(theme::surface())
                        .border_1()
                        .border_color(theme::border())
                        .shadow_lg()
                        .flex()
                        .flex_wrap()
                        .gap_x(px(28.))
                        .gap_y(px(16.))
                        .children(groups.iter().map(|(name, keys)| {
                            div()
                                .w(px(240.))
                                .flex()
                                .flex_col()
                                .gap(px(5.))
                                .child(div().pb(px(2.)).child(key(name)))
                                .children(keys.iter().map(|(k, w)| row(k, w)))
                        })),
                )
                .into_any_element(),
        )
    }

    /// Set the camera turning by itself, or with the same again stop it.
    pub fn toggle_spin(&mut self, cx: &mut Context<Self>) {
        self.toggle_turn(Turn::Spin, cx);
    }

    pub fn toggle_rock(&mut self, cx: &mut Context<Self>) {
        self.toggle_turn(Turn::Rock, cx);
    }

    fn toggle_turn(&mut self, turn: Turn, cx: &mut Context<Self>) {
        self.open = None;
        if self.turning == Some(turn) {
            self.halt();
            cx.notify();
            return;
        }
        self.halt();
        self.turning = Some(turn);
        let start = std::time::Instant::now();
        self.turn_task = Some(cx.spawn(async move |this, cx| {
            let mut last = 0.0f32;
            loop {
                cx.background_executor()
                    .timer(Duration::from_millis(16))
                    .await;
                let t = start.elapsed().as_secs_f32();
                // Spinning: a turn every 14 seconds. Rocking: 40 degrees each
                // way, there and back in six.
                let step = match turn {
                    Turn::Spin => (t - last) * std::f32::consts::TAU / 14.0,
                    Turn::Rock => {
                        let w = std::f32::consts::TAU / 6.0;
                        0.7 * ((w * t).sin() - (w * last).sin())
                    }
                };
                last = t;
                let going = this.update(cx, |this, cx| {
                    if this.turning != Some(turn) {
                        return false;
                    }
                    this.camera.rotate(step, 0.0);
                    cx.notify();
                    true
                });
                if !matches!(going, Ok(true)) {
                    return;
                }
            }
        }));
        cx.notify();
    }

    fn coast(&mut self, velocity: [f32; 2], cx: &mut Context<Self>) {
        self.spin = Some(velocity);
        self.hover = None;
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
        self.fly_task = None;
        self.turning = None;
        self.turn_task = None;
        self.last_step = [0.0; 2];
    }

    fn on_viewport_down(&mut self, e: &MouseDownEvent, _: &mut Window, cx: &mut Context<Self>) {
        self.halt();
        self.open = None;
        self.hover = None;
        self.press = (e.button == MouseButton::Left && !e.modifiers.shift && e.click_count < 2)
            .then_some(e.position);
        if e.click_count >= 2 {
            let all: Vec<&Structure> = self.visible().into_iter().map(|s| s.as_ref()).collect();
            if !all.is_empty() {
                self.camera = Camera {
                    perspective: self.camera.perspective,
                    slab: self.camera.slab,
                    ..Camera::frame_all(&all)
                };
            }
        } else {
            self.drag = Some(if e.modifiers.shift || e.button == MouseButton::Right {
                Drag::Pan
            } else if e.modifiers.secondary() {
                Drag::Roll
            } else {
                Drag::Orbit
            });
            self.last_mouse = e.position;
        }
        cx.notify();
    }

    fn on_mouse_move(&mut self, e: &MouseMoveEvent, window: &mut Window, cx: &mut Context<Self>) {
        if self.resizing {
            if theme::vertical() {
                let bottom = window.viewport_size().height;
                theme::set_band_height(
                    f32::from(bottom - e.position.y).min(f32::from(bottom) - 240.),
                );
            } else {
                let right = f32::from(window.viewport_size().width);
                theme::set_side_width(f32::from(e.position.x).min(right - 320.));
            }
            cx.notify();
            return;
        }
        if let Some((species, knob)) = self.sliding.clone() {
            self.slide(species, knob, e.position.x, cx);
            return;
        }
        if self.scrubbing {
            self.scrub_to(e.position.x, cx);
            return;
        }
        let Some(drag) = self.drag else {
            self.set_hover(Some(e.position), cx);
            return;
        };
        if let Some(p) = self.press {
            let (dx, dy) = (f32::from(e.position.x - p.x), f32::from(e.position.y - p.y));
            if dx * dx + dy * dy > CLICK_SLOP * CLICK_SLOP {
                self.press = None;
            }
        }
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
            Drag::Roll => {
                // The angle the pointer sweeps about the middle of the view,
                // so the molecule turns with the hand.
                let bounds = self.cache.borrow().bounds;
                let cx_ = f32::from(bounds.origin.x + bounds.size.width / 2.0) + self.camera.pan[0];
                let cy_ =
                    f32::from(bounds.origin.y + bounds.size.height / 2.0) + self.camera.pan[1];
                let (x, y) = (f32::from(e.position.x) - cx_, f32::from(e.position.y) - cy_);
                let (px_, py_) = (x - dx, y - dy);
                if x * x + y * y > 4.0 && px_ * px_ + py_ * py_ > 4.0 {
                    let a = (px_ * y - py_ * x).atan2(px_ * x + py_ * y);
                    self.camera.roll(a);
                }
            }
        }
        cx.notify();
    }

    fn on_mouse_up(&mut self, e: &MouseUpEvent, _: &mut Window, cx: &mut Context<Self>) {
        if std::mem::take(&mut self.resizing) {
            self.persist();
        }
        self.scrubbing = false;
        self.sliding = None;
        let was = self.drag.take();
        if let Some(p) = self.press.take() {
            // A click, not a drag: the view has not moved, so nothing coasts.
            self.last_step = [0.0; 2];
            let target = self.pick_at(p);
            if self.measuring {
                self.measure_click(target, cx);
            } else if e.modifiers.secondary() {
                // ⌘-click: pick it and fly to it. Never a second click that
                // lets go, and a miss leaves the pick as it was.
                if target.is_some() {
                    if target != self.selected {
                        self.selected = target;
                        self.compose(cx);
                    }
                    self.focus_selected(cx);
                }
            } else {
                self.select(target, cx);
            }
            self.set_hover(Some(e.position), cx);
        }
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
        self.hover = None;
        // With a slab, the option key moves it through the structure, and
        // with shift too, thickens or thins it.
        if let (true, Some([near, far])) = (e.modifiers.alt, self.camera.slab) {
            let dy = match e.delta {
                ScrollDelta::Pixels(p) => f32::from(p.y) + f32::from(p.x),
                ScrollDelta::Lines(p) => (p.y + p.x) * 18.0,
            };
            let r = self.camera.radius;
            let step = dy * r * 0.004;
            self.camera.slab = Some(if e.modifiers.shift {
                let mid = (near + far) * 0.5;
                let half = ((near - far) * 0.5 + step * 0.5).clamp(0.25, 2.0 * r);
                [mid + half, mid - half]
            } else {
                let shift = step.clamp(-2.0 * r - far, 2.0 * r - near);
                [near + shift, far + shift]
            });
            cx.notify();
            return;
        }
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

/// Whether the structure's cell is a periodic box worth working across: a
/// simulation's, a trajectory's, or a small molecule's or crystal's, and not
/// the unit cell a protein crystallised in.
fn periodic(s: &Structure) -> bool {
    s.cell.is_some()
        && (!s.has_polymer()
            || s.frame_count() > 1
            || matches!(
                s.format,
                Some("GROMACS" | "LAMMPS data" | "LAMMPS dump" | "XYZ" | "XSF" | "VASP")
            ))
}

/// Text pinned to an atom of the structure shown.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Label {
    pub atom: u32,
    pub text: String,
}

/// How a volume is drawn: its surface at a level, solid or as a mesh.
#[derive(Clone, Copy, PartialEq, Debug)]
pub(crate) struct Iso {
    pub shown: bool,
    pub level: f32,
    pub mesh: bool,
}

impl Default for Iso {
    fn default() -> Self {
        Iso {
            shown: true,
            level: 1.0,
            mesh: false,
        }
    }
}

/// How atoms are placed across the cell's walls: as the file has them, or
/// one of three ways, each undoing the others.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub(crate) enum Place {
    #[default]
    AsRead,
    Wrap,
    Whole,
    InTime,
}

/// The periodic-boundary changes on now.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub(crate) struct Boundary {
    pub place: Place,
    pub fitted: bool,
    pub replicate: Option<[usize; 3]>,
}

/// What can be done to a structure across its periodic boundaries.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Periodic {
    Wrap,
    Whole,
    InTime,
    Replicate([usize; 3]),
    /// Every frame laid over the first.
    Fit,
    Original,
}

/// A value as a key shows it: no more figures than it needs.
pub(crate) fn figure(v: f32) -> String {
    let a = v.abs();
    if a == 0.0 {
        return "0".into();
    }
    if (0.01..1e5).contains(&a) {
        let digits = if a >= 100.0 {
            0
        } else if a >= 1.0 {
            2
        } else {
            3
        };
        format!("{v:.digits$}")
    } else {
        format!("{v:.2e}")
    }
}

pub(crate) fn pill(
    id: &'static str,
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
/// "Chain B, C and O keep their own colours", or nothing when none do.
fn kept_note(names: &[String]) -> Option<String> {
    let list = match names {
        [] => return None,
        [one] => return Some(format!("{one} keeps its own colour")),
        [rest @ .., last] if names.len() <= 4 => format!("{} and {last}", rest.join(", ")),
        _ => format!("{} and {} more", names[..3].join(", "), names.len() - 3),
    };
    Some(format!("{list} keep their own colours"))
}

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

/// A menu card hung from the button it belongs to, just under it or, `up`,
/// just over it, at its own width; moved along to stay inside the window
/// where it would run past an edge, as in a narrow viewport it can.
fn hang(list: gpui::Div, up: bool) -> impl IntoElement {
    let (corner, below) = if up {
        (gpui::Corner::BottomLeft, px(-4.))
    } else {
        (gpui::Corner::TopLeft, px(36.))
    };
    // Held at the button's top-left corner, the point it is measured from;
    // left in the flow it would be under the button, wherever that ends.
    div().absolute().top_0().left_0().child(
        deferred(
            gpui::anchored()
                .anchor(corner)
                .offset(gpui::point(px(-3.), below))
                .snap_to_window_with_margin(px(8.))
                // In the flow of the anchored element, which does the placing.
                .child(list.relative().top_0().left_0()),
        )
        .with_priority(2),
    )
}

/// Keep a press on a floating card to itself. A menu hangs out over the
/// viewport, which closes menus on mouse down; without this the menu would
/// be gone before the click on its item finished.
pub(crate) fn swallow(_: &MouseDownEvent, _: &mut Window, cx: &mut App) {
    cx.stop_propagation();
}

/// One line of a menu: a name, a tick when it is the one in use.
pub(crate) fn item(
    id: impl Into<gpui::ElementId>,
    text: String,
    chosen: bool,
) -> gpui::Stateful<gpui::Div> {
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
            d.child(
                div()
                    .text_size(px(10.))
                    .text_color(theme::muted())
                    .child("✓"),
            )
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

/// The middle of some points, and how far the furthest is from it.
fn spread(points: impl Iterator<Item = [f32; 3]> + Clone) -> ([f32; 3], f32) {
    let (sum, n) = points.clone().fold(([0.0f32; 3], 0usize), |(s, n), p| {
        ([s[0] + p[0], s[1] + p[1], s[2] + p[2]], n + 1)
    });
    if n == 0 {
        return ([0.0; 3], 0.0);
    }
    let c = sum.map(|v| v / n as f32);
    let reach = points
        .map(|p| ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2)).sqrt())
        .fold(0.0, f32::max);
    (c, reach)
}

/// A whole residue picked: by its alpha carbon, or else its first atom.
fn residue_pick(s: &Structure, residue: usize) -> Option<Picked> {
    let r = s.residues.get(residue)?;
    let atom = match r.ca.or(r.nucleic) {
        Some(ca) => ca as usize,
        None => s.atoms.iter().position(|a| a.residue as usize == residue)?,
    };
    Some(Picked {
        atom,
        residue: true,
        depth: 0.0,
    })
}

/// The letter a residue goes by in the sequence: its amino acid's, or a
/// nucleotide's base.
fn residue_letter(r: &crate::protein::Residue) -> char {
    crate::protein::residue_code(r)
}

/// Residues per block of the sequence, each block numbered.
const SEQUENCE_BLOCK: usize = 10;

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
                // With nothing loaded there is nothing to overlay, and Open
                // opens straight away.
                d.when(!loaded, |d| {
                    d.child(
                        tab("open-file", "Open", false)
                            .on_click(cx.listener(|this, _, _, cx| this.prompt_open(cx))),
                    )
                })
                .when(loaded, |d| {
                    d.child(self.tab_with_menu(
                        Dropdown::Open,
                        "open",
                        "Open",
                        self.open_menu(cx).into_any_element(),
                        cx,
                    ))
                })
                .child(self.tab_with_menu(
                    Dropdown::Samples,
                    "samples",
                    "Samples",
                    self.samples_menu(cx).into_any_element(),
                    cx,
                ))
            })
            .when(loaded, |d| {
                d.child(self.tab_with_menu(
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
                .child(self.tab_with_menu(
                    Dropdown::Film,
                    "film",
                    "Film",
                    self.film_menu(cx).into_any_element(),
                    cx,
                ))
                .when(viewing, |d| {
                    d.child(self.tab_with_menu(
                        Dropdown::Save,
                        "save",
                        "Save",
                        self.save_menu(cx).into_any_element(),
                        cx,
                    ))
                })
            })
            .child(div().flex_1())
            .children(self.recording_badge(cx))
            .when(loaded && !self.narrow, |d| d.child(self.angle_readout(cx)))
            .when_some(badge.filter(|_| !self.narrow), |d, badge| {
                d.child(
                    div()
                        .font_family(theme::MONO)
                        .text_size(px(10.))
                        .text_color(theme::faint())
                        .child(badge),
                )
            })
    }

    /// How far the view is turned about each axis from square on, and how
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
                    .on_click(cx.listener(move |this, _, _, cx| this.toggle(which, cx))),
            )
            .when(open, |d| d.child(deferred(card).with_priority(2)))
    }

    /// The structures built into the app, to open at any time.
    fn samples_menu(&self, cx: &mut Context<Self>) -> impl IntoElement {
        // What was opened lately, by the name it would be known by.
        let recent: Vec<(String, String)> = crate::settings::recent()
            .into_iter()
            .take(6)
            .map(|o| {
                let name = match o.strip_prefix("sample:") {
                    Some(file) => SAMPLES
                        .iter()
                        .find(|s| s.file == file)
                        .map_or(file.to_string(), |s| s.name.to_string()),
                    None => std::path::Path::new(&o)
                        .file_name()
                        .map_or(o.clone(), |n| n.to_string_lossy().into_owned()),
                };
                (o, name)
            })
            .collect();
        menu(260.)
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.open = None;
                cx.notify();
            }))
            .when(!recent.is_empty(), |d| {
                d.child(
                    div()
                        .px(px(8.))
                        .pt(px(3.))
                        .pb(px(2.))
                        .child(label("RECENT")),
                )
                .children(recent.into_iter().enumerate().map(|(i, (origin, name))| {
                    div()
                        .id(("recent", i))
                        .debug_selector(move || format!("recent-{i}"))
                        .px(px(8.))
                        .py(px(4.))
                        .rounded(px(5.))
                        .cursor_pointer()
                        .hover(|s| s.bg(theme::sunken()))
                        .text_size(px(11.5))
                        .text_color(theme::text())
                        .child(clip(&name, 34))
                        .on_click(cx.listener(move |this, _, _, cx| {
                            this.open = None;
                            this.reopen(&origin, cx);
                        }))
                }))
                .child(divide())
            })
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

        // A group too wide for a phone-shaped window wraps onto a second
        // line rather than running off both sides.
        let group = |children: Vec<gpui::AnyElement>| {
            div()
                .flex()
                .flex_wrap()
                .justify_center()
                .max_w_full()
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

        // In a viewport as narrow as a phone the rows will not fit across,
        // so each becomes one button naming the choice made, with the rest
        // in a menu under it.
        if self.compact {
            let reprs: Vec<gpui::AnyElement> = Repr::ALL
                .iter()
                .map(|r| {
                    let r = *r;
                    let enabled = !r.needs_polymer() || polymer;
                    item(("repr-item", r as usize), r.label().to_string(), r == repr)
                        .debug_selector(move || format!("repr-item-{}", r.label()))
                        .when(!enabled, |d| d.opacity(0.4))
                        .on_click(cx.listener(move |this, _, _, cx| {
                            this.open = None;
                            this.set_repr(r, cx);
                        }))
                        .into_any_element()
                })
                .collect();
            let styles: Vec<gpui::AnyElement> = Style::ALL
                .iter()
                .map(|st| {
                    let st = *st;
                    item(
                        ("style-item", st as usize),
                        st.label().to_string(),
                        st == style,
                    )
                    .on_click(cx.listener(move |this, _, _, cx| {
                        this.open = None;
                        this.set_style(st, cx);
                    }))
                    .into_any_element()
                })
                .collect();
            let colour_name = |c: Coloring| match c {
                Coloring::Chain if !polymer => "Molecule",
                _ => c.label(),
            };
            let colourings: Vec<gpui::AnyElement> = options
                .iter()
                .filter(|c| **c != Coloring::Confidence || has_confidence)
                .filter(|c| {
                    **c != Coloring::BFactor || structure.is_some_and(|s| s.b_span().is_some())
                })
                .filter(|c| **c != Coloring::Structure || polymer)
                .filter(|c| {
                    **c != Coloring::Property || structure.is_some_and(|s| !s.props.is_empty())
                })
                .map(|c| {
                    let c = *c;
                    item(
                        ("colour-item", c as usize),
                        colour_name(c).to_string(),
                        c == coloring,
                    )
                    .debug_selector(move || format!("colour-item-{}", c.label()))
                    .on_click(cx.listener(move |this, _, _, cx| {
                        this.open = None;
                        this.set_coloring(c, cx);
                    }))
                    .into_any_element()
                })
                .collect();
            let top = bar()
                .debug_selector(|| "top-controls".into())
                .top(px(12.))
                .child(self.choice(Dropdown::Repr, "repr-menu", repr.label(), false, reprs, cx))
                .child(self.choice(
                    Dropdown::Style,
                    "style-menu",
                    style.label(),
                    false,
                    styles,
                    cx,
                ))
                .child(group(vec![pill("measure", "Measure", self.measuring)
                    .debug_selector(|| "measure".into())
                    .cursor_pointer()
                    .on_click(cx.listener(|this, _, _, cx| this.toggle_measuring(cx)))
                    .into_any_element()]))
                .when(has_cell, |d| {
                    d.child(group(vec![pill("cell", "Cell", self.show_cell)
                        .cursor_pointer()
                        .on_click(cx.listener(|this, _, _, cx| {
                            this.show_cell = !this.show_cell;
                            this.rebuild(cx);
                        }))
                        .into_any_element()]))
                });
            let bottom = bar()
                .bottom(px(bottom))
                .child(self.palette_button(cx))
                .child(self.choice(
                    Dropdown::Coloring,
                    "colour-menu",
                    colour_name(coloring),
                    true,
                    colourings,
                    cx,
                ));
            return vec![top.into_any_element(), bottom.into_any_element()];
        }

        // Up top, how it is drawn; at the foot, what it is coloured by
        // beside the palette it is coloured from.
        let top = bar()
            .debug_selector(|| "top-controls".into())
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
            .child(group(vec![pill("measure", "Measure", self.measuring)
                .debug_selector(|| "measure".into())
                .cursor_pointer()
                .on_click(cx.listener(|this, _, _, cx| this.toggle_measuring(cx)))
                .into_any_element()]))
            .when(has_cell, |d| {
                d.child(group(vec![pill("cell", "Cell", self.show_cell)
                    .cursor_pointer()
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.show_cell = !this.show_cell;
                        this.rebuild(cx);
                    }))
                    .into_any_element()]))
            });

        let bottom = bar()
            .bottom(px(bottom))
            .child(self.palette_button(cx))
            .child(group(
                options
                    .iter()
                    // pLDDT means nothing without a predicted model, so it
                    // is offered only for one; B-factors only where the file
                    // has them, and helices only in a protein.
                    .filter(|c| **c != Coloring::Confidence || has_confidence)
                    .filter(|c| {
                        **c != Coloring::BFactor || structure.is_some_and(|s| s.b_span().is_some())
                    })
                    .filter(|c| **c != Coloring::Structure || polymer)
                    .filter(|c| {
                        **c != Coloring::Property || structure.is_some_and(|s| !s.props.is_empty())
                    })
                    .map(|c| {
                        let c = *c;
                        let name = match c {
                            Coloring::Chain if !polymer => "Molecule",
                            _ => c.label(),
                        };
                        pill_keyed(("colour", c as usize), name, c == coloring)
                            .cursor_pointer()
                            .on_click(cx.listener(move |this, _, _, cx| this.set_coloring(c, cx)))
                            .into_any_element()
                    })
                    .collect(),
            ));

        vec![top.into_any_element(), bottom.into_any_element()]
    }

    /// A choice made from a menu: a button naming what is chosen, the size
    /// of the palette button beside it, with the menu opening away from the
    /// viewport's edge it sits on.
    fn choice(
        &self,
        which: Dropdown,
        id: &'static str,
        label: &'static str,
        up: bool,
        items: Vec<gpui::AnyElement>,
        cx: &mut Context<Self>,
    ) -> impl IntoElement {
        let open = self.open == Some(which);
        let list = menu(156.)
            .debug_selector(move || format!("{id}-list"))
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.open = None;
                cx.notify();
            }))
            .children(items);
        let arrow = if open == up { "▾" } else { "▴" };
        div()
            .relative()
            .child(
                div()
                    .id(id)
                    .debug_selector(move || id.into())
                    .flex()
                    .items_center()
                    .gap(px(7.))
                    .h(px(32.))
                    .pl(px(12.))
                    .pr(px(10.))
                    .rounded(px(8.))
                    .bg(theme::surface())
                    .border_1()
                    .border_color(theme::border())
                    .shadow_sm()
                    .cursor_pointer()
                    .when(open, |d| d.bg(theme::sunken()))
                    .hover(|s| s.bg(theme::sunken()))
                    .child(
                        div()
                            .text_size(px(11.))
                            .text_color(theme::ink())
                            .child(label),
                    )
                    .child(
                        div()
                            .text_size(px(10.))
                            .text_color(theme::muted())
                            .child(arrow),
                    )
                    .on_click(cx.listener(move |this, _, _, cx| this.toggle(which, cx))),
            )
            .when(open, |d| d.child(hang(list, up)))
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
                    d.child(
                        div()
                            .text_size(px(10.))
                            .text_color(theme::muted())
                            .child("✓"),
                    )
                })
                .on_click(cx.listener(move |this, _, _, cx| {
                    this.palette = p;
                    this.open = None;
                    this.rebuild(cx);
                    this.say_kept(cx);
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
        let top =
            div()
                .flex()
                .flex_col()
                .gap(px(1.))
                .child(heading(PaletteGroup::Classic, true))
                .child(div().flex().gap(px(4.)).children(
                    members(PaletteGroup::Classic).map(|p| div().flex_1().child(item(p))),
                ));
        let left = column(&[PaletteGroup::Perceptual]);
        let right = column(&[PaletteGroup::Thematic, PaletteGroup::Muted]);
        let list = menu(372.)
            .top_auto()
            .bottom(px(34.))
            .debug_selector(|| "palette-menu".into())
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.open = None;
                cx.notify();
            }))
            .child(top)
            .child(
                div()
                    .flex()
                    .items_start()
                    .gap(px(4.))
                    .child(left)
                    .child(right),
            )
            .when_some(kept_note(&self.own_colors()), |d, note| {
                // Picked colours outlast a change of palette; say which, so
                // the ones that did not move are not a surprise.
                d.child(
                    div()
                        .mt(px(6.))
                        .mx(px(8.))
                        .pt(px(8.))
                        .pb(px(2.))
                        .border_t_1()
                        .border_color(theme::border())
                        .flex()
                        .items_center()
                        .gap(px(8.))
                        .child(
                            div()
                                .flex_1()
                                .text_size(px(10.5))
                                .text_color(theme::muted())
                                .child(note),
                        )
                        .child(
                            div()
                                .id("reset-colours")
                                .debug_selector(|| "reset-colours".into())
                                .text_size(px(10.5))
                                .text_color(theme::muted())
                                .cursor_pointer()
                                .hover(|s| s.text_color(theme::ink()))
                                .child("Reset")
                                .on_click(cx.listener(|this, _, _, cx| this.reset_colors(cx))),
                        ),
                )
            });
        // In a narrow viewport the menu is moved along, if it must be, to
        // stay inside the window.
        let list = open.then(|| {
            if self.compact {
                hang(list, true).into_any_element()
            } else {
                deferred(list).with_priority(2).into_any_element()
            }
        });
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
                    .cursor_pointer()
                    .hover(|s| s.bg(theme::sunken()))
                    .tooltip(|_, cx| cx.new(|_| Tip("Palette".into())).into())
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
                    .on_click(cx.listener(|this, _, _, cx| this.toggle(Dropdown::Palette, cx))),
            )
            .children(list)
            .when_some(self.kept.clone().filter(|_| !open), |d, note| {
                d.child(
                    deferred(
                        div()
                            .id("kept-notice")
                            .debug_selector(|| "kept-notice".into())
                            .absolute()
                            .bottom(px(34.))
                            .left_0()
                            .flex()
                            .items_center()
                            .gap(px(10.))
                            .px(px(10.))
                            .py(px(6.))
                            .rounded(px(8.))
                            .bg(theme::surface())
                            .border_1()
                            .border_color(theme::border())
                            .shadow_md()
                            .whitespace_nowrap()
                            .on_mouse_down(MouseButton::Left, swallow)
                            .child(
                                div()
                                    .text_size(px(11.))
                                    .text_color(theme::text())
                                    .child(note),
                            )
                            .child(
                                div()
                                    .id("kept-reset")
                                    .debug_selector(|| "kept-reset".into())
                                    .text_size(px(11.))
                                    .text_color(theme::ink())
                                    .cursor_pointer()
                                    .hover(|s| s.underline())
                                    .child("Reset")
                                    .on_click(cx.listener(|this, _, _, cx| this.reset_colors(cx))),
                            ),
                    )
                    .with_priority(2),
                )
            })
    }

    /// A small button on a row with a colour of its own, giving it back
    /// the palette's. It stands in the margin left of the row's dot, so the
    /// buttons of the chains and of the atoms line up whatever their names.
    fn unpick(&self, row: usize, target: Pick, cx: &mut Context<Self>) -> impl IntoElement {
        let id = match target {
            Pick::Chain(_) => "chain-unpick",
            Pick::Ligand(_) => "ligand-unpick",
            Pick::Species(_) => "species-unpick",
        };
        let button = div()
            .id((id, row))
            .debug_selector(move || format!("{id}-{row}"))
            .flex_none()
            .size(px(16.))
            .flex()
            .items_center()
            .justify_center()
            .rounded(px(4.))
            .text_size(px(11.))
            .line_height(px(16.))
            .text_color(theme::faint())
            .cursor_pointer()
            .hover(|s| s.text_color(theme::ink()).bg(theme::sunken()))
            .child("↺")
            .tooltip(|_, cx| cx.new(|_| Tip("Reset to the palette colour".into())).into())
            .on_click(cx.listener(move |this, _, _, cx| this.paint(target.clone(), None, cx)));
        div()
            .absolute()
            .left(px(-20.))
            .top_0()
            .bottom_0()
            .flex()
            .items_center()
            .child(button)
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
            .child(
                item("spin", "Spin".into(), self.turning == Some(Turn::Spin))
                    .debug_selector(|| "spin".into())
                    .on_click(cx.listener(|this, _, _, cx| this.toggle_turn(Turn::Spin, cx))),
            )
            .child(
                item("rock", "Rock".into(), self.turning == Some(Turn::Rock))
                    .debug_selector(|| "rock".into())
                    .on_click(cx.listener(|this, _, _, cx| this.toggle_turn(Turn::Rock, cx))),
            )
            .child(divide())
            .child(
                item("perspective", "Perspective".into(), self.camera.perspective)
                    .debug_selector(|| "perspective".into())
                    .on_click(cx.listener(|this, _, _, cx| this.toggle_perspective(cx))),
            )
            .child(
                item("slab", "Clipping slab".into(), self.camera.slab.is_some())
                    .debug_selector(|| "slab".into())
                    .on_click(cx.listener(|this, _, _, cx| this.toggle_slab(cx))),
            )
            .child(divide())
            .child(
                item("save-view", "Save this viewpoint".into(), false)
                    .debug_selector(|| "save-view".into())
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.views.push(this.camera);
                        this.note(format!("Saved as view {}", this.views.len()));
                        cx.notify();
                    })),
            )
            .children(self.views.iter().enumerate().map(|(i, &view)| {
                item(("view", i), format!("View {}", i + 1), false)
                    .debug_selector(move || format!("view-{i}"))
                    .on_click(cx.listener(move |this, _, _, cx| {
                        this.open = None;
                        this.fly_to_view(view, cx);
                    }))
            }))
            .when(!self.views.is_empty(), |d| {
                d.child(
                    item("clear-views", "Clear saved views".into(), false)
                        .debug_selector(|| "clear-views".into())
                        .on_click(cx.listener(|this, _, _, cx| {
                            this.views.clear();
                            cx.notify();
                        })),
                )
            })
            .child(divide())
            .child(
                item("recentre", "Centre".into(), false)
                    .on_click(cx.listener(|this, _, _, cx| this.recentre(cx))),
            )
            .when(self.selected.is_some(), |d| {
                d.child(
                    item("focus-menu", "Focus on the pick".into(), false)
                        .debug_selector(|| "focus-menu".into())
                        .on_click(cx.listener(|this, _, _, cx| {
                            this.open = None;
                            this.focus_selected(cx);
                        })),
                )
            })
            .child(
                item("reset-view", "Reset view".into(), false)
                    .debug_selector(|| "reset-view".into())
                    .on_click(cx.listener(|this, _, _, cx| this.reset_view(cx))),
            )
    }

    /// Start recording a camera path, or stop the one being recorded.
    pub fn toggle_recording(&mut self, cx: &mut Context<Self>) {
        self.open = None;
        if self.recording() {
            self.stop_recording(cx);
        } else {
            self.start_recording(cx);
        }
    }

    /// A film of the view, or the camera path to make one along.
    fn film_menu(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let film = item("movie-open", "Export film…".into(), false)
            .debug_selector(|| "movie-open".into())
            .on_click(cx.listener(|this, _, _, cx| this.open_movie(cx)));
        let record = item(
            "movie-record-menu",
            if self.recording() {
                "Stop recording".into()
            } else {
                "Record a camera path".into()
            },
            false,
        )
        .debug_selector(|| "movie-record-menu".into())
        .on_click(cx.listener(|this, _, _, cx| this.toggle_recording(cx)));
        menu(200.)
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.open = None;
                cx.notify();
            }))
            .child(film)
            .child(record)
    }

    /// How far a ray trace has got, the view's or a saved image's, in the
    /// corner of the viewport while one runs.
    fn tracing_label(&self) -> Option<impl IntoElement> {
        let (what, part) = match &self.saving {
            Some((done, rows)) => (
                "Ray tracing the image",
                done.load(std::sync::atomic::Ordering::Relaxed) as f32 / (*rows).max(1) as f32,
            ),
            None => ("Ray tracing", self.cache.borrow().tracing()?),
        };
        Some(
            div()
                .absolute()
                .bottom(px(24.))
                .right(px(20.))
                .px(px(8.))
                .py(px(3.))
                .rounded(px(5.))
                .bg(theme::surface())
                .border_1()
                .border_color(theme::border())
                .font_family(theme::MONO)
                .text_size(px(10.))
                .text_color(theme::muted())
                .child(format!("{what} · {}%", (part * 100.0).floor() as u32)),
        )
    }

    /// Write the viewport out as a picture, at a multiple of its size.
    fn image_menu(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let (w, h) = self.cache.borrow().size;
        let transparent = item(
            "transparent",
            "Transparent background".into(),
            self.transparent,
        )
        .debug_selector(|| "transparent".into())
        .on_click(cx.listener(|this, _, _, cx| {
            this.transparent = !this.transparent;
            cx.notify();
        }));
        let traced = item("traced", "Ray traced · slower".into(), self.traced)
            .debug_selector(|| "traced".into())
            .tooltip(|_, cx| {
                cx.new(|_| {
                    Tip("Shadows, ambient occlusion and smooth edges, in the view once it holds still and in saved images. Off by default: a traced picture takes seconds where a plain one takes milliseconds, a traced film minutes".into())
                })
                .into()
            })
            .on_click(cx.listener(|this, _, _, cx| {
                this.traced = !this.traced;
                cx.notify();
            }));
        menu(200.)
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.open = None;
                cx.notify();
            }))
            .child(traced)
            .child(divide())
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
            .child(divide())
            .child(transparent)
    }

    /// The logo and what the tab is for, at the head of the panel.
    fn masthead(&self) -> impl IntoElement {
        // Each name at the size that spans the panel beside the logo.
        let (name, size, purpose) = match self.mode {
            Mode::View => ("profold", 66., "Proteins, molecules and crystals"),
            Mode::Fold => (
                "Fold a protein",
                34.,
                "Structure prediction with Apple's SimpleFold",
            ),
        };
        div()
            .flex()
            .items_center()
            .gap(px(12.))
            .px(px(20.))
            .child(gpui::img(self.mark.clone()).size(px(88.)).flex_none())
            .child(
                div()
                    .flex_1()
                    .min_w_0()
                    .flex()
                    .flex_col()
                    .gap(px(8.))
                    .child(
                        div()
                            .font_family(theme::SERIF)
                            .font_weight(gpui::FontWeight::BOLD)
                            .text_size(px(size))
                            .line_height(px(size + 2.))
                            .whitespace_nowrap()
                            .text_color(theme::ink())
                            .child(name),
                    )
                    .child(
                        div()
                            .text_size(px(12.))
                            .font_weight(gpui::FontWeight::MEDIUM)
                            .text_color(theme::text())
                            .child(purpose),
                    ),
            )
    }

    fn search(&self, focused: bool, cx: &mut Context<Self>) -> impl IntoElement {
        let detected = self.detected.clone();
        let folding = self.mode == Mode::Fold;
        let ready = if folding {
            matches!(
                detected,
                Input::Sequence(_) | Input::Pdb(_) | Input::UniProt(_)
            )
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

        let empty = detected == Input::Empty;
        // What the field holds, named inside it, beside the key that acts on it.
        let kind = match &detected {
            Input::Sequence(s) => format!("Sequence · {} aa", s.len()),
            Input::Pdb(_) => "PDB ID".into(),
            Input::UniProt(_) => "UniProt".into(),
            Input::Cid(_) => "PubChem".into(),
            Input::File(_) | Input::Missing(_) => "File".into(),
            Input::Query(_) if folding => "Not foldable".into(),
            Input::Query(_) => "Search".into(),
            Input::Unknown => "Unrecognised".into(),
            Input::Empty => String::new(),
        };
        let run = div()
            .id("run")
            .debug_selector(|| "run".into())
            .flex_none()
            .h(px(26.))
            .px(px(9.))
            .flex()
            .items_center()
            .gap(px(5.))
            .rounded(px(6.))
            .text_size(px(11.))
            .when(ready && !busy, |d| {
                d.bg(theme::ink())
                    .text_color(theme::on_ink())
                    .hover(|s| s.opacity(0.86))
                    .cursor_pointer()
            })
            .when(!ready || busy, |d| {
                d.border_1()
                    .border_color(theme::border())
                    .text_color(theme::faint())
            })
            .when(!empty, |d| d.child(div().text_size(px(11.)).child(action)))
            .child("↵")
            .on_click(cx.listener(|this, _, _, cx| this.submit(cx)));
        // A name search says so on its button; the chip would only repeat it.
        let chip = (!empty && !(matches!(detected, Input::Query(_)) && !folding)).then(|| {
            div()
                .flex_none()
                .px(px(7.))
                .py(px(2.))
                .rounded(px(4.))
                .bg(theme::sunken())
                .text_size(px(10.5))
                .text_color(if ready { theme::text() } else { theme::muted() })
                .child(kind)
        });

        let hint = |text: String| {
            div()
                .min_w_0()
                .truncate()
                .text_size(px(10.5))
                .text_color(theme::faint())
                .child(text)
                .into_any_element()
        };
        let choosing = self.dropdown && !self.hits.is_empty();
        let status = match &detected {
            Input::Empty if folding => Some(hint("A sequence, PDB ID or UniProt accession".into())),
            Input::Empty => Some(hint(
                "Names, genes, PDB and UniProt IDs, CIDs, sequences or files".into(),
            )),
            Input::Query(_) if folding => Some(hint("Not a sequence, PDB ID or accession".into())),
            Input::Query(_) if choosing => {
                Some(hint("↑↓ to choose · ↵ to open · esc to close".into()))
            }
            Input::Query(_) => Some(hint("Searches the PDB, AlphaFold DB and PubChem".into())),
            Input::Sequence(_) if !folding => Some(hint("Looked up in the PDB by sequence".into())),
            Input::Unknown | Input::Missing(_) | Input::File(_) => Some(hint(detected.detail())),
            _ if comparable => Some(hint("Fold it, or fold and compare with the deposit".into())),
            _ => None,
        }
        .map(|h| {
            div()
                .flex()
                .items_center()
                .gap(px(6.))
                .min_w_0()
                .child(h)
                .when(self.searching.is_some(), |d| {
                    d.child(
                        div()
                            .flex_none()
                            .text_size(px(10.))
                            .text_color(theme::faint())
                            .child("searching…"),
                    )
                })
        });

        div()
            .relative()
            .flex()
            .flex_col()
            .gap(px(10.))
            .px(px(20.))
            .child(
                div().flex().items_center().gap(px(8.)).child(
                    div()
                        .id("search-field")
                        .flex_1()
                        .min_w_0()
                        .flex()
                        .items_center()
                        .gap(px(6.))
                        .h(px(40.))
                        .pl(px(12.))
                        .pr(px(7.))
                        .rounded(px(9.))
                        .bg(if focused {
                            theme::surface()
                        } else {
                            theme::sunken()
                        })
                        .border_1()
                        .border_color(if focused {
                            theme::border_focus()
                        } else {
                            theme::border()
                        })
                        .text_size(px(12.5))
                        .text_color(theme::text())
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
                        .child(
                            div()
                                .flex_1()
                                .min_w_0()
                                // Typed text is set in mono, for sequences and
                                // IDs; the placeholder reads as prose.
                                .when(!empty, |d| d.font_family(theme::MONO))
                                .child(self.input.clone()),
                        )
                        .when_some(chip, |d, c| d.child(c))
                        .child(run),
                ),
            )
            // While the next entry is to be laid over this one, say so.
            .when(self.overlay_next && !busy, |d| {
                let onto = self
                    .structure
                    .as_ref()
                    .map(|s| s.id.to_uppercase())
                    .unwrap_or_default();
                d.child(
                    div()
                        .flex()
                        .items_center()
                        .justify_between()
                        .gap(px(8.))
                        .text_size(px(10.5))
                        .child(div().text_color(theme::text()).child(format!(
                            "The next structure opened is superimposed onto {onto}"
                        )))
                        .child(
                            div()
                                .id("overlay-cancel")
                                .debug_selector(|| "overlay-cancel".into())
                                .px(px(6.))
                                .rounded(px(4.))
                                .text_color(theme::muted())
                                .cursor_pointer()
                                .hover(|d| d.bg(theme::sunken()).text_color(theme::ink()))
                                .child("Cancel")
                                .on_click(cx.listener(|this, _, _, cx| {
                                    this.overlay_next = false;
                                    cx.notify();
                                })),
                        ),
                )
            })
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
        div().px(px(20.)).when_some(self.error.clone(), |d, err| {
            d.child(
                div()
                    .text_size(px(11.))
                    .text_color(theme::danger())
                    .child(err),
            )
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
            .when(!folding, |d| {
                d.child(self.sample_list(cx)).child(div().h(px(8.)))
            })
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
                )
                .child(
                    div()
                        .id("comparison-close")
                        .debug_selector(|| "comparison-close".into())
                        .px(px(5.))
                        .rounded(px(4.))
                        .text_size(px(12.))
                        .text_color(theme::faint())
                        .cursor_pointer()
                        .hover(|d| d.bg(theme::sunken()).text_color(theme::ink()))
                        .child("×")
                        .tooltip(|_, cx| cx.new(|_| Tip("Show the reference only".into())).into())
                        .on_mouse_down(MouseButton::Left, |_, _, cx| cx.stop_propagation())
                        .on_click(cx.listener(|this, _, _, cx| this.remove_overlay(cx))),
                ),
        )
    }

    /// The deposited structure, and each assembly that is more than it.
    fn assembly_choice(
        &self,
        shown: &Structure,
        cx: &mut Context<Self>,
    ) -> Option<gpui::AnyElement> {
        let base: &Structure = self.deposited.as_deref().unwrap_or(shown);
        let offered: Vec<(usize, String)> = base
            .assemblies
            .iter()
            .enumerate()
            .filter(|(_, a)| !a.is_deposited(base))
            .map(|(i, a)| (i, a.label(base)))
            .collect();
        if offered.is_empty() || self.compare.is_some() {
            return None;
        }
        // Two assemblies of one name are told apart by number.
        let named = |i: usize, label: &str| {
            if offered.iter().filter(|(_, l)| l == label).count() > 1 {
                format!("{label} {}", base.assemblies[i].id)
            } else {
                label.to_string()
            }
        };
        let option = |k: Option<usize>, label: String| {
            let chosen = self.assembly == k;
            div()
                .id(("assembly", k.map_or(0, |k| k + 1)))
                .debug_selector(move || {
                    format!(
                        "assembly-{}",
                        k.map_or("deposited".to_string(), |k| k.to_string())
                    )
                })
                .px(px(7.))
                .py(px(1.))
                .rounded(px(5.))
                .text_size(px(10.))
                .when(chosen, |d| {
                    d.bg(theme::surface()).text_color(theme::ink()).shadow_sm()
                })
                .when(!chosen, |d| {
                    d.text_color(theme::muted())
                        .cursor_pointer()
                        .hover(|s| s.text_color(theme::ink()))
                })
                .child(label)
                .on_click(cx.listener(move |this, _, _, cx| this.show_assembly(k, cx)))
        };
        Some(
            div()
                .pt(px(14.))
                .flex()
                .items_center()
                .justify_between()
                .gap(px(8.))
                .child(key("assembly"))
                .child(
                    div()
                        .flex()
                        .flex_wrap()
                        .p(px(2.))
                        .gap(px(2.))
                        .rounded(px(7.))
                        .bg(theme::sunken())
                        .child(option(None, "Deposited".into()))
                        .children(offered.iter().map(|(i, l)| option(Some(*i), named(*i, l)))),
                )
                .into_any_element(),
        )
    }

    /// Bring the plots' worked-out parts up to the structure shown.
    fn plots_for(
        &self,
        structure: &Structure,
    ) -> (Arc<Vec<Option<(f32, f32)>>>, Option<Arc<gpui::RenderImage>>) {
        let key = (structure as *const Structure as usize, structure.frame);
        let mut cache = self.plot_cache.borrow_mut();
        if cache.key != Some(key) {
            cache.key = Some(key);
            cache.torsions = Arc::new(structure.torsions());
            cache.pae = structure.pae.as_ref().map(|pae| {
                // AlphaFold DB's colours: dark green where the error is
                // small, through to white at the top of the scale.
                let n = pae.n.clamp(1, 600);
                let step = pae.n as f32 / n as f32;
                let low = [0.05, 0.30, 0.15];
                picture(n, n, |x, y| {
                    let v = pae.at(
                        ((y as f32 * step) as usize).min(pae.n - 1),
                        ((x as f32 * step) as usize).min(pae.n - 1),
                    );
                    let t = (v / pae.max).clamp(0.0, 1.0);
                    [0, 1, 2].map(|i| low[i] + (1.0 - low[i]) * t)
                })
            });
        }
        (cache.torsions.clone(), cache.pae.clone())
    }

    /// φ against ψ, a dot per residue in its secondary structure's colour.
    /// A dot pointed at glows its residue; clicked, it picks it out; what is
    /// picked is ringed.
    fn ramachandran(
        &self,
        structure: &Structure,
        cx: &mut Context<Self>,
    ) -> Option<gpui::AnyElement> {
        if !structure.residues.iter().any(|r| r.ca.is_some()) {
            return None;
        }
        let (torsions, _) = self.plots_for(structure);
        if torsions.iter().all(Option::is_none) {
            return None;
        }
        if !self.rama_open {
            return Some(
                self.plot_row("ramachandran", false, None, cx)
                    .into_any_element(),
            );
        }
        let selected = self.selected_residue();
        let pointed = self
            .rama_hover
            .or(self.seq_hover)
            .or(self.hovered_residue());
        let dots: Vec<(usize, f32, f32, [f32; 3])> = torsions
            .iter()
            .enumerate()
            .filter_map(|(r, t)| {
                let (phi, psi) = (*t)?;
                let res = &structure.residues[r];
                let color = if res.name == "GLY" {
                    [0.55, 0.56, 0.58]
                } else {
                    render::SS_COLORS
                        .iter()
                        .find(|c| c.0 == res.ss)
                        .map_or([0.6; 3], |c| c.2)
                };
                Some((r, (phi + 180.0) / 360.0, (180.0 - psi) / 360.0, color))
            })
            .collect();
        let store = self.rama_bounds.clone();
        let ground = ramachandran_ground();
        let caption = pointed.and_then(|r| {
            let (phi, psi) = torsions.get(r).copied().flatten()?;
            let res = &structure.residues[r];
            Some(format!(
                "{} {} · φ {phi:.0}° ψ {psi:.0}°",
                res.name, res.seq
            ))
        });
        let ink = theme::ink();
        let plot = div()
            .id("ramachandran")
            .debug_selector(|| "ramachandran".into())
            .relative()
            .w(px(PLOT))
            .h(px(PLOT))
            .rounded(px(4.))
            .overflow_hidden()
            .border_1()
            .border_color(theme::border())
            .child(
                gpui::img(gpui::ImageSource::Render(ground))
                    .absolute()
                    .inset_0()
                    .size_full(),
            )
            .child({
                let dots = dots.clone();
                canvas(
                    move |bounds, _, _| store.set(bounds),
                    move |bounds, _, window, _| {
                        let (w, h) = (f32::from(bounds.size.width), f32::from(bounds.size.height));
                        // The axes through zero.
                        let line = gpui::Hsla { a: 0.25, ..ink };
                        window.paint_quad(gpui::fill(
                            Bounds::new(
                                bounds.origin + gpui::point(px(w * 0.5), px(0.)),
                                gpui::size(px(1.), px(h)),
                            ),
                            line,
                        ));
                        window.paint_quad(gpui::fill(
                            Bounds::new(
                                bounds.origin + gpui::point(px(0.), px(h * 0.5)),
                                gpui::size(px(w), px(1.)),
                            ),
                            line,
                        ));
                        for &(r, x, y, c) in &dots {
                            let big = Some(r) == selected || Some(r) == pointed;
                            let d = if big { 7.0 } else { 4.0 };
                            let at = bounds.origin
                                + gpui::point(px(x * w - d * 0.5), px(y * h - d * 0.5));
                            let mut q =
                                gpui::fill(Bounds::new(at, gpui::size(px(d), px(d))), swatch(c))
                                    .corner_radii(px(d * 0.5));
                            if big {
                                q = q.border_widths(px(1.5)).border_color(ink);
                            }
                            window.paint_quad(q);
                        }
                    },
                )
                .absolute()
                .inset_0()
                .size_full()
            })
            .on_mouse_move(cx.listener(move |this, e: &MouseMoveEvent, _, cx| {
                let b = this.rama_bounds.get();
                let (w, h) = (f32::from(b.size.width), f32::from(b.size.height));
                let (x, y) = (
                    f32::from(e.position.x - b.origin.x),
                    f32::from(e.position.y - b.origin.y),
                );
                let near = dots
                    .iter()
                    .map(|&(r, dx, dy, _)| (r, (dx * w - x).powi(2) + (dy * h - y).powi(2)))
                    .filter(|(_, d)| *d < 36.0)
                    .min_by(|a, b| a.1.total_cmp(&b.1))
                    .map(|(r, _)| r);
                if near != this.rama_hover {
                    this.rama_hover = near;
                    this.set_seq_hover(near, cx);
                    cx.notify();
                }
            }))
            .on_hover(cx.listener(|this, hovered: &bool, _, cx| {
                if !*hovered && this.rama_hover.take().is_some() {
                    this.set_seq_hover(None, cx);
                    cx.notify();
                }
            }))
            .on_click(cx.listener(|this, _, _, cx| {
                let Some(r) = this.rama_hover else { return };
                let Some(s) = this.structure.as_ref() else {
                    return;
                };
                let target = residue_pick(s, r).map(|picked| Target { which: 0, picked });
                this.select(target, cx);
            }));
        Some(
            div()
                .flex()
                .flex_col()
                .child(self.plot_row(
                    "ramachandran",
                    true,
                    Some(caption.unwrap_or_else(|| "φ across, ψ up".into())),
                    cx,
                ))
                .child(div().pt(px(10.)).child(plot))
                .into_any_element(),
        )
    }

    /// A predicted model's expected error between each pair of residues, as
    /// AlphaFold DB draws it; pointed at, it names the pair.
    fn pae_plot(&self, structure: &Structure, cx: &mut Context<Self>) -> Option<gpui::AnyElement> {
        let pae = structure.pae.clone()?;
        let (_, image) = self.plots_for(structure);
        let image = image?;
        if !self.pae_open {
            return Some(
                self.plot_row("predicted aligned error", false, None, cx)
                    .into_any_element(),
            );
        }
        let store = self.pae_bounds.clone();
        let polymer: Vec<usize> = (0..structure.residues.len())
            .filter(|&r| !structure.residues[r].outside_polymer())
            .collect();
        let caption = self.pae_hover.and_then(|(i, j)| {
            let (a, b) = (
                &structure.residues[*polymer.get(i)?],
                &structure.residues[*polymer.get(j)?],
            );
            Some(format!(
                "{}{} on {}{} · {:.1} Å",
                b.name,
                b.seq,
                a.name,
                a.seq,
                pae.at(i, j)
            ))
        });
        let n = pae.n;
        Some(
            div()
                .flex()
                .flex_col()
                .child(self.plot_row(
                    "predicted aligned error",
                    true,
                    Some(caption.unwrap_or_else(|| format!("0 – {:.0} Å", pae.max))),
                    cx,
                ))
                .child(
                    div().pt(px(10.)).child(
                        div()
                            .id("pae")
                            .debug_selector(|| "pae".into())
                            .relative()
                            .w(px(PLOT))
                            .h(px(PLOT))
                            .rounded(px(4.))
                            .overflow_hidden()
                            .child(
                                gpui::img(gpui::ImageSource::Render(image))
                                    .absolute()
                                    .inset_0()
                                    .size_full(),
                            )
                            .child(
                                canvas(move |bounds, _, _| store.set(bounds), |_, _, _, _| {})
                                    .absolute()
                                    .inset_0()
                                    .size_full(),
                            )
                            .on_mouse_move(cx.listener(move |this, e: &MouseMoveEvent, _, cx| {
                                let b = this.pae_bounds.get();
                                let (w, h) = (f32::from(b.size.width), f32::from(b.size.height));
                                if w <= 0.0 || h <= 0.0 {
                                    return;
                                }
                                let x =
                                    (f32::from(e.position.x - b.origin.x) / w * n as f32) as usize;
                                let y =
                                    (f32::from(e.position.y - b.origin.y) / h * n as f32) as usize;
                                let at = Some((y.min(n - 1), x.min(n - 1)));
                                if at != this.pae_hover {
                                    this.pae_hover = at;
                                    // The residue the error is of glows.
                                    let r = at.and_then(|(_, j)| {
                                        this.structure.as_ref().and_then(|s| {
                                            (0..s.residues.len())
                                                .filter(|&r| !s.residues[r].outside_polymer())
                                                .nth(j)
                                        })
                                    });
                                    this.set_seq_hover(r, cx);
                                    cx.notify();
                                }
                            }))
                            .on_hover(cx.listener(|this, hovered: &bool, _, cx| {
                                if !*hovered && this.pae_hover.take().is_some() {
                                    this.set_seq_hover(None, cx);
                                    cx.notify();
                                }
                            })),
                    ),
                )
                .into_any_element(),
        )
    }

    /// A plot's row in the facts table: its name, what the pointer is over
    /// while it is open, and a button that opens or closes it.
    fn plot_row(
        &self,
        label: &'static str,
        open: bool,
        caption: Option<String>,
        cx: &mut Context<Self>,
    ) -> impl IntoElement {
        let rama = label == "ramachandran";
        let id = if rama {
            "open-ramachandran"
        } else {
            "open-pae"
        };
        div()
            .flex()
            .items_center()
            .justify_between()
            .gap(px(10.))
            .py(px(8.))
            .border_b_1()
            .border_color(theme::border())
            .child(key(label))
            .child(
                div()
                    .flex()
                    .items_center()
                    .gap(px(8.))
                    .children(
                        caption
                            .map(|c| div().text_size(px(9.5)).text_color(theme::faint()).child(c)),
                    )
                    .child(
                        div()
                            .id(id)
                            .debug_selector(move || id.into())
                            .px(px(7.))
                            .py(px(2.))
                            .rounded(px(4.))
                            .text_size(px(10.))
                            .cursor_pointer()
                            .text_color(theme::muted())
                            .hover(|s| s.bg(theme::sunken()).text_color(theme::text()))
                            .child(if open { "Hide" } else { "Show" })
                            .on_click(cx.listener(move |this, _, _, cx| {
                                if rama {
                                    this.rama_open = !this.rama_open;
                                    this.rama_hover = None;
                                } else {
                                    this.pae_open = !this.pae_open;
                                    this.pae_hover = None;
                                }
                                cx.notify();
                            })),
                    ),
            )
    }

    /// The selections: a field to type one in, and a row for each made, with
    /// how many atoms it takes, how it is drawn, its colour, and a ×.
    fn layer_list(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let rows = self.layers.iter().enumerate().map(|(i, layer)| {
            let n = layer.atoms.iter().filter(|a| **a).count();
            let style_pill = |k: render::LayerStyle| {
                let chosen = layer.style == k;
                div()
                    .id(("layer-style", i * 3 + k as usize))
                    .debug_selector(move || format!("layer-{i}-{}", k.label()))
                    .px(px(6.))
                    .py(px(1.))
                    .rounded(px(5.))
                    .text_size(px(10.))
                    .when(chosen, |d| {
                        d.bg(theme::surface()).text_color(theme::ink()).shadow_sm()
                    })
                    .when(!chosen, |d| {
                        d.text_color(theme::muted())
                            .cursor_pointer()
                            .hover(|s| s.text_color(theme::ink()))
                    })
                    .child(k.label())
                    .on_click(cx.listener(move |this, _, _, cx| {
                        if let Some(l) = this.layers.get_mut(i) {
                            l.style = k;
                            this.rebuild(cx);
                        }
                    }))
            };
            let color = layer.color;
            div()
                .flex()
                .flex_col()
                .gap(px(3.))
                .child(
                    div()
                        .flex()
                        .items_center()
                        .gap(px(8.))
                        .child(
                            div()
                                .id(("layer-color", i))
                                .debug_selector(move || format!("layer-{i}-color"))
                                .flex_none()
                                .size(px(13.))
                                .flex()
                                .items_center()
                                .justify_center()
                                .rounded_full()
                                .cursor_pointer()
                                .hover(|d| d.bg(theme::sunken()))
                                .child(
                                    div()
                                        .size(px(8.))
                                        .rounded_full()
                                        .when_some(color, |d, c| d.bg(swatch(c)))
                                        .when(color.is_none(), |d| {
                                            d.border_1().border_color(theme::faint())
                                        }),
                                )
                                .tooltip(|_, cx| {
                                    cx.new(|_| Tip("Element colours or a uniform colour".into()))
                                        .into()
                                })
                                .on_click(cx.listener(move |this, _, _, cx| {
                                    if let Some(l) = this.layers.get_mut(i) {
                                        let at = LAYER_COLORS
                                            .iter()
                                            .position(|c| *c == l.color)
                                            .unwrap_or(0);
                                        l.color = LAYER_COLORS[(at + 1) % LAYER_COLORS.len()];
                                        this.rebuild(cx);
                                    }
                                })),
                        )
                        .child(
                            div()
                                .flex_1()
                                .min_w_0()
                                .font_family(theme::MONO)
                                .text_size(px(10.5))
                                .text_color(theme::text())
                                .child(clip(layer.selection.text(), 30)),
                        )
                        .child(
                            div()
                                .font_family(theme::MONO)
                                .text_size(px(10.))
                                .text_color(theme::faint())
                                .child(thousands(n)),
                        )
                        .child(
                            div()
                                .id(("layer-remove", i))
                                .debug_selector(move || format!("layer-{i}-remove"))
                                .px(px(4.))
                                .rounded(px(3.))
                                .text_size(px(11.))
                                .text_color(theme::faint())
                                .cursor_pointer()
                                .hover(|d| d.bg(theme::sunken()).text_color(theme::ink()))
                                .child("×")
                                .on_click(cx.listener(move |this, _, _, cx| {
                                    if i < this.layers.len() {
                                        this.layers.remove(i);
                                        this.rebuild(cx);
                                    }
                                })),
                        ),
                )
                .child(
                    div().pl(px(21.)).flex().child(
                        div()
                            .flex()
                            .p(px(2.))
                            .gap(px(2.))
                            .rounded(px(7.))
                            .bg(theme::sunken())
                            .children(render::LayerStyle::ALL.map(style_pill)),
                    ),
                )
        });
        div()
            .pt(px(14.))
            .flex()
            .flex_col()
            .gap(px(6.))
            .child(key("selections"))
            .child(
                div()
                    .id("select-field")
                    .debug_selector(|| "select-field".into())
                    // Anywhere in the box puts the cursor in it.
                    .on_mouse_down(
                        MouseButton::Left,
                        cx.listener(|this, _, window, cx| {
                            window.focus(&this.select_input.read(cx).focus_handle(cx))
                        }),
                    )
                    .cursor(gpui::CursorStyle::IBeam)
                    .h(px(28.))
                    .px(px(8.))
                    .flex()
                    .items_center()
                    .rounded(px(6.))
                    .bg(theme::sunken())
                    .border_1()
                    .border_color(theme::border())
                    .font_family(theme::MONO)
                    .text_size(px(11.))
                    .text_color(theme::text())
                    .child(self.select_input.clone()),
            )
            .when_some(self.select_error.clone(), |d, e| {
                d.child(div().text_size(px(10.)).text_color(theme::muted()).child(e))
            })
            .children(rows)
    }

    /// What the colours mean, for the colourings whose colours are not the
    /// chains': a ramp with its ends for B-factors, a swatch each for
    /// helix, strand and coil, and pLDDT's bands.
    fn colour_key(
        &self,
        structure: &Structure,
        cx: &mut Context<Self>,
    ) -> Option<gpui::AnyElement> {
        let row = |color: [f32; 3], label: String| {
            div()
                .flex()
                .items_center()
                .gap(px(6.))
                .child(div().size(px(8.)).rounded_full().bg(swatch(color)))
                .child(
                    div()
                        .text_size(px(10.5))
                        .text_color(theme::text())
                        .child(label),
                )
        };
        let body = match self.coloring {
            Coloring::BFactor => {
                let (lo, hi) = structure.b_span()?;
                let palette = self.palette;
                div()
                    .flex()
                    .items_center()
                    .gap(px(6.))
                    .font_family(theme::MONO)
                    .text_size(px(10.))
                    .text_color(theme::muted())
                    .child(format!("{lo:.0}"))
                    .child(
                        div()
                            .flex_1()
                            .h(px(8.))
                            .flex()
                            .rounded(px(4.))
                            .overflow_hidden()
                            .children((0..24).map(move |i| {
                                div()
                                    .flex_1()
                                    .h_full()
                                    .bg(swatch(palette.ramp(i as f32 / 23.0)))
                            })),
                    )
                    .child(format!("{hi:.0} Å²"))
                    .into_any_element()
            }
            Coloring::Property => {
                let chosen = structure.prop_color.map_or(0, |p| p.index);
                let (lo, hi) = structure.prop_span(chosen)?;
                let palette = self.palette;
                // Which column, as pills that wrap, then its ramp.
                let names = div().flex().flex_wrap().gap(px(4.)).children(
                    structure.props.iter().enumerate().map(|(i, p)| {
                        pill_keyed(("prop", i), p.name.clone(), i == chosen)
                            .cursor_pointer()
                            .on_click(cx.listener(move |this, _, _, cx| {
                                this.colour_by_property(i, None, cx)
                            }))
                    }),
                );
                let ramp = div()
                    .flex()
                    .items_center()
                    .gap(px(6.))
                    .font_family(theme::MONO)
                    .text_size(px(10.))
                    .text_color(theme::muted())
                    .child(figure(lo))
                    .child(
                        div()
                            .flex_1()
                            .h(px(8.))
                            .flex()
                            .rounded(px(4.))
                            .overflow_hidden()
                            .children((0..24).map(move |i| {
                                div()
                                    .flex_1()
                                    .h_full()
                                    .bg(swatch(palette.ramp(i as f32 / 23.0)))
                            })),
                    )
                    .child(figure(hi));
                div()
                    .flex()
                    .flex_col()
                    .gap(px(8.))
                    .child(names)
                    .child(ramp)
                    .into_any_element()
            }
            Coloring::Structure => div()
                .flex()
                .gap(px(14.))
                .children(
                    render::SS_COLORS
                        .iter()
                        .map(|(_, name, c)| row(*c, name.to_string())),
                )
                .into_any_element(),
            Coloring::Hydrophobic | Coloring::Charge => {
                let hydro = self.coloring == Coloring::Hydrophobic;
                let (low, high) = if hydro {
                    (render::HYDROPHILIC, render::HYDROPHOBIC)
                } else {
                    (render::NEGATIVE, render::POSITIVE)
                };
                let ends = if hydro {
                    ("hydrophilic", "hydrophobic")
                } else if self.repr == Repr::Surface {
                    ("− potential", "+")
                } else {
                    ("acidic", "basic")
                };
                div()
                    .flex()
                    .items_center()
                    .gap(px(6.))
                    .text_size(px(10.))
                    .text_color(theme::muted())
                    .child(ends.0)
                    .child(
                        div()
                            .flex_1()
                            .h(px(8.))
                            .flex()
                            .rounded(px(4.))
                            .overflow_hidden()
                            .children((0..24).map(move |i| {
                                div().flex_1().h_full().bg(swatch(render::diverging(
                                    i as f32 / 11.5 - 1.0,
                                    low,
                                    high,
                                )))
                            })),
                    )
                    .child(ends.1)
                    .into_any_element()
            }
            Coloring::Confidence if structure.has_confidence => div()
                .flex()
                .flex_wrap()
                .gap_x(px(14.))
                .gap_y(px(4.))
                .children(
                    render::PLDDT_BANDS
                        .iter()
                        .map(|(name, v)| row(render::plddt_band_color(*v), name.to_string())),
                )
                .into_any_element(),
            _ => return None,
        };
        Some(
            div()
                .pt(px(14.))
                .flex()
                .flex_col()
                .gap(px(6.))
                .child(key(match self.coloring {
                    Coloring::BFactor => "b-factor",
                    Coloring::Structure => "secondary structure",
                    Coloring::Hydrophobic => "hydrophobicity",
                    Coloring::Charge if self.repr == Repr::Surface => "electrostatic potential",
                    Coloring::Charge => "charge",
                    Coloring::Property => "property",
                    _ => "pLDDT",
                }))
                .child(body)
                .into_any_element(),
        )
    }

    /// One row per chain, in the colour that chain is painted. The dots go
    /// hollow when the model is coloured by something other than chain, so the
    /// key never claims a colour the viewport is not using. A dot opens a
    /// picker for that chain's colour.
    fn legend(&self, structure: &Structure, cx: &mut Context<Self>) -> impl IntoElement {
        let by_chain = self.coloring == Coloring::Chain;
        let palette = self.palette;
        let total = structure.chains.len();
        let shown = if self.all_chains { total } else { total.min(6) };
        let rest = total - shown;
        let changed = structure.chains.iter().any(|c| c.color.is_some());

        div()
            .pt(px(14.))
            .flex()
            .flex_col()
            .gap(px(5.))
            .children(
                structure
                    .chains
                    .iter()
                    .take(shown)
                    .enumerate()
                    .map(|(i, chain)| {
                        let name = if chain.id.trim().is_empty() {
                            format!("Chain {}", i + 1)
                        } else {
                            format!("Chain {}", chain.id)
                        };
                        let picking = self.picking == Some(Pick::Chain(i));
                        let default = palette.chain(i);
                        div()
                            .relative()
                            .flex()
                            .items_center()
                            .gap(px(8.))
                            .child(
                                div()
                                    .id(("chain-dot", i))
                                    .debug_selector(move || format!("chain-dot-{i}"))
                                    .flex_none()
                                    .size(px(13.))
                                    .flex()
                                    .items_center()
                                    .justify_center()
                                    .rounded_full()
                                    .cursor_pointer()
                                    .when(chain.color.is_some() && !picking, |d| {
                                        d.border_1().border_color(theme::faint())
                                    })
                                    .when(picking, |d| d.border_1().border_color(theme::ink()))
                                    .hover(|d| d.bg(theme::sunken()))
                                    .child(
                                        div()
                                            .size(px(8.))
                                            .rounded_full()
                                            .when(by_chain, |d| {
                                                d.bg(swatch(structure.chain_color(i, palette)))
                                            })
                                            .when(!by_chain, |d| {
                                                d.border_1().border_color(theme::faint())
                                            }),
                                    )
                                    .tooltip(|_, cx| {
                                        cx.new(|_| Tip("Choose a colour".into())).into()
                                    })
                                    .on_click(cx.listener(move |this, _, _, cx| {
                                        this.picking = Some(Pick::Chain(i));
                                        cx.notify();
                                    })),
                            )
                            .child(
                                div()
                                    .flex_1()
                                    .text_size(px(11.5))
                                    .text_color(theme::text())
                                    .child(name.clone()),
                            )
                            .child(
                                div()
                                    .font_family(theme::MONO)
                                    .text_size(px(10.))
                                    .text_color(theme::faint())
                                    .child({
                                        let own = chain
                                            .residues
                                            .iter()
                                            .map(|r| &structure.residues[*r as usize])
                                            .filter(|r| !r.outside_polymer());
                                        let (n, nucleic) = own.fold((0, 0), |(n, k), r| {
                                            (n + 1, k + r.nucleic.is_some() as usize)
                                        });
                                        // Nucleotides, where most of the chain is them.
                                        format!("{n} {}", if nucleic * 2 > n { "nt" } else { "aa" })
                                    }),
                            )
                            .when(chain.color.is_some(), |d| {
                                d.child(self.unpick(i, Pick::Chain(i), cx))
                            })
                            .when(picking, |d| {
                                d.child(
                                    deferred(self.color_picker(
                                        Pick::Chain(i),
                                        name,
                                        default,
                                        chain.color,
                                        cx,
                                    ))
                                    .with_priority(3),
                                )
                            })
                    }),
            )
            .when(rest > 0 || changed, |d| {
                d.child(
                    div()
                        .pl(px(21.))
                        .flex()
                        .items_center()
                        .text_size(px(10.5))
                        .child(div().flex_1().when(rest > 0, |d| {
                            d.child(
                                div()
                                    .id("more-chains")
                                    .debug_selector(|| "more-chains".into())
                                    .text_color(theme::faint())
                                    .cursor_pointer()
                                    .hover(|s| s.text_color(theme::ink()))
                                    .child(format!("and {rest} more"))
                                    .on_click(cx.listener(|this, _, _, cx| {
                                        this.all_chains = true;
                                        cx.notify();
                                    })),
                            )
                        }))
                        .when(changed, |d| {
                            d.child(
                                div()
                                    .id("reset-chains")
                                    .debug_selector(|| "reset-chains".into())
                                    .text_size(px(9.5))
                                    .text_color(theme::muted())
                                    .cursor_pointer()
                                    .hover(|s| s.text_color(theme::ink()))
                                    .child("Reset")
                                    .on_click(cx.listener(|this, _, _, cx| {
                                        this.reset_chain_colors(cx);
                                    })),
                            )
                        }),
                )
            })
    }

    /// One row per kind of ligand beside a polymer: its colour, which opens
    /// a picker, how many there are, and whether the colour goes on its
    /// carbons, on every atom, or nowhere, leaving its elements' colours.
    fn ligand_list(&self, structure: &Structure, cx: &mut Context<Self>) -> impl IntoElement {
        let ligands = structure.ligands();
        let carbon = render::element_colors(structure, self.palette, self.standard_colors)[6];
        // An ion or a phosphate has no carbons to colour, so is offered only
        // its elements' colours or a solid one.
        let carbonless: Vec<bool> = ligands
            .iter()
            .map(|(name, _)| {
                !structure.atoms.iter().any(|a| {
                    a.element.0 == 6 && structure.residues[a.residue as usize].name == *name
                })
            })
            .collect();
        div()
            .when(!ligands.is_empty(), |d| d.pt(px(14.)))
            .flex()
            .flex_col()
            .gap(px(5.))
            .when(!ligands.is_empty(), |d| d.child(key("ligands")))
            .children(ligands.into_iter().enumerate().map(|(i, (name, n))| {
                let paint = structure.ligand_paint.get(&name).copied();
                let carbonless = carbonless[i];
                let picking = self.picking == Some(Pick::Ligand(name.clone()));
                let mode =
                    |k: usize, label: &'static str, chosen: bool, to: Option<LigandPaint>| {
                        let name = name.clone();
                        div()
                            .id(("ligand-mode", i * 3 + k))
                            .debug_selector({
                                let name = name.clone();
                                move || format!("ligand-{name}-{label}")
                            })
                            .px(px(7.))
                            .py(px(1.))
                            .rounded(px(5.))
                            .text_size(px(10.))
                            .when(chosen, |d| {
                                d.bg(theme::surface()).text_color(theme::ink()).shadow_sm()
                            })
                            .when(!chosen, |d| {
                                d.text_color(theme::muted())
                                    .cursor_pointer()
                                    .hover(|s| s.text_color(theme::ink()))
                            })
                            .child(label)
                            .on_click(cx.listener(move |this, _, _, cx| {
                                this.set_ligand_paint(name.clone(), to, cx)
                            }))
                    };
                let color = paint.map(|p| p.color()).unwrap_or(LIGAND_GOLD);
                let pick = Pick::Ligand(name.clone());
                div()
                    .relative()
                    .flex()
                    .items_center()
                    .gap(px(8.))
                    .child(
                        div()
                            .id(("ligand-dot", i))
                            .debug_selector({
                                let name = name.clone();
                                move || format!("ligand-dot-{name}")
                            })
                            .flex_none()
                            .size(px(13.))
                            .flex()
                            .items_center()
                            .justify_center()
                            .rounded_full()
                            .cursor_pointer()
                            .when(paint.is_some() && !picking, |d| {
                                d.border_1().border_color(theme::faint())
                            })
                            .when(picking, |d| d.border_1().border_color(theme::ink()))
                            .hover(|d| d.bg(theme::sunken()))
                            .child(
                                div()
                                    .size(px(8.))
                                    .rounded_full()
                                    .bg(swatch(paint.map(|p| p.color()).unwrap_or(carbon))),
                            )
                            .tooltip(|_, cx| cx.new(|_| Tip("Choose a colour".into())).into())
                            .on_click(cx.listener(move |this, _, _, cx| {
                                this.picking = Some(pick.clone());
                                cx.notify();
                            })),
                    )
                    .child(
                        div()
                            .flex_1()
                            .flex()
                            .items_baseline()
                            .gap(px(6.))
                            .child(
                                div()
                                    .text_size(px(11.5))
                                    .text_color(theme::text())
                                    .child(name.clone()),
                            )
                            .when(n > 1, |d| {
                                d.child(
                                    div()
                                        .font_family(theme::MONO)
                                        .text_size(px(10.))
                                        .text_color(theme::faint())
                                        .child(format!("×{n}")),
                                )
                            })
                            .child({
                                let name = name.clone();
                                div()
                                    .id(("ligand-site", i))
                                    .debug_selector({
                                        let name = name.clone();
                                        move || format!("ligand-site-{name}")
                                    })
                                    .px(px(5.))
                                    .rounded(px(4.))
                                    .text_size(px(10.))
                                    .text_color(theme::faint())
                                    .cursor_pointer()
                                    .hover(|d| d.bg(theme::sunken()).text_color(theme::ink()))
                                    .child("Site")
                                    .tooltip(|_, cx| {
                                        cx.new(|_| Tip("Its binding site".into())).into()
                                    })
                                    .on_click(
                                        cx.listener(move |this, _, _, cx| {
                                            this.next_site(&name, cx)
                                        }),
                                    )
                            }),
                    )
                    .child(
                        div()
                            .flex()
                            .p(px(2.))
                            .gap(px(2.))
                            .rounded(px(7.))
                            .bg(theme::sunken())
                            .child(mode(0, "Element", paint.is_none(), None))
                            .when(!carbonless, |d| {
                                d.child(mode(
                                    1,
                                    "Carbons",
                                    matches!(paint, Some(LigandPaint::Carbons(_))),
                                    Some(LigandPaint::Carbons(color)),
                                ))
                            })
                            .child(mode(
                                2,
                                "Solid",
                                matches!(paint, Some(LigandPaint::Solid(_))),
                                Some(LigandPaint::Solid(color)),
                            )),
                    )
                    .when(paint.is_some(), |d| {
                        d.child(self.unpick(i, Pick::Ligand(name.clone()), cx))
                    })
                    .when(picking, |d| {
                        d.child(
                            deferred(self.color_picker(
                                Pick::Ligand(name.clone()),
                                name.clone(),
                                carbon,
                                paint.map(|p| p.color()),
                                cx,
                            ))
                            .with_priority(3),
                        )
                    })
            }))
    }

    /// One row per element — or per type, for atoms of no element — with
    /// its colour, its size and opacity as sliders and as figures to type
    /// over, and its count. The colour opens a picker.
    fn species_table(&self, structure: &Structure, cx: &mut Context<Self>) -> impl IntoElement {
        let by_element = self.coloring == Coloring::Element
            || (self.compare.is_none() && structure.has_ligands());
        let palette = self.palette;
        let colors = render::element_colors(structure, palette, self.standard_colors);
        let mut rows: Vec<(Species, String, [f32; 3], usize)> = structure
            .composition()
            .into_iter()
            .filter(|(e, _)| e.is_known())
            .map(|(e, n)| {
                (
                    Species::Element(e.0),
                    e.symbol().to_string(),
                    colors[e.0 as usize],
                    n,
                )
            })
            .collect();
        let mut by_kind = vec![0usize; structure.kinds.len()];
        for a in structure.atoms.iter().filter(|a| !a.element.is_known()) {
            if let Some(n) = by_kind.get_mut(a.kind as usize) {
                *n += 1;
            }
        }
        rows.extend(
            by_kind
                .iter()
                .enumerate()
                .filter(|(_, n)| **n > 0)
                .map(|(k, n)| {
                    let name = structure.kinds[k].clone();
                    (
                        Species::Kind(name.clone()),
                        format!("Type {name}"),
                        palette.chain(k),
                        *n,
                    )
                }),
        );
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
            .children(
                rows.into_iter()
                    .enumerate()
                    .map(|(i, (species, name, default, n))| {
                        let look = self.appearance.get(&species);
                        let shown = look.color.unwrap_or(default);
                        let picking =
                            self.picking.as_ref() == Some(&Pick::Species(species.clone()));
                        let pick = Pick::Species(species.clone());
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
                                            .when(look.color.is_some() && !picking, |d| {
                                                d.border_1().border_color(theme::faint())
                                            })
                                            .when(picking, |d| {
                                                d.border_1().border_color(theme::ink())
                                            })
                                            .hover(|d| d.bg(theme::sunken()))
                                            .child(
                                                div()
                                                    .size(px(8.))
                                                    .rounded_full()
                                                    .when(by_element, |d| {
                                                        d.bg(swatch(shown))
                                                            .opacity(look.alpha.max(0.2))
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
                                            .flex_1()
                                            .overflow_hidden()
                                            .text_size(px(11.5))
                                            .text_color(theme::text())
                                            .whitespace_nowrap()
                                            .child(name.clone()),
                                    ),
                            )
                            .when(look.color.is_some(), |d| {
                                d.child(self.unpick(i, Pick::Species(species.clone()), cx))
                            })
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
                                    deferred(self.color_picker(
                                        Pick::Species(species.clone()),
                                        name,
                                        default,
                                        look.color,
                                        cx,
                                    ))
                                    .with_priority(3),
                                )
                            })
                    }),
            )
            .child(self.standard_colors_toggle(cx))
            .children(self.visibility_toggles(structure, cx))
    }

    /// Water and hydrogens on or off at a click, where the structure has them.
    fn visibility_toggles(
        &self,
        structure: &Structure,
        cx: &mut Context<Self>,
    ) -> Option<gpui::AnyElement> {
        let water = structure
            .residues
            .iter()
            .any(|r| crate::protein::is_solvent(&r.name));
        let hydrogens = structure
            .atoms
            .iter()
            .any(|a| a.element == crate::protein::Element::H);
        if !water && !hydrogens {
            return None;
        }
        let water_shown = !self.appearance.hide_solvent;
        let h_shown = !self.appearance.get(&Species::Element(1)).hidden();
        let toggle = |id: &'static str, label: &'static str, shown: bool| {
            div()
                .id(id)
                .debug_selector(move || id.into())
                .flex()
                .items_center()
                .gap(px(6.))
                .px(px(8.))
                .py(px(2.))
                .rounded(px(5.))
                .text_size(px(10.5))
                .cursor_pointer()
                .when(shown, |d| {
                    d.bg(theme::surface()).text_color(theme::ink()).shadow_sm()
                })
                .when(!shown, |d| {
                    d.text_color(theme::muted())
                        .hover(|s| s.text_color(theme::ink()))
                })
                .child(
                    div()
                        .size(px(6.))
                        .rounded_full()
                        .when(shown, |d| d.bg(theme::ink()))
                        .when(!shown, |d| d.border_1().border_color(theme::faint())),
                )
                .child(label)
        };
        Some(
            div()
                .mt(px(6.))
                .flex()
                .items_center()
                .gap(px(10.))
                .child(
                    div()
                        .text_size(px(11.))
                        .text_color(theme::muted())
                        .child("Show"),
                )
                .child(
                    div()
                        .flex()
                        .p(px(2.))
                        .gap(px(2.))
                        .rounded(px(7.))
                        .bg(theme::sunken())
                        .when(water, |d| {
                            d.child(toggle("show-water", "Water", water_shown).on_click(
                                cx.listener(|this, _, _, cx| {
                                    this.appearance.hide_solvent = !this.appearance.hide_solvent;
                                    this.rebuild(cx);
                                }),
                            ))
                        })
                        .when(hydrogens, |d| {
                            d.child(toggle("show-hydrogens", "Hydrogens", h_shown).on_click(
                                cx.listener(|this, _, _, cx| {
                                    let h = Species::Element(1);
                                    let look = this.appearance.get(&h);
                                    let look = if look.hidden() {
                                        Look {
                                            size: 1.0,
                                            alpha: 1.0,
                                            ..look
                                        }
                                    } else {
                                        Look { size: 0.0, ..look }
                                    };
                                    this.appearance.set(h, look);
                                    this.rebuild(cx);
                                }),
                            ))
                        }),
                )
                .into_any_element(),
        )
    }

    /// Whether elements wear the standard colours or the palette's. Only
    /// elements: chains and ramps keep the palette either way.
    fn standard_colors_toggle(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let standard = self.standard_colors;
        let option = |id: &'static str, label: &'static str, chosen: bool, set: bool| {
            div()
                .id(id)
                .debug_selector(move || id.into())
                .px(px(8.))
                .py(px(2.))
                .rounded(px(5.))
                .text_size(px(10.5))
                .when(chosen, |d| {
                    d.bg(theme::surface()).text_color(theme::ink()).shadow_sm()
                })
                .when(!chosen, |d| {
                    d.text_color(theme::muted())
                        .cursor_pointer()
                        .hover(|s| s.text_color(theme::ink()))
                })
                .child(label)
                .on_click(cx.listener(move |this, _, _, cx| {
                    if this.standard_colors != set {
                        this.standard_colors = set;
                        this.rebuild(cx);
                    }
                }))
        };
        div()
            .mt(px(4.))
            .flex()
            .items_center()
            .gap(px(10.))
            .child(
                div()
                    .text_size(px(11.))
                    .text_color(theme::muted())
                    .child("Element colours"),
            )
            .child(
                div()
                    .flex()
                    .p(px(2.))
                    .gap(px(2.))
                    .rounded(px(7.))
                    .bg(theme::sunken())
                    .child(
                        option("element-colours-standard", "Standard", standard, true).tooltip(
                            |_, cx| {
                                cx.new(|_| Tip("Carbon grey, oxygen red, nitrogen blue".into()))
                                    .into()
                            },
                        ),
                    )
                    .child(
                        option("element-colours-palette", "Palette", !standard, false).tooltip(
                            |_, cx| {
                                cx.new(|_| Tip("From the palette, most common first".into()))
                                    .into()
                            },
                        ),
                    ),
            )
    }

    /// Presets, then the colour as hex and as RGB to type over. The first
    /// preset is the species' or chain's own colour, which clears any choice.
    fn color_picker(
        &self,
        target: Pick,
        name: String,
        default: [f32; 3],
        chosen: Option<[f32; 3]>,
        cx: &mut Context<Self>,
    ) -> impl IntoElement {
        let current = chosen.unwrap_or(default);
        let swatch_button =
            |id: usize, c: [f32; 3], chosen: bool, set: Option<[f32; 3]>, target: Pick| {
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
                    .on_click(
                        cx.listener(move |this, _, _, cx| this.paint(target.clone(), set, cx)),
                    )
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
                        swatch_button(0, default, chosen.is_none(), None, target.clone())
                            .tooltip(|_, cx| cx.new(|_| Tip("Custom colour".into())).into()),
                    )
                    .children(SWATCHES.iter().enumerate().map(|(k, c)| {
                        swatch_button(k + 1, *c, chosen == Some(*c), Some(*c), target.clone())
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
                        Field::Hex(target.clone()),
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
                        Field::Rgb(target),
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
            .id((
                if knob == Knob::Size {
                    "size"
                } else {
                    "opacity"
                },
                row,
            ))
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
                    .child(
                        div()
                            .h_full()
                            .w(relative(t))
                            .rounded_full()
                            .bg(theme::muted()),
                    )
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

    fn info_atoms(&self, structure: &Structure, cx: &mut Context<Self>) -> Vec<gpui::AnyElement> {
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

        let head = div()
            .flex()
            .flex_col()
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
            .children(self.colour_key(structure, cx))
            .children(self.volume_tools(structure, cx));
        let atoms = self.species_table(structure, cx);
        let facts = div()
            .flex()
            .flex_col()
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
                    .children(self.periodic_tools(structure, cx))
                    .when(structure.frame_count() > 1, |d| {
                        d.child(spec("frames", thousands(structure.frame_count()), false))
                    })
                    .children(self.trajectory_analysis(structure, cx)),
            )
            .children(self.source_line(structure));
        vec![
            head.into_any_element(),
            atoms.into_any_element(),
            facts.into_any_element(),
        ]
    }

    /// What is on screen, set like the opening of an article: the category in
    /// small capitals, the name large, the sentence beneath it, then the
    /// figures in a ruled table that reads straight down the left edge.
    /// The sequence, chain by chain, in numbered blocks of ten. A letter
    /// glows its residue in the viewport while pointed at and picks it out
    /// when clicked; what is picked or pointed at in the viewport is marked
    /// here in turn.
    fn sequence_view(&self, structure: &Structure, cx: &mut Context<Self>) -> impl IntoElement {
        let selected = self.selected_residue();
        let pointed = self.hovered_residue().or(self.seq_hover);
        let total = structure.chains.len();
        let shown = if self.all_chains { total } else { total.min(6) };
        let several = total > 1;
        div().pt(px(9.)).flex().flex_col().gap(px(8.)).children(
            structure
                .chains
                .iter()
                .take(shown)
                .enumerate()
                .filter_map(|(ci, chain)| {
                    let residues: Vec<usize> = chain
                        .residues
                        .iter()
                        .map(|&r| r as usize)
                        .filter(|&r| !structure.residues[r].outside_polymer())
                        .collect();
                    if residues.is_empty() {
                        return None;
                    }
                    let blocks = residues.chunks(SEQUENCE_BLOCK).map(|block| {
                        let first = structure.residues[block[0]].seq;
                        div()
                            .flex()
                            .flex_col()
                            .child(
                                div()
                                    .h(px(11.))
                                    .font_family(theme::MONO)
                                    .text_size(px(8.))
                                    .text_color(theme::faint())
                                    .child(first.to_string()),
                            )
                            .child(div().flex().children(block.iter().map(|&ri| {
                                let r = &structure.residues[ri];
                                let tip: SharedString = format!("{} {}", r.name, r.seq).into();
                                let on = selected == Some(ri);
                                let near = pointed == Some(ri);
                                div()
                                    .id(("residue", ri))
                                    .debug_selector(move || format!("residue-{ri}"))
                                    .w(px(7.5))
                                    .h(px(14.))
                                    .flex()
                                    .items_center()
                                    .justify_center()
                                    .rounded(px(2.))
                                    .cursor_pointer()
                                    .font_family(theme::MONO)
                                    .text_size(px(10.5))
                                    .when(on, |d| d.bg(theme::ink()).text_color(theme::on_ink()))
                                    .when(!on && near, |d| {
                                        d.bg(theme::sunken()).text_color(theme::ink())
                                    })
                                    .when(!on && !near, |d| d.text_color(theme::muted()))
                                    .child(residue_letter(r).to_string())
                                    .tooltip(move |_, cx| cx.new(|_| Tip(tip.clone())).into())
                                    .on_hover(cx.listener(move |this, hovered: &bool, _, cx| {
                                        if *hovered {
                                            this.set_seq_hover(Some(ri), cx);
                                        } else if this.seq_hover == Some(ri) {
                                            this.set_seq_hover(None, cx);
                                        }
                                    }))
                                    .on_click(cx.listener(move |this, _, _, cx| {
                                        let Some(s) = this.structure.as_ref() else {
                                            return;
                                        };
                                        let target = residue_pick(s, ri)
                                            .map(|picked| Target { which: 0, picked });
                                        this.select(target, cx);
                                    }))
                            })))
                    });
                    Some(
                        div()
                            .flex()
                            .flex_col()
                            .when(several, |d| {
                                let id = chain.id.trim();
                                d.child(div().text_size(px(9.)).text_color(theme::faint()).child(
                                    spaced(&if id.is_empty() {
                                        format!("CHAIN {}", ci + 1)
                                    } else {
                                        format!("CHAIN {}", id.to_uppercase())
                                    }),
                                ))
                            })
                            .child(
                                div()
                                    .flex()
                                    .flex_wrap()
                                    .gap_x(px(6.))
                                    .gap_y(px(2.))
                                    .children(blocks),
                            ),
                    )
                }),
        )
    }

    /// What is known of the structure, in sections the panel stacks in one
    /// column or the vertical layout's band sets side by side.
    fn info(&self, structure: &Structure, cx: &mut Context<Self>) -> Vec<gpui::AnyElement> {
        if !structure.has_polymer() {
            return self.info_atoms(structure, cx);
        }
        let sequence = structure.sequence();
        let (name, blurb) = headline(structure);
        // The numbering is worth showing only when the file does not start at
        // one, which is the sign of a fragment or a trimmed construct.
        // Water is no residue of the protein's.
        let own: Vec<&crate::protein::Residue> = structure
            .residues
            .iter()
            .filter(|r| !crate::protein::is_solvent(&r.name))
            .collect();
        let residues = match (own.first(), own.last()) {
            (Some(a), Some(b)) if a.seq != 1 => {
                format!("{} aa  ·  {}–{}", thousands(own.len()), a.seq, b.seq)
            }
            _ => format!("{} aa", thousands(own.len())),
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

        let head = div()
            .flex()
            .flex_col()
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
            .children(self.assembly_choice(structure, cx))
            .child(self.legend(structure, cx))
            .children(self.colour_key(structure, cx))
            .children(self.volume_tools(structure, cx))
            .child(self.ligand_list(structure, cx))
            .child(self.layer_list(cx));
        let atoms = self.species_table(structure, cx);
        let facts = div()
            .flex()
            .flex_col()
            .child(
                div()
                    .pt(px(16.))
                    .flex()
                    .flex_col()
                    .child(spec("structure", structure.id.to_uppercase(), true))
                    .child(spec("residues", residues, false))
                    .child(spec("chains", chains, false))
                    .child(spec("atoms", thousands(structure.atoms.len()), false))
                    .child(spec(
                        "mass",
                        format!("{:.1} kDa", structure.mass_kda()),
                        false,
                    ))
                    // A deposited crystal's cell says nothing a reader of a
                    // protein needs; a simulation box does.
                    .when_some(structure.cell.filter(|_| periodic(structure)), |d, cell| {
                        let [a, b, c] = cell.lengths();
                        d.child(spec("cell", format!("{a:.2} × {b:.2} × {c:.2} Å"), false))
                    })
                    .children(self.periodic_tools(structure, cx))
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
                                        d.text_color(theme::muted()).hover(|s| {
                                            s.bg(theme::sunken()).text_color(theme::text())
                                        })
                                    })
                                    .child(if self.copied { "Copied" } else { "Copy" })
                                    .on_click(cx.listener(move |this, _, _, cx| {
                                        this.copy_sequence(sequence.clone(), cx)
                                    })),
                            ),
                    ),
            )
            .child(self.sequence_view(structure, cx));
        // The plots, each closed under a row of its own like the facts.
        let plots = div()
            .flex()
            .flex_col()
            .child(
                div()
                    .pt(px(8.))
                    .flex()
                    .flex_col()
                    .children(self.ramachandran(structure, cx))
                    .children(self.pae_plot(structure, cx)),
            )
            .children(self.trajectory_analysis(structure, cx))
            .when(structure.frame_count() > 1, |d| {
                d.child(
                    div()
                        .pt(px(8.))
                        .text_size(px(10.5))
                        .text_color(theme::faint())
                        .child(format!("{} frames", thousands(structure.frame_count()))),
                )
            })
            .children(self.source_line(structure));
        vec![
            head.into_any_element(),
            atoms.into_any_element(),
            facts.into_any_element(),
            plots.into_any_element(),
        ]
    }

    /// The search, the status and what is known of the structure: a column
    /// down the left, or, in the vertical layout, a band under the viewport
    /// with the search at its left and the sections side by side beside it,
    /// each as wide as the column would be. The masthead is then over the
    /// viewport rather than in the band.
    fn panel(
        &self,
        folding: bool,
        focused: bool,
        structure: Option<Arc<Structure>>,
        vertical: bool,
        cx: &mut Context<Self>,
    ) -> impl IntoElement {
        let sections: Vec<gpui::AnyElement> = match &structure {
            Some(s) => self.info(s, cx),
            None => vec![self.examples(cx).into_any_element()],
        };
        // Beside the view, or under it in a window too narrow for columns:
        // one column, under the masthead beside the view. Under it, the
        // masthead is over the viewport instead, at the top of the window.
        if !vertical || self.narrow {
            return div()
                .id("side")
                .flex()
                .flex_col()
                .flex_none()
                .gap(px(16.))
                .pt(px(if vertical { 14. } else { 20. }))
                .map(|d| {
                    if vertical {
                        d.w_full().h(px(theme::band_height()))
                    } else {
                        d.w(px(theme::side_width())).h_full()
                    }
                })
                .bg(theme::surface())
                .when(!vertical, |d| d.child(self.masthead()))
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
                                .map(|d| match structure {
                                    Some(_) => d.child(
                                        div()
                                            .flex()
                                            .flex_col()
                                            .px(px(20.))
                                            .pb(px(20.))
                                            .child(rule())
                                            .children(sections),
                                    ),
                                    None => d.children(sections),
                                }),
                        ),
                );
        }
        let has_structure = structure.is_some();
        div()
            .id("side")
            .flex()
            .flex_none()
            .w_full()
            .h(px(theme::band_height()))
            .bg(theme::surface())
            .child(
                div()
                    .id("band-search")
                    .flex()
                    .flex_col()
                    .flex_none()
                    .gap(px(14.))
                    .w(px(PANEL_WIDTH))
                    .h_full()
                    .pt(px(18.))
                    .pb(px(16.))
                    .border_r_1()
                    .border_color(theme::border())
                    .overflow_y_scroll()
                    .child(self.search(focused, cx))
                    .when(folding, |d| d.children(self.plugin_options(cx)))
                    .child(self.status(cx)),
            )
            .child(
                div()
                    .id("panel")
                    .flex_1()
                    .min_w_0()
                    .h_full()
                    .overflow_y_scroll()
                    .child(div().flex().flex_wrap().pb(px(20.)).children(
                        sections.into_iter().map(|section| {
                            div()
                                .flex_grow()
                                .flex_basis(px(PANEL_WIDTH))
                                .max_w(px(PANEL_WIDTH + 120.))
                                .min_w(px(PANEL_WIDTH))
                                .when(has_structure, |d| d.px(px(20.)).pt(px(4.)))
                                .when(!has_structure, |d| d.pt(px(18.)))
                                .child(section)
                        }),
                    )),
            )
    }

    /// The line between the panel and the viewport, which can be dragged to
    /// give one more room and the other less. A double-click puts it back.
    fn divider(&self, vertical: bool, cx: &mut Context<Self>) -> impl IntoElement {
        let grip = div()
            .id("divider")
            .absolute()
            .when(vertical, |d| {
                d.left_0()
                    .right_0()
                    .top(px(-3.))
                    .h(px(7.))
                    .cursor_row_resize()
            })
            .when(!vertical, |d| {
                d.top_0()
                    .bottom_0()
                    .left(px(-3.))
                    .w(px(7.))
                    .cursor_col_resize()
            })
            .on_mouse_down(
                MouseButton::Left,
                cx.listener(|this, e: &MouseDownEvent, _, cx| {
                    cx.stop_propagation();
                    if e.click_count >= 2 {
                        theme::set_side_width(theme::SIDE_WIDTH);
                        theme::set_band_height(theme::BAND_HEIGHT);
                        this.persist();
                        cx.notify();
                        return;
                    }
                    this.resizing = true;
                }),
            );
        div()
            .relative()
            .flex_none()
            .bg(theme::border())
            .when(vertical, |d| d.w_full().h(px(1.)))
            .when(!vertical, |d| d.h_full().w(px(1.)))
            .when(self.resizing, |d| d.bg(theme::faint()))
            .child(grip)
    }

    fn stage(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let has_structure = self.structure.is_some();
        let busy: SharedString = self.busy.clone().unwrap_or_default();

        div()
            .id("stage")
            .relative()
            .flex_1()
            .min_w_0()
            .min_h_0()
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
                self.traced,
            ))
            .children(self.comparison(cx))
            .when(self.busy.is_none(), |d| {
                let cards: Vec<gpui::AnyElement> = [
                    self.measure_card(cx),
                    self.site_card(cx),
                    self.selected_card(cx),
                ]
                .into_iter()
                .flatten()
                .collect();
                d.children(self.site_labels())
                    .children(self.atom_labels())
                    .children(self.slab_gauge())
                    .children(self.measure_lines())
                    .children(self.measure_labels())
                    .when(!cards.is_empty(), |d| {
                        d.child(
                            div()
                                .absolute()
                                .top(px(56.))
                                .right(px(20.))
                                .flex()
                                .flex_col()
                                .items_end()
                                .gap(px(8.))
                                // Clicks on the cards are theirs, not the viewport's.
                                .on_mouse_down(MouseButton::Left, |_, _, cx| cx.stop_propagation())
                                .on_mouse_move(cx.listener(Self::over_card))
                                .children(cards),
                        )
                    })
                    .children(self.hover_tip())
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
                        .children(self.empty_message().into_iter().enumerate().map(|(i, line)| {
                            div()
                                .text_size(px(if i == 0 { 12. } else { 10.5 }))
                                .text_color(theme::faint())
                                .child(line)
                        })),
                )
            })
            .when(self.busy.is_some(), |d| d.child(self.working(busy.clone())))
            .children(self.tracing_label())
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
            .when(has_structure && self.busy.is_none() && self.frame_count() <= 1 && !self.narrow, |d| {
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
                        .child(format!(
                            "drag to orbit · shift-drag to pan · {}-drag to roll · scroll to zoom · double-click to reset · {} for shortcuts",
                            theme::keys("⌘"),
                            theme::keys("⌘/")
                        )),
                )
            })
            .children(self.console_view())
            .children(self.movie_dialog(cx))
            .children(self.shortcut_sheet(cx))
    }
}

impl Viewer {
    /// What is said about an atom or residue: a coloured badge and its
    /// name, over a few labelled facts.
    fn atom_card(&self, t: Target) -> Option<gpui::Div> {
        let s = self.shown(t.which)?;
        let card = s.card(t.picked.atom, t.picked.residue)?;
        // Which of a superimposed pair, when both are drawn.
        let whose = (self.compare.is_some() && self.visible().len() > 1).then(|| {
            let s = if t.which == 0 {
                self.structure.as_ref()
            } else {
                self.compare.as_ref()
            };
            let tint = if t.which == 0 {
                align::REFERENCE_TINT
            } else {
                align::PREDICTION_TINT
            };
            (s.map(|s| s.id.to_uppercase()).unwrap_or_default(), tint)
        });
        let on_ink = theme::on_ink();
        let soft = |a: f32| gpui::Hsla { a, ..on_ink };
        // A residue shows its letter in its kind's colour; an atom, its
        // element or type in the colour the atoms list gives it.
        let atom = &s.atoms[t.picked.atom];
        let species_color = |species: &Species, fallback: [f32; 3]| {
            self.appearance.get(species).color.unwrap_or(fallback)
        };
        let element_color = |c: [f32; 3]| {
            let z = atom.element.0;
            let c = if self.standard_colors {
                c
            } else {
                render::element_colors(s, self.palette, false)[z as usize]
            };
            species_color(&Species::Element(z), c)
        };
        let (badge, color, round) = match (card.class, card.element, &card.kind) {
            (_, Some((e, c)), _) => (e.to_string(), element_color(c), true),
            (_, None, Some((k, name))) => (
                name.clone(),
                species_color(&Species::Kind(name.clone()), self.palette.chain(*k)),
                true,
            ),
            // A ligand has no letter; its code stands in.
            (Some(('X', _, c)), None, None) => (
                card.tag.split(' ').next().unwrap_or("").to_string(),
                c,
                false,
            ),
            (Some((code, _, c)), None, None) => (code.to_string(), c, false),
            (None, None, None) => ("·".to_string(), [0.55, 0.55, 0.58], true),
        };
        // Three letters is all that fits.
        let badge: String = badge.chars().take(3).collect();
        let light = 0.299 * color[0] + 0.587 * color[1] + 0.114 * color[2] > 0.62;
        let badge = div()
            .flex_none()
            .size(px(26.))
            .flex()
            .items_center()
            .justify_center()
            .when(round, |d| d.rounded_full())
            .when(!round, |d| d.rounded(px(6.)))
            .bg(swatch(color))
            // A pale badge, hydrogen's, would melt into a pale card.
            .border_1()
            .border_color(gpui::Hsla {
                a: 0.12,
                ..gpui::black()
            })
            .text_size(px(match badge.chars().count() {
                0 | 1 => 13.,
                2 => 11.,
                _ => 9.5,
            }))
            .font_weight(gpui::FontWeight::BOLD)
            .text_color(if light {
                gpui::rgb(0x1d1d22)
            } else {
                gpui::rgb(0xffffff)
            })
            .child(badge);
        let mut sub = card.tag.clone();
        if let Some(c) = &card.chain {
            sub.push_str(&format!(" · chain {c}"));
        }
        let header = div().flex().items_center().gap(px(8.)).child(badge).child(
            div()
                .flex()
                .flex_col()
                .child(
                    div()
                        .text_size(px(12.5))
                        .font_weight(gpui::FontWeight::SEMIBOLD)
                        .child(card.title.clone()),
                )
                .child(
                    div()
                        .font_family(theme::MONO)
                        .text_size(px(9.5))
                        .text_color(soft(0.6))
                        .child(sub),
                ),
        );

        let dot = |c: [f32; 3]| div().flex_none().size(px(7.)).rounded_full().bg(swatch(c));
        let row = |label: &'static str, value: gpui::AnyElement| {
            div()
                .flex()
                .items_center()
                .gap(px(8.))
                .child(
                    div()
                        .w(px(64.))
                        .flex_none()
                        .text_size(px(9.))
                        .text_color(soft(0.5))
                        .child(if label == "pLDDT" {
                            label.to_string()
                        } else {
                            spaced(&label.to_uppercase())
                        }),
                )
                .child(value)
        };
        let mut rows: Vec<gpui::AnyElement> = Vec::new();
        if let Some(name) = &card.atom {
            let mut v = div().flex().items_center().gap(px(5.)).child(
                div()
                    .font_family(theme::MONO)
                    .text_size(px(10.5))
                    .child(name.clone()),
            );
            if let Some((e, c)) = card
                .element
                .filter(|(e, _)| e != name)
                .map(|(e, c)| (e, element_color(c)))
            {
                v = v.child(
                    div()
                        .flex()
                        .items_center()
                        .gap(px(4.))
                        .px(px(5.))
                        .rounded(px(4.))
                        .bg(soft(0.1))
                        .child(dot(c))
                        .child(div().text_size(px(9.5)).text_color(soft(0.75)).child(e)),
                );
            }
            rows.push(row("atom", v.into_any_element()).into_any_element());
        }
        if let Some((_, kind, c)) = card.class {
            let v = div()
                .flex()
                .items_center()
                .gap(px(5.))
                .child(dot(c))
                .child(div().text_size(px(10.5)).child(kind));
            rows.push(row("kind", v.into_any_element()).into_any_element());
        }
        if let Some(ss) = card.ss {
            let (_, name, c) = crate::render::SS_COLORS
                .iter()
                .find(|(s, _, _)| *s == ss)
                .copied()
                .unwrap_or((ss, "Coil", [0.78; 3]));
            let v = div()
                .flex()
                .items_center()
                .gap(px(5.))
                .child(dot(c))
                .child(div().text_size(px(10.5)).child(name.to_lowercase()));
            rows.push(row("fold", v.into_any_element()).into_any_element());
        }
        if let Some((label, b)) = card.b {
            let mut v = div().flex().items_center().gap(px(6.)).child(
                div()
                    .font_family(theme::MONO)
                    .text_size(px(10.5))
                    .child(format!("{b:.1}")),
            );
            if label == "pLDDT" {
                let c = crate::render::plddt_band_color(b);
                // AlphaFold's cut-offs, at 90, 70 and 50.
                let band = match b {
                    b if b >= 90.0 => "very high",
                    b if b >= 70.0 => "confident",
                    b if b >= 50.0 => "low",
                    _ => "very low",
                };
                v = v
                    .child(
                        div()
                            .w(px(40.))
                            .h(px(4.))
                            .rounded_full()
                            .bg(soft(0.15))
                            .child(
                                div()
                                    .h_full()
                                    .rounded_full()
                                    .w(px(40. * (b / 100.).clamp(0.04, 1.)))
                                    .bg(swatch(c)),
                            ),
                    )
                    .child(div().text_size(px(9.5)).text_color(soft(0.6)).child(band));
            }
            rows.push(row(label, v.into_any_element()).into_any_element());
        }
        if let Some([x, y, z]) = card.pos {
            let v = div()
                .font_family(theme::MONO)
                .text_size(px(10.))
                .text_color(soft(0.75))
                .child(format!("{x:.2}  {y:.2}  {z:.2} Å"));
            rows.push(row("xyz", v.into_any_element()).into_any_element());
        }

        Some(
            div()
                .flex()
                .flex_col()
                .gap(px(7.))
                .min_w(px(150.))
                .px(px(10.))
                .py(px(8.))
                .rounded(px(9.))
                .bg(theme::ink())
                .text_color(on_ink)
                .shadow_lg()
                .whitespace_nowrap()
                // Which of the pair, in its colour from the pills above.
                .children(whose.map(|(w, tint)| {
                    div()
                        .flex()
                        .items_center()
                        .gap(px(5.))
                        .child(div().size(px(6.)).rounded_full().bg(swatch(tint)))
                        .child(
                            div()
                                .text_size(px(9.))
                                .text_color(soft(0.6))
                                .child(spaced(&w)),
                        )
                }))
                .child(header)
                .when(!rows.is_empty(), |d| {
                    d.child(div().h(px(1.)).w_full().bg(soft(0.12)))
                        .child(div().flex().flex_col().gap(px(4.)).children(rows))
                }),
        )
    }

    /// The card beside the pointer, naming what it is over.
    fn hover_tip(&self) -> Option<gpui::AnyElement> {
        let (t, at) = self.hover?;
        // Nothing to add over the thing already named in the corner.
        if self
            .selected
            .is_some_and(|s| (s.which, s.picked.atom) == (t.which, t.picked.atom))
        {
            return None;
        }
        let size = self.cache.borrow().bounds.size;
        let card = self.atom_card(t)?;
        // Kept on the near side of the pointer to the middle, so it does not
        // run off the edge.
        let right = f32::from(at.x) > f32::from(size.width) * 0.6;
        let below = f32::from(at.y) < f32::from(size.height) * 0.7;
        let mut tip = div().absolute();
        tip = if right {
            tip.right(size.width - at.x + px(14.))
        } else {
            tip.left(at.x + px(14.))
        };
        tip = if below {
            tip.top(at.y + px(14.))
        } else {
            tip.bottom(size.height - at.y + px(14.))
        };
        Some(tip.child(card).into_any_element())
    }

    /// The lines of every measurement, and of the one under way, as drawn
    /// in the viewport from its top left corner.
    fn measure_segments(&self) -> Vec<(Point<Pixels>, Point<Pixels>)> {
        let partial = Measurement {
            kind: self.measure_kind,
            which: self.pending_which,
            atoms: self.pending.clone(),
        };
        self.measurements
            .iter()
            .chain(std::iter::once(&partial))
            .filter_map(|m| Some((m, self.shown(m.which)?)))
            .flat_map(|(m, s)| m.segments(s))
            .filter_map(|(a, b)| Some((self.project(a)?, self.project(b)?)))
            .collect()
    }

    /// Measurements drawn on top of the molecule rather than in it, so a line
    /// through the middle of a protein is not lost inside it: dashes, dark
    /// under yellow, to stand out on either page.
    fn measure_lines(&self) -> Option<gpui::AnyElement> {
        // Measurements in yellow, and the binding site's contacts each in
        // their kind's colour.
        let mut segments: Vec<(Point<Pixels>, Point<Pixels>, u32)> = self
            .measure_segments()
            .into_iter()
            .map(|(a, b)| (a, b, 0xffd23f))
            .collect();
        segments.extend(self.contact_segments());
        if segments.is_empty() {
            return None;
        }
        Some(
            canvas(
                |_, _, _| {},
                move |bounds, _, window, _| {
                    let dash = [px(5.), px(4.)];
                    let mut colors: Vec<u32> = segments.iter().map(|s| s.2).collect();
                    colors.sort_unstable();
                    colors.dedup();
                    // All the dark undersides first, so no colour is
                    // crossed out by another line's.
                    let stroke =
                        |width: f32, color: gpui::Hsla, only: Option<u32>, window: &mut Window| {
                            let mut path = gpui::PathBuilder::stroke(px(width)).dash_array(&dash);
                            for &(a, b, c) in &segments {
                                if only.is_none_or(|o| o == c) {
                                    path.move_to(bounds.origin + a);
                                    path.line_to(bounds.origin + b);
                                }
                            }
                            if let Ok(path) = path.build() {
                                window.paint_path(path, color);
                            }
                        };
                    stroke(3.5, gpui::Hsla::from(gpui::rgba(0x1a1a1acc)), None, window);
                    for c in colors {
                        stroke(1.75, gpui::Hsla::from(gpui::rgb(c)), Some(c), window);
                    }
                },
            )
            .absolute()
            .inset_0()
            .into_any_element(),
        )
    }

    /// The binding site's contacts, as drawn in the viewport, with their
    /// colours; the kinds hidden from its card are left out.
    fn contact_segments(&self) -> Vec<(Point<Pixels>, Point<Pixels>, u32)> {
        let (Some(site), Some(s)) = (&self.site, self.shown(0)) else {
            return Vec::new();
        };
        site.contacts
            .iter()
            .filter(|c| !site.hidden.contains(&c.kind))
            .filter_map(|c| {
                let a = self.project(s.atoms[c.a].pos)?;
                let b = self.project(s.atoms[c.b].pos)?;
                Some((a, b, c.kind.color()))
            })
            .collect()
    }

    fn toggle_contacts(&mut self, kind: contacts::Kind, cx: &mut Context<Self>) {
        if let Some(site) = &mut self.site {
            match site.hidden.iter().position(|k| *k == kind) {
                Some(i) => {
                    site.hidden.remove(i);
                }
                None => site.hidden.push(kind),
            }
            cx.notify();
        }
    }

    /// Each measurement's figure, written beside it in the viewport.
    fn measure_labels(&self) -> Vec<gpui::AnyElement> {
        self.measurements
            .iter()
            .filter_map(|m| {
                let s = self.shown(m.which)?;
                let at = self.project(m.anchor(s)?)?;
                let text = m.kind.format(m.value(s)?);
                Some(
                    div()
                        .absolute()
                        .left(at.x + px(6.))
                        .top(at.y - px(9.))
                        .px(px(5.))
                        .py(px(1.))
                        .rounded(px(4.))
                        .bg(theme::ink())
                        .text_color(theme::on_ink())
                        .font_family(theme::MONO)
                        .text_size(px(10.))
                        .whitespace_nowrap()
                        .child(text)
                        .into_any_element(),
                )
            })
            .collect()
    }

    /// While measuring: what kind, what to click next, and every measurement
    /// made, each with its figure and a way to take it away.
    fn measure_card(&self, cx: &mut Context<Self>) -> Option<gpui::AnyElement> {
        if !self.measuring {
            return None;
        }
        let kind = self.measure_kind;
        let need = kind.atoms();
        let hint = match self.pending.len() {
            0 => format!("Click {need} atoms"),
            n => format!("{n} of {need} picked · click empty space to reset"),
        };
        let kinds = div()
            .flex()
            .gap(px(2.))
            .children(measure::Kind::ALL.iter().map(|&k| {
                pill_keyed(("measure-kind", k as usize), k.label(), k == kind)
                    .debug_selector(move || format!("measure-kind-{}", k as usize))
                    .cursor_pointer()
                    .on_click(cx.listener(move |this, _, _, cx| this.set_measure_kind(k, cx)))
                    .into_any_element()
            }));
        let rows = self.measurements.iter().enumerate().filter_map(|(i, m)| {
            let s = self.shown(m.which).or_else(|| {
                // A hidden structure's measurements stay listed.
                if m.which == 0 {
                    self.structure.as_ref()
                } else {
                    self.compare.as_ref()
                }
            })?;
            let tags = m
                .atoms
                .iter()
                .map(|&a| s.atom_tag(a))
                .collect::<Vec<_>>()
                .join(" – ");
            let value = m.kind.format(m.value(s)?);
            Some(
                div()
                    .flex()
                    .items_center()
                    .gap(px(8.))
                    .py(px(3.))
                    .border_t_1()
                    .border_color(theme::border())
                    .child(
                        div()
                            .flex_1()
                            .min_w_0()
                            .text_size(px(10.))
                            .text_color(theme::muted())
                            .child(clip(&tags, 44)),
                    )
                    .child(
                        div()
                            .font_family(theme::MONO)
                            .text_size(px(10.5))
                            .text_color(theme::ink())
                            .child(value),
                    )
                    .child(
                        div()
                            .id(("remove-measurement", i))
                            .debug_selector(move || format!("remove-measurement-{i}"))
                            .px(px(4.))
                            .rounded(px(3.))
                            .text_size(px(11.))
                            .text_color(theme::faint())
                            .cursor_pointer()
                            .hover(|d| d.bg(theme::sunken()).text_color(theme::ink()))
                            .child("×")
                            .on_click(
                                cx.listener(move |this, _, _, cx| this.remove_measurement(i, cx)),
                            ),
                    ),
            )
        });
        Some(
            div()
                .w(px(300.))
                .flex()
                .flex_col()
                .gap(px(6.))
                .p(px(8.))
                .rounded(px(8.))
                .bg(theme::surface())
                .border_1()
                .border_color(theme::border())
                .shadow_sm()
                // Clicks on the card are its own, not the viewport's.
                .on_mouse_down(MouseButton::Left, |_, _, cx| cx.stop_propagation())
                .on_mouse_move(cx.listener(Self::over_card))
                .child(kinds)
                .child(
                    div()
                        .text_size(px(10.))
                        .text_color(theme::faint())
                        .child(hint),
                )
                .child(div().flex().flex_col().children(rows))
                .when(!self.measurements.is_empty(), |d| {
                    d.child(
                        div().flex().justify_end().child(
                            div()
                                .id("clear-measurements")
                                .debug_selector(|| "clear-measurements".into())
                                .text_size(px(9.5))
                                .text_color(theme::muted())
                                .cursor_pointer()
                                .hover(|d| d.text_color(theme::ink()))
                                .child("Clear")
                                .on_click(cx.listener(|this, _, _, cx| {
                                    this.measurements.clear();
                                    this.compose(cx);
                                })),
                        ),
                    )
                })
                .into_any_element(),
        )
    }

    /// What was clicked on, named in the corner until it is let go of.
    fn selected_card(&self, cx: &mut Context<Self>) -> Option<gpui::AnyElement> {
        if self.measuring {
            return None;
        }
        let t = self.selected?;
        let card = self.atom_card(t)?;
        let s = self.shown(t.which)?;
        let residue = s.atoms[t.picked.atom].residue as usize;
        // A ligand's own atom offers its binding site.
        let ligand = t.which == 0
            && s.has_polymer()
            && s.residues[residue].outside_polymer()
            && !crate::protein::is_solvent(&s.residues[residue].name);
        let on_ink = theme::on_ink();
        let soft = |a: f32| gpui::Hsla { a, ..on_ink };
        let action = |id: &'static str, label: &'static str| {
            div()
                .id(id)
                .debug_selector(move || id.into())
                .px(px(8.))
                .py(px(3.))
                .rounded(px(5.))
                .border_1()
                .border_color(soft(0.18))
                .bg(soft(0.06))
                .text_size(px(10.))
                .text_color(soft(0.85))
                .cursor_pointer()
                .hover(|d| d.bg(soft(0.14)).border_color(soft(0.3)).text_color(on_ink))
                .child(label)
        };
        Some(
            card.child(div().h(px(1.)).w_full().bg(soft(0.12)))
                .child(
                    div()
                        .flex()
                        .items_center()
                        .gap(px(5.))
                        .child(
                            action("focus-selected", "Focus")
                                .on_click(cx.listener(|this, _, _, cx| this.focus_selected(cx))),
                        )
                        .when(t.which == 0, |d| {
                            let labelled =
                                self.labels.iter().any(|l| l.atom as usize == t.picked.atom);
                            d.child(
                                action(
                                    "label-selected",
                                    if labelled { "Unlabel" } else { "Label" },
                                )
                                .on_click(cx.listener(|this, _, _, cx| this.toggle_label(cx))),
                            )
                        })
                        .when(ligand, |d| {
                            d.child(action("site-selected", "Binding site").on_click(
                                cx.listener(move |this, _, _, cx| this.show_site(residue, cx)),
                            ))
                        }),
                )
                .into_any_element(),
        )
    }

    /// The binding site looked at: which ligand, how many residues around
    /// it, a way on to the next copy, and a way out.
    fn site_card(&self, cx: &mut Context<Self>) -> Option<gpui::AnyElement> {
        let site = self.site.as_ref()?;
        let s = self.shown(0)?;
        let r = &s.residues[site.ligand];
        let chain = s
            .chains
            .get(r.chain as usize)
            .map(|c| c.id.trim().to_string())
            .unwrap_or_default();
        let copies = s.copies(&r.name);
        let at = copies.iter().position(|&c| c == site.ligand).unwrap_or(0);
        let name = r.name.clone();
        let link = |id: &'static str, label: String| {
            div()
                .id(id)
                .debug_selector(move || id.into())
                .px(px(6.))
                .py(px(1.))
                .rounded(px(4.))
                .text_size(px(10.))
                .text_color(theme::muted())
                .cursor_pointer()
                .hover(|d| d.bg(theme::sunken()).text_color(theme::ink()))
                .child(label)
        };
        Some(
            div()
                .w(px(300.))
                .flex()
                .flex_col()
                .gap(px(3.))
                .p(px(8.))
                .rounded(px(8.))
                .bg(theme::surface())
                .border_1()
                .border_color(theme::border())
                .shadow_sm()
                .on_mouse_down(MouseButton::Left, |_, _, cx| cx.stop_propagation())
                .on_mouse_move(cx.listener(Self::over_card))
                .child(
                    div()
                        .flex()
                        .items_center()
                        .gap(px(6.))
                        .child(
                            div()
                                .flex_1()
                                .text_size(px(11.5))
                                .text_color(theme::ink())
                                .child(if chain.is_empty() {
                                    format!("{} {}", r.name, r.seq)
                                } else {
                                    format!("{} {} · chain {chain}", r.name, r.seq)
                                }),
                        )
                        .when(copies.len() > 1, |d| {
                            d.child(
                                link("next-site", format!("{} of {} →", at + 1, copies.len()))
                                    .on_click(
                                        cx.listener(move |this, _, _, cx| {
                                            this.next_site(&name, cx)
                                        }),
                                    ),
                            )
                        })
                        .child(
                            link("close-site", "×".into())
                                .on_click(cx.listener(|this, _, _, cx| this.close_site(cx))),
                        ),
                )
                .child(
                    div()
                        .text_size(px(10.))
                        .text_color(theme::faint())
                        .child(format!(
                            "{} residues within {SITE_CUTOFF:.0} Å",
                            site.residues.len()
                        )),
                )
                .child(
                    div()
                        .flex()
                        .flex_wrap()
                        .gap_x(px(10.))
                        .gap_y(px(2.))
                        .children(contacts::Kind::ALL.iter().filter_map(|&kind| {
                            let n = site.contacts.iter().filter(|c| c.kind == kind).count();
                            if n == 0 {
                                return None;
                            }
                            let shown = !site.hidden.contains(&kind);
                            Some(
                                div()
                                    .id(("contact-kind", kind as usize))
                                    .debug_selector(move || format!("contacts-{kind:?}"))
                                    .flex()
                                    .items_center()
                                    .gap(px(5.))
                                    .cursor_pointer()
                                    .text_size(px(10.))
                                    .text_color(if shown { theme::text() } else { theme::faint() })
                                    .child(
                                        div()
                                            .w(px(12.))
                                            .h(px(3.))
                                            .rounded(px(1.5))
                                            .when(shown, |d| d.bg(gpui::rgb(kind.color())))
                                            .when(!shown, |d| {
                                                d.border_1().border_color(theme::faint())
                                            }),
                                    )
                                    .child(format!("{n} {}", kind.label(n)))
                                    .tooltip(move |_, cx| {
                                        cx.new(|_| {
                                            Tip(if shown { "Hide these" } else { "Show these" }
                                                .into())
                                        })
                                        .into()
                                    })
                                    .on_click(cx.listener(move |this, _, _, cx| {
                                        this.toggle_contacts(kind, cx)
                                    })),
                            )
                        })),
                )
                .into_any_element(),
        )
    }

    /// The binding site's residues, each named beside its side chain.
    /// The labels pinned to atoms, where the atoms are drawn.
    fn atom_labels(&self) -> Vec<gpui::AnyElement> {
        let Some(s) = self.shown(0) else {
            return Vec::new();
        };
        self.labels
            .iter()
            .filter_map(|l| {
                let at = self.project(s.atoms.get(l.atom as usize)?.pos)?;
                Some(
                    div()
                        .absolute()
                        .left(at.x + px(4.))
                        .top(at.y - px(8.))
                        .px(px(4.))
                        .rounded(px(3.))
                        .bg(theme::surface())
                        .opacity(0.9)
                        .text_color(theme::ink())
                        .text_size(px(10.5))
                        .whitespace_nowrap()
                        .child(l.text.clone())
                        .into_any_element(),
                )
            })
            .collect()
    }

    /// The words a residue is labelled with: its name and number, and its
    /// chain where there are several.
    fn residue_label(s: &Structure, r: usize) -> String {
        let res = &s.residues[r];
        let chain = s
            .chains
            .get(res.chain as usize)
            .map(|c| c.id.trim())
            .filter(|c| !c.is_empty() && s.chains.len() > 1);
        match chain {
            Some(c) => format!("{} {}:{c}", res.name, res.seq),
            None => format!("{} {}", res.name, res.seq),
        }
    }

    /// Label what `selection` takes: a residue each, at its alpha carbon or
    /// first atom; with `atoms`, each atom by name; with `text`, the first
    /// atom with that. How many labels it added.
    pub(super) fn add_labels(
        &mut self,
        selection: &Selection,
        atoms: bool,
        text: Option<String>,
        cx: &mut Context<Self>,
    ) -> usize {
        let Some(s) = self.structure.clone() else {
            return 0;
        };
        let taken = selection.atoms(&s);
        let chosen: Vec<usize> = taken
            .iter()
            .enumerate()
            .filter(|(_, t)| **t)
            .map(|(i, _)| i)
            .collect();
        let before = self.labels.len();
        if let Some(text) = text {
            if let Some(&a) = chosen.first() {
                self.labels.push(Label {
                    atom: a as u32,
                    text,
                });
            }
        } else if atoms || !s.has_polymer() {
            for a in chosen {
                let name = s.atoms[a].name();
                let text = if name.is_empty() {
                    format!("{} {}", s.atoms[a].element.symbol(), a + 1)
                } else {
                    name.to_string()
                };
                self.labels.push(Label {
                    atom: a as u32,
                    text,
                });
            }
        } else {
            let mut seen = std::collections::HashSet::new();
            for a in chosen {
                let r = s.atoms[a].residue as usize;
                if seen.insert(r) {
                    let at = s.residues[r].ca.map_or(a, |c| c as usize);
                    self.labels.push(Label {
                        atom: at as u32,
                        text: Self::residue_label(&s, r),
                    });
                }
            }
        }
        // One label to an atom: a new one replaces what was there.
        let mut seen = std::collections::HashSet::new();
        let mut kept: Vec<Label> = self
            .labels
            .drain(..)
            .rev()
            .filter(|l| seen.insert(l.atom))
            .collect();
        kept.reverse();
        let added = kept.len().saturating_sub(before.min(kept.len()));
        self.labels = kept;
        cx.notify();
        added
    }

    /// Label the picked residue or atom, or take its label off again.
    fn toggle_label(&mut self, cx: &mut Context<Self>) {
        let Some(t) = self.selected.filter(|t| t.which == 0) else {
            return;
        };
        let Some(s) = self.structure.clone() else {
            return;
        };
        let atom = t.picked.atom as u32;
        if let Some(i) = self.labels.iter().position(|l| l.atom == atom) {
            self.labels.remove(i);
        } else {
            let r = s.atoms[t.picked.atom].residue as usize;
            let text = if t.picked.residue || s.residues[r].ca == Some(atom) {
                Self::residue_label(&s, r)
            } else {
                format!(
                    "{} {}",
                    s.atoms[t.picked.atom].name(),
                    Self::residue_label(&s, r)
                )
            };
            self.labels.push(Label { atom, text });
        }
        cx.notify();
    }

    fn site_labels(&self) -> Vec<gpui::AnyElement> {
        let (Some(site), Some(s)) = (&self.site, self.shown(0)) else {
            return Vec::new();
        };
        site.residues
            .iter()
            .filter_map(|&r| {
                let at = self.project(s.side_chain_centre(r)?)?;
                let res = &s.residues[r];
                Some(
                    div()
                        .absolute()
                        .left(at.x + px(4.))
                        .top(at.y - px(7.))
                        .px(px(4.))
                        .rounded(px(3.))
                        .bg(theme::surface())
                        .opacity(0.9)
                        .text_color(theme::ink())
                        .text_size(px(9.5))
                        .whitespace_nowrap()
                        .child(format!("{}{}", res.name, res.seq))
                        .into_any_element(),
                )
            })
            .collect()
    }

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
        if self.turning.is_some() {
            self.turning = None;
            self.turn_task = None;
            cx.notify();
        }
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
            .on_mouse_move(cx.listener(Self::over_card))
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
                pill(
                    "play",
                    if self.playing { "Pause" } else { "Play" },
                    self.playing,
                )
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
        let vertical = theme::vertical();
        let width = f32::from(window.viewport_size().width);
        self.narrow = width < 2. * PANEL_WIDTH + 40.;
        let stage = if vertical {
            width
        } else {
            width - theme::side_width()
        };
        self.compact = stage < 640.;

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
            .when(!folding || structure.is_some(), |d| {
                d.child(self.toolbar(cx))
            })
            .child(
                div()
                    .flex()
                    .when(vertical, |d| d.flex_col())
                    .flex_1()
                    .min_h_0()
                    .overflow_hidden()
                    .children({
                        let panel = self
                            .panel(folding, focused, structure, vertical, cx)
                            .into_any_element();
                        let divider = self.divider(vertical, cx).into_any_element();
                        let stage = self.stage(cx).into_any_element();
                        if vertical {
                            // The name of the app over the viewport, as a
                            // page's masthead is over its picture.
                            let masthead = div()
                                .flex_none()
                                .w_full()
                                .pt(px(14.))
                                .pb(px(12.))
                                .bg(theme::surface())
                                .border_b_1()
                                .border_color(theme::border())
                                .child(self.masthead())
                                .into_any_element();
                            vec![masthead, stage, divider, panel]
                        } else {
                            vec![panel, divider, stage]
                        }
                    }),
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
        view.update(cx, |this, _| {
            this.scene.spheres.iter().map(|s| s.color).collect()
        })
    }

    /// Where in the window the atom nearest the viewer is drawn.
    fn front_atom(view: &Entity<Viewer>, cx: &mut VisualTestContext) -> Point<Pixels> {
        view.update(cx, |this, _| {
            let s = this.structure.as_ref().unwrap();
            let cam = this.camera;
            let m = cam.basis();
            let cache = this.cache.borrow();
            let (w, h) = cache.size;
            let k = w as f32 / f32::from(cache.bounds.size.width);
            let view = |p: [f32; 3]| {
                let q = [
                    p[0] - cam.center[0],
                    p[1] - cam.center[1],
                    p[2] - cam.center[2],
                ];
                let d = |r: [f32; 3]| r[0] * q[0] + r[1] * q[1] + r[2] * q[2];
                [d(m[0]), d(m[1]), d(m[2])]
            };
            let v = s
                .atoms
                .iter()
                .map(|a| view(a.pos))
                .max_by(|a, b| a[2].total_cmp(&b[2]))
                .unwrap();
            let scale = w.min(h) as f32 / (2.2 * cam.radius) * cam.zoom;
            let x = w as f32 * 0.5 + cam.pan[0] + v[0] * scale;
            let y = h as f32 * 0.5 + cam.pan[1] - v[1] * scale;
            cache.bounds.origin + gpui::point(px(x / k), px(y / k))
        })
    }

    #[gpui::test]
    fn pointing_at_an_atom_names_it_and_clicking_picks_it(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        let at = front_atom(&view, cx);
        cx.simulate_mouse_move(at, None, Modifiers::none());
        cx.run_until_parked();
        let hover = view.update(cx, |this, _| this.hover.map(|(t, _)| t));
        let hover = hover.expect("nothing under the pointer");
        assert_eq!(hover.which, 0);

        let glows = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| {
                this.scene.spheres.iter().filter(|s| s.alpha < 1.0).count()
            })
        };
        assert_eq!(glows(&view, cx), 0);
        cx.simulate_click(at, Modifiers::none());
        cx.run_until_parked();
        assert_eq!(view.update(cx, |this, _| this.selected), Some(hover));
        assert_eq!(glows(&view, cx), 1, "the picked atom does not glow");

        // A drag turns the view and picks nothing.
        cx.simulate_mouse_down(at, MouseButton::Left, Modifiers::none());
        cx.simulate_mouse_move(
            at + gpui::point(px(40.), px(0.)),
            MouseButton::Left,
            Modifiers::none(),
        );
        cx.simulate_mouse_up(
            at + gpui::point(px(40.), px(0.)),
            MouseButton::Left,
            Modifiers::none(),
        );
        cx.run_until_parked();
        assert_eq!(view.update(cx, |this, _| this.selected), Some(hover));

        // Clicking empty space lets go.
        let corner = view.update(cx, |this, _| this.cache.borrow().bounds.origin)
            + gpui::point(px(3.), px(120.));
        cx.simulate_click(corner, Modifiers::none());
        cx.run_until_parked();
        assert_eq!(view.update(cx, |this, _| this.selected), None);
        assert_eq!(glows(&view, cx), 0);
    }

    #[gpui::test]
    fn command_click_picks_an_atom_and_flies_to_it(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        let at = front_atom(&view, cx);
        let before = view.update(cx, |this, _| this.camera.center);
        cx.simulate_click(at, Modifiers::secondary_key());
        cx.executor().advance_clock(Duration::from_millis(1500));
        cx.run_until_parked();
        let (picked, centre) = view.update(cx, |this, _| {
            let t = this.selected.expect("nothing picked");
            (
                this.structure.as_ref().unwrap().atoms[t.picked.atom].pos,
                this.camera.center,
            )
        });
        assert_ne!(centre, before, "the view did not move");
        let off = (0..3)
            .map(|i| (centre[i] - picked[i]).powi(2))
            .sum::<f32>()
            .sqrt();
        assert!(off < 0.5, "the view is {off} Å from the atom");

        // Again on the same atom keeps it picked, rather than letting go.
        cx.simulate_click(at, Modifiers::secondary_key());
        cx.run_until_parked();
        assert!(view.update(cx, |this, _| this.selected.is_some()));
    }

    #[gpui::test]
    fn the_recording_light_covers_none_of_the_controls(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        click(cx, "tab-film");
        click(cx, "movie-record-menu");
        assert!(view.update(cx, |this, _| this.recording()));
        let light = cx.debug_bounds("recording").expect("no recording light");
        let controls = cx.debug_bounds("top-controls").expect("no controls");
        assert!(
            !light.intersects(&controls),
            "the light {light:?} sits on the controls {controls:?}"
        );
        click(cx, "recording");
        assert!(!view.update(cx, |this, _| this.recording()));
    }

    #[gpui::test]
    fn a_box_is_replicated_wrapped_and_put_back_from_the_console(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        let melt = SAMPLES
            .iter()
            .find(|s| s.file.ends_with("lammpstrj.gz"))
            .unwrap();
        view.update(cx, |this, cx| {
            let s = formats::load_bytes(melt.file, melt.bytes.to_vec()).unwrap();
            this.show(s, cx);
        });
        cx.run_until_parked();
        let atoms = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| {
                this.structure.as_ref().map_or(0, |s| s.atoms.len())
            })
        };
        let n = atoms(&view, cx);
        view.update(cx, |this, cx| this.run_command("replicate 2 1 1", cx))
            .unwrap();
        assert_eq!(atoms(&view, cx), 2 * n);
        let boundary = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| this.boundary)
        };
        // Placing atoms one way stands in for another, and the copies stay.
        view.update(cx, |this, cx| this.run_command("unwrap", cx))
            .unwrap();
        assert_eq!(boundary(&view, cx).place, Place::InTime);
        view.update(cx, |this, cx| this.run_command("wrap", cx))
            .unwrap();
        assert_eq!(boundary(&view, cx).place, Place::Wrap);
        assert_eq!(atoms(&view, cx), 2 * n);
        // Asked again, what is on goes off.
        assert_eq!(boundary(&view, cx).replicate, Some([2, 1, 1]));
        view.update(cx, |this, cx| this.run_command("wrap", cx))
            .unwrap();
        assert_eq!(boundary(&view, cx).place, Place::AsRead);
        // Wrapped atoms are inside the cell, and back where the file had
        // them once wrapping is off.
        view.update(cx, |this, cx| this.run_command("replicate 2 1 1", cx))
            .unwrap();
        assert_eq!(atoms(&view, cx), n);
        let as_read = view.update(cx, |this, _| this.structure.as_ref().unwrap().atoms[0].pos);
        view.update(cx, |this, cx| this.run_command("wrap", cx))
            .unwrap();
        view.update(cx, |this, cx| this.run_command("wrap", cx))
            .unwrap();
        assert_eq!(
            view.update(cx, |this, _| this.structure.as_ref().unwrap().atoms[0].pos),
            as_read
        );
        view.update(cx, |this, cx| this.run_command("replicate 3 1 1", cx))
            .unwrap();
        view.update(cx, |this, cx| this.run_command("original", cx))
            .unwrap();
        assert_eq!(atoms(&view, cx), n);
        assert_eq!(boundary(&view, cx), Boundary::default());
        assert!(view
            .update(cx, |this, cx| this.run_command("original", cx))
            .is_err());
        assert!(view
            .update(cx, |this, cx| this.run_command("replicate 2 2", cx))
            .is_err());
        // A box of plain beads has no properties to colour by.
        assert!(view
            .update(cx, |this, cx| this.run_command("color property", cx))
            .is_err());
    }

    /// A cube file of one carbon with a p orbital about it: positive above,
    /// negative below.
    fn p_orbital_cube() -> String {
        let n = 20;
        let mut text = String::from("p orbital\ntest\n    1   -4.75   -4.75   -4.75\n");
        for k in 0..3 {
            let mut v = [0.0; 3];
            v[k] = 0.5;
            text += &format!("   {n}   {:.3}   {:.3}   {:.3}\n", v[0], v[1], v[2]);
        }
        text += "    6    6.0    0.0    0.0    0.0\n";
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let (x, y, z) = [i, j, k].map(|q| -4.75 + q as f32 * 0.5).into();
                    let r2 = x * x + y * y + z * z;
                    text += &format!("{:.5e} ", z * (-r2 / 4.0).exp());
                }
                text += "\n";
            }
        }
        text
    }

    #[gpui::test]
    fn an_orbital_shows_two_lobes_and_a_map_opens_alone(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        let s = formats::load_bytes("orbital.cube", p_orbital_cube().into_bytes()).unwrap();
        let volume = s.volume.clone().expect("the cube's grid is read");
        assert_eq!(volume.dims, [20, 20, 20]);
        assert!(volume.signed());
        view.update(cx, |this, cx| this.show(s, cx));
        cx.run_until_parked();
        let colours = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| {
                let mut c: Vec<[u32; 3]> = this
                    .base
                    .tris
                    .iter()
                    .map(|t| t.c.map(|x| (x * 100.0) as u32))
                    .collect();
                c.dedup();
                c.sort();
                c.dedup();
                c.len()
            })
        };
        assert_eq!(colours(&view, cx), 2, "a blue lobe and a red one");
        view.update(cx, |this, cx| this.run_command("iso mesh", cx))
            .unwrap();
        assert!(view.update(cx, |this, _| this.base.tris.is_empty()
            && this.base.capsules.len() > 100));
        view.update(cx, |this, cx| this.run_command("iso off", cx))
            .unwrap();
        // The grid's box is all that is left of rods.
        assert!(view.update(cx, |this, _| this.base.capsules.len() <= 12
            && this.base.tris.is_empty()));
        assert!(view
            .update(cx, |this, cx| this.run_command("iso high", cx))
            .is_err());

        // A map with no atoms opens on its own, as its surface in its box.
        let mut map = vec![0u8; 1024];
        for (o, v) in [
            (0, 8i32),
            (4, 8),
            (8, 8),
            (12, 2),
            (28, 8),
            (32, 8),
            (36, 8),
            (64, 1),
            (68, 2),
            (72, 3),
        ] {
            map[o..o + 4].copy_from_slice(&v.to_le_bytes());
        }
        for (o, v) in [(40, 8.0f32), (44, 8.0), (48, 8.0)] {
            map[o..o + 4].copy_from_slice(&v.to_le_bytes());
        }
        for k in 0..8 {
            for j in 0..8 {
                for i in 0..8 {
                    let d = [i, j, k].map(|q| q as f32 - 3.5);
                    map.extend(
                        (4.0 - (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()).to_le_bytes(),
                    );
                }
            }
        }
        let alone = formats::load_bytes("blob.mrc", map).unwrap();
        assert!(alone.atoms.is_empty() && alone.volume.is_some());
        view.update(cx, |this, cx| this.finish(Ok(alone), cx));
        cx.run_until_parked();
        assert!(view.update(cx, |this, _| !this.base.tris.is_empty()));
    }

    #[gpui::test]
    fn a_measure_is_worked_out_only_when_picked_and_of_the_atoms_chosen(cx: &mut TestAppContext) {
        use crate::analysis::{Measure, Plot};
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        let crambin = formats::load_bytes(SAMPLES[2].file, SAMPLES[2].bytes.to_vec()).unwrap();
        let frames: Vec<crate::protein::Snapshot> = (0..6)
            .map(|k| crate::protein::Snapshot {
                step: None,
                pos: crambin
                    .atoms
                    .iter()
                    .map(|a| [a.pos[0] * (1.0 + 0.02 * k as f32), a.pos[1], a.pos[2]])
                    .collect(),
                cell: None,
            })
            .collect();
        let moving = crambin.with_frames(frames).unwrap();
        view.update(cx, |this, cx| this.show(moving, cx));
        cx.run_until_parked();
        let settle = |cx: &mut VisualTestContext| {
            for _ in 0..20 {
                cx.executor().advance_clock(Duration::from_millis(300));
                cx.run_until_parked();
            }
        };
        // Opened, nothing is worked out until a measure is picked.
        view.update(cx, |this, cx| this.run_command("analyse", cx))
            .unwrap();
        settle(cx);
        assert!(view.update(cx, |this, _| this.analysis.is_empty()
            && this.analysis_task.is_none()));
        // A protein is measured by its alpha carbons to begin with.
        let (selection, count) = view.update(cx, |this, _| {
            let s = this.structure.clone().unwrap();
            (this.analysis_selection(&s), this.analysis_atoms(&s).len())
        });
        assert_eq!(selection, "(protein and name CA)");
        assert_eq!(count, 46);
        view.update(cx, |this, cx| this.run_command("analyse rmsd", cx))
            .unwrap();
        settle(cx);
        let rmsd = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| {
                let s = this.structure.clone().unwrap();
                this.analysis
                    .get(&this.analysis_key(&s, Measure::Rmsd))
                    .cloned()
            })
        };
        // Stretched further each frame, it drifts further each frame.
        let Some(plot) = rmsd(&view, cx) else {
            panic!("worked out")
        };
        let Plot::Time(points) = &*plot else { panic!() };
        assert!(points.windows(2).all(|w| w[1].1 > w[0].1));
        assert_eq!(view.update(cx, |this, _| this.analysis.len()), 1);
        // Of other atoms, it is worked out again for them.
        view.update(cx, |this, cx| {
            this.run_command("analyse rmsd of elem S", cx)
        })
        .unwrap();
        settle(cx);
        assert_eq!(view.update(cx, |this, _| this.analysis.len()), 2);
        assert_eq!(
            view.update(cx, |this, _| this.analysis_custom.clone()),
            Some("elem S".into())
        );
        assert!(view
            .update(cx, |this, cx| this.run_command("analyse rg of resn", cx))
            .is_err());
        // The plot is laid out, and a click near its right end goes to the
        // last frame.
        assert!(cx.debug_bounds("analysis-plot").is_some());
        view.update(cx, |this, cx| this.go_to_plotted(0.98, cx));
        assert_eq!(
            view.update(cx, |this, _| this.structure.as_ref().unwrap().frame),
            5
        );
        // Fitting changes the frames, so what is shown is worked out again.
        view.update(cx, |this, cx| this.run_command("fit", cx))
            .unwrap();
        assert!(view.update(cx, |this, _| this.analysis_task.is_some()));
    }

    #[test]
    fn only_a_periodic_box_offers_the_cell_tools() {
        let load = |name: &str| {
            let sample = SAMPLES.iter().find(|s| s.file.starts_with(name)).unwrap();
            formats::load_bytes(sample.file, sample.bytes.to_vec()).unwrap()
        };
        // Haemoglobin's crystal cell is not a box to work across.
        let haemoglobin = load("4hhb");
        assert!(haemoglobin.cell.is_some() && !super::periodic(&haemoglobin));
        assert!(super::periodic(&load("lj-melt")));
        assert!(super::periodic(&load("rocksalt")));
        assert!(!super::periodic(&load("caffeine")));
    }

    #[gpui::test]
    fn residues_and_atoms_are_labelled_from_the_console(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        view.update(cx, |this, cx| {
            let s = formats::load_bytes(SAMPLES[2].file, SAMPLES[2].bytes.to_vec()).unwrap();
            this.show(s, cx);
        });
        cx.run_until_parked();
        let labels = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| this.labels.clone())
        };
        view.update(cx, |this, cx| this.run_command("label resi 1-3", cx))
            .unwrap();
        let l = labels(&view, cx);
        assert_eq!(l.len(), 3);
        assert_eq!(l[0].text, "THR 1");
        view.update(cx, |this, cx| {
            this.run_command("label atoms resi 5 and name CA+CB", cx)
        })
        .unwrap();
        assert_eq!(labels(&view, cx).len(), 5);
        // Text of one's own, on the first atom taken; it replaces a label there.
        view.update(cx, |this, cx| {
            this.run_command("label resi 1 and name CA \"the start\"", cx)
        })
        .unwrap();
        let l = labels(&view, cx);
        assert_eq!(l.len(), 5);
        assert!(l.iter().any(|l| l.text == "the start"));
        view.update(cx, |this, cx| this.run_command("label clear", cx))
            .unwrap();
        assert!(labels(&view, cx).is_empty());
        assert!(view
            .update(cx, |this, cx| this.run_command("label", cx))
            .is_err());
    }

    #[gpui::test]
    fn the_sequence_and_the_viewport_point_at_each_other(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        view.update(cx, |this, cx| {
            let s = formats::load_bytes(SAMPLES[2].file, SAMPLES[2].bytes.to_vec()).unwrap();
            this.show(s, cx);
        });
        cx.run_until_parked();
        let glows = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            // In Cartoon a residue is painted on the ribbon: what differs
            // from the scene with nothing marked.
            view.update(cx, |this, _| {
                let (now, base) = (&this.scene, &this.base);
                now.tris
                    .iter()
                    .zip(&base.tris)
                    .filter(|(a, b)| a.c != b.c)
                    .count()
                    + now.spheres.iter().filter(|s| s.alpha < 1.0).count()
            })
        };
        assert_eq!(glows(&view, cx), 0);

        // Pointing at a letter glows its residue, and leaving it lets go.
        let letter = cx.debug_bounds("residue-3").expect("no sequence on screen");
        cx.simulate_mouse_move(letter.center(), None, Modifiers::none());
        cx.run_until_parked();
        assert_eq!(view.update(cx, |this, _| this.seq_hover), Some(3));
        assert!(glows(&view, cx) > 0, "the residue pointed at does not glow");
        let away = view.update(cx, |this, _| this.cache.borrow().bounds.origin)
            + gpui::point(px(3.), px(120.));
        cx.simulate_mouse_move(away, None, Modifiers::none());
        cx.run_until_parked();
        assert_eq!(view.update(cx, |this, _| this.seq_hover), None);

        // Clicking it picks the residue out, as a click on its ribbon would.
        click(cx, "residue-3");
        let picked = view.update(cx, |this, _| this.selected_residue());
        assert_eq!(picked, Some(3));
        assert!(view.update(cx, |this, _| this.selected.unwrap().picked.residue));
        assert!(glows(&view, cx) > 0);

        // And what is picked in the viewport is the one marked here.
        let at = front_atom(&view, cx);
        cx.simulate_click(at, Modifiers::none());
        cx.run_until_parked();
        let (picked, hit) = view.update(cx, |this, _| {
            let s = this.structure.as_ref().unwrap();
            (
                this.selected_residue(),
                this.selected
                    .map(|t| s.atoms[t.picked.atom].residue as usize),
            )
        });
        assert_eq!(picked, hit);
    }

    /// Where in the window atom `i` is drawn.
    fn atom_on_screen(
        view: &Entity<Viewer>,
        cx: &mut VisualTestContext,
        i: usize,
    ) -> Point<Pixels> {
        view.update(cx, |this, _| {
            let p = this.structure.as_ref().unwrap().atoms[i].pos;
            this.cache.borrow().bounds.origin + this.project(p).unwrap()
        })
    }

    #[gpui::test]
    fn measuring_takes_two_clicks_for_a_distance(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        click(cx, "measure");
        assert!(view.update(cx, |this, _| this.measuring));

        // Two atoms that are each the nearest under their own point.
        let first = front_atom(&view, cx);
        let a = view
            .update(cx, |this, _| this.pick_at(first))
            .unwrap()
            .picked
            .atom;
        let atoms = view.update(cx, |this, _| this.structure.as_ref().unwrap().atoms.len());
        let mut found = None;
        for i in (0..atoms).filter(|&i| i != a) {
            let at = atom_on_screen(&view, cx, i);
            if view
                .update(cx, |this, _| this.pick_at(at))
                .is_some_and(|t| t.picked.atom == i)
            {
                found = Some((i, at));
                break;
            }
        }
        let (b, second) = found.expect("no second atom to click");

        cx.simulate_click(first, Modifiers::none());
        cx.run_until_parked();
        assert_eq!(view.update(cx, |this, _| this.pending.clone()), vec![a]);
        assert!(
            view.update(cx, |this, _| this.selected.is_none()),
            "measuring picked the atom out"
        );
        cx.simulate_click(second, Modifiers::none());
        cx.run_until_parked();

        let (made, value, expected) = view.update(cx, |this, _| {
            let s = this.structure.as_ref().unwrap();
            let m = this.measurements.first().cloned();
            let (p, q) = (s.atoms[a].pos, s.atoms[b].pos);
            let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
            (m.clone(), m.and_then(|m| m.value(s)), d)
        });
        let made = made.expect("no measurement made");
        assert_eq!(made.atoms, vec![a, b]);
        assert!((value.unwrap() - expected).abs() < 1e-4);
        assert!(view.update(cx, |this, _| this.pending.is_empty()));
        // Its line is drawn, and its figure is written.
        assert_eq!(view.update(cx, |this, _| this.measure_segments().len()), 1);
        assert_eq!(view.update(cx, |this, _| this.measure_labels().len()), 1);

        click(cx, "remove-measurement-0");
        assert!(view.update(cx, |this, _| this.measurements.is_empty()));
        assert!(view.update(cx, |this, _| this.measure_segments().is_empty()));

        // An angle wants a third.
        click(cx, "measure-kind-1");
        assert_eq!(
            view.update(cx, |this, _| this.measure_kind),
            measure::Kind::Angle
        );
        cx.simulate_click(first, Modifiers::none());
        cx.simulate_click(second, Modifiers::none());
        cx.run_until_parked();
        assert!(view.update(cx, |this, _| this.measurements.is_empty()
            && this.pending.len() == 2));

        // Done measuring, a click picks atoms out again.
        click(cx, "measure");
        cx.simulate_click(first, Modifiers::none());
        cx.run_until_parked();
        assert!(view.update(cx, |this, _| this.selected.is_some()));
    }

    #[gpui::test]
    fn a_ligand_shows_its_binding_site(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        view.update(cx, |this, cx| {
            let s = formats::load_bytes(SAMPLES[5].file, SAMPLES[5].bytes.to_vec()).unwrap();
            this.show(s, cx);
        });
        cx.run_until_parked();
        let zoom = view.update(cx, |this, _| this.camera.zoom);
        let extra = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| {
                this.scene.capsules.len() - this.base.capsules.len()
            })
        };

        click(cx, "ligand-site-HEM");
        let (first, labels) =
            view.update(cx, |this, _| (this.site.clone(), this.site_labels().len()));
        let first = first.expect("no binding site shown");
        assert!(!first.residues.is_empty());
        assert_eq!(labels, first.residues.len());
        assert!(extra(&view, cx) > 0, "its side chains are not drawn");
        // The camera has come in on it, once it has flown there.
        cx.executor().advance_clock(Duration::from_secs(1));
        cx.run_until_parked();
        let (zoom_now, centre) = view.update(cx, |this, _| (this.camera.zoom, this.camera.center));
        assert!(zoom_now > zoom, "the camera did not come closer");
        let haem = view.update(cx, |this, _| {
            let s = this.structure.as_ref().unwrap();
            spread(
                s.residue_atoms(first.ligand)
                    .iter()
                    .map(|&i| s.atoms[i].pos),
            )
            .0
        });
        assert!((0..3).all(|k| (centre[k] - haem[k]).abs() < 1e-3));

        // What holds it: the iron on its histidine among them, drawn, and
        // hidden from the card.
        let (kinds, lines) = view.update(cx, |this, _| {
            let site = this.site.as_ref().unwrap();
            (
                site.contacts.iter().map(|c| c.kind).collect::<Vec<_>>(),
                this.contact_segments().len(),
            )
        });
        assert!(kinds.contains(&contacts::Kind::Metal), "{kinds:?}");
        assert_eq!(lines, kinds.len());
        click(cx, "contacts-Metal");
        let metals = kinds
            .iter()
            .filter(|k| **k == contacts::Kind::Metal)
            .count();
        assert_eq!(
            view.update(cx, |this, _| this.contact_segments().len()),
            lines - metals
        );
        click(cx, "contacts-Metal");
        assert_eq!(
            view.update(cx, |this, _| this.contact_segments().len()),
            lines
        );

        // On to the next haem, and back out.
        click(cx, "next-site");
        let second = view.update(cx, |this, _| this.site.clone()).unwrap();
        assert_ne!(second.ligand, first.ligand);
        click(cx, "close-site");
        assert!(view.update(cx, |this, _| this.site.is_none()));
        assert_eq!(extra(&view, cx), 0);
    }

    const DIMER: &str = "REMARK 350 BIOMOLECULE: 1\nREMARK 350 APPLY THE FOLLOWING TO CHAINS: A\nREMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000\nREMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000\nREMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000\nREMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000       30.00000\nREMARK 350   BIOMT2   2  0.000000 -1.000000  0.000000        0.00000\nREMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000\nATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00 10.00           N\nATOM      2  CA  GLY A   1       2.450   2.000   3.000  1.00 10.00           C\nATOM      3  C   GLY A   1       3.000   3.400   3.000  1.00 10.00           C\nATOM      4  O   GLY A   1       2.300   4.400   3.000  1.00 10.00           O\nATOM      5  N   GLY A   2       4.800   2.000   3.000  1.00 10.00           N\nATOM      6  CA  GLY A   2       6.250   2.000   3.000  1.00 10.00           C\nATOM      7  C   GLY A   2       6.800   3.400   3.000  1.00 10.00           C\nATOM      8  O   GLY A   2       6.100   4.400   3.000  1.00 10.00           O\nATOM      9  N   GLY A   3       8.600   2.000   3.000  1.00 10.00           N\nATOM     10  CA  GLY A   3      10.050   2.000   3.000  1.00 10.00           C\nATOM     11  C   GLY A   3      10.600   3.400   3.000  1.00 10.00           C\nATOM     12  O   GLY A   3       9.900   4.400   3.000  1.00 10.00           O\n";

    #[gpui::test]
    fn an_assembly_is_shown_in_place_of_the_deposited_chains(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        view.update(cx, |this, cx| {
            let s =
                formats::parse(DIMER, formats::Format::Pdb, "dimer".into(), Source::File).unwrap();
            this.show(s, cx);
        });
        cx.run_until_parked();
        let chains = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| this.structure.as_ref().unwrap().chains.len())
        };
        assert_eq!(chains(&view, cx), 1);
        click(cx, "assembly-0");
        assert_eq!(chains(&view, cx), 2);
        assert_eq!(view.update(cx, |this, _| this.assembly), Some(0));
        // Both copies are drawn: twice the ribbon.
        let tris = view.update(cx, |this, _| this.base.tris.len());
        click(cx, "assembly-deposited");
        assert_eq!(chains(&view, cx), 1);
        assert!(view.update(cx, |this, _| this.deposited.is_none()));
        assert_eq!(view.update(cx, |this, _| this.base.tris.len()) * 2, tris);
    }

    #[gpui::test]
    fn another_structure_is_laid_over_the_one_shown(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        let haemoglobin = formats::load_bytes(SAMPLES[5].file, SAMPLES[5].bytes.to_vec()).unwrap();
        view.update(cx, |this, cx| this.show(haemoglobin.clone(), cx));
        cx.run_until_parked();

        click(cx, "tab-open");
        click(cx, "overlay-entry");
        assert!(view.update(cx, |this, _| this.overlay_next));
        assert!(
            cx.debug_bounds("overlay-cancel").is_some(),
            "nothing says the next is overlaid"
        );

        // What comes next, turned and moved, lands back on the first.
        let mut moved = haemoglobin.clone();
        moved.id = "moved".into();
        for a in &mut moved.atoms {
            let [x, y, z] = a.pos;
            a.pos = [z + 30.0, y, -x];
        }
        view.update(cx, |this, cx| this.finish(Ok(moved), cx));
        cx.run_until_parked();
        let (fit, comparing, armed) = view.update(cx, |this, _| {
            (this.fit, this.compare.is_some(), this.overlay_next)
        });
        let (rmsd, matched) = fit.expect("nothing was laid over");
        assert!(comparing && !armed);
        assert!(rmsd < 0.01, "{rmsd}");
        assert!(matched > 500);
        let p = view.update(cx, |this, _| {
            (
                this.structure.as_ref().unwrap().atoms[10].pos,
                this.compare.as_ref().unwrap().atoms[10].pos,
            )
        });
        assert!((0..3).all(|k| (p.0[k] - p.1[k]).abs() < 0.01));

        click(cx, "comparison-close");
        let (comparing, tint) = view.update(cx, |this, _| {
            (
                this.compare.is_some(),
                this.structure.as_ref().unwrap().tint,
            )
        });
        assert!(!comparing && tint.is_none());
    }

    #[gpui::test]
    fn a_session_opens_as_it_was_left(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        view.update(cx, |this, cx| {
            let s = formats::load_bytes(SAMPLES[5].file, SAMPLES[5].bytes.to_vec()).unwrap();
            this.show(s, cx);
            this.repr = Repr::BallAndStick;
            this.coloring = Coloring::BFactor;
            this.palette = Palette::Viridis;
            this.camera.rotate(0.7, -0.3);
            this.camera.zoom = 2.5;
            Arc::make_mut(this.structure.as_mut().unwrap()).chains[1].color = Some([0.1, 0.2, 0.3]);
            this.measurements.push(Measurement {
                kind: measure::Kind::Distance,
                which: 0,
                atoms: vec![10, 400],
            });
            this.selected = Some(Target {
                which: 0,
                picked: Picked {
                    atom: 55,
                    residue: true,
                    depth: 0.0,
                },
            });
            this.rebuild(cx);
            this.next_site("HEM", cx);
        });
        cx.executor().advance_clock(Duration::from_secs(1));
        cx.run_until_parked();
        let json = view.update(cx, |this, _| this.session_json()).unwrap();
        let before = view.update(cx, |this, _| {
            (
                this.camera.orientation,
                this.camera.zoom,
                this.site.as_ref().map(|s| s.ligand),
                this.structure.as_ref().unwrap().atoms.len(),
            )
        });

        let (other, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        other
            .update(cx, |this, cx| this.open_session(&json, cx))
            .unwrap();
        cx.run_until_parked();
        other.update(cx, |this, _| {
            let s = this.structure.as_ref().unwrap();
            assert_eq!(s.atoms.len(), before.3);
            assert_eq!(s.source, Source::File);
            assert_eq!(this.repr, Repr::BallAndStick);
            assert_eq!(this.coloring, Coloring::BFactor);
            assert_eq!(this.palette, Palette::Viridis);
            assert_eq!(s.chains[1].color, Some([0.1, 0.2, 0.3]));
            assert!((this.camera.zoom - before.1).abs() < 1e-5);
            assert!((0..3)
                .all(|i| (0..3)
                    .all(|j| (this.camera.orientation[i][j] - before.0[i][j]).abs() < 1e-5)));
            assert_eq!(this.measurements.len(), 1);
            assert_eq!(this.measurements[0].atoms, vec![10, 400]);
            assert_eq!(this.selected.map(|t| t.picked.atom), Some(55));
            assert_eq!(this.site.as_ref().map(|s| s.ligand), before.2);
        });
        // A session from a newer profold is refused, not misread.
        let newer = json.replacen("\"profold\": 1", "\"profold\": 99", 1);
        assert!(other
            .update(cx, |this, cx| this.open_session(&newer, cx))
            .is_err());
    }

    #[gpui::test]
    fn water_is_kept_but_hidden_beside_a_protein(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        view.update(cx, |this, cx| {
            let s = formats::load_bytes(SAMPLES[5].file, SAMPLES[5].bytes.to_vec()).unwrap();
            this.show(s, cx);
        });
        cx.run_until_parked();
        let (waters, hidden) = view.update(cx, |this, _| {
            let s = this.structure.as_ref().unwrap();
            (
                s.residues.iter().filter(|r| r.name == "HOH").count(),
                this.appearance.hide_solvent,
            )
        });
        assert!(waters > 100, "haemoglobin's waters were dropped");
        assert!(hidden);
        let spheres = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| this.base.spheres.len())
        };
        let without = spheres(&view, cx);
        click(cx, "show-water");
        assert!(view.update(cx, |this, _| !this.appearance.hide_solvent));
        assert_eq!(
            spheres(&view, cx),
            without + waters,
            "each water is one oxygen drawn"
        );
        click(cx, "show-water");
        assert_eq!(spheres(&view, cx), without);
    }

    #[gpui::test]
    fn the_camera_spins_and_rocks_until_touched(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        let angle = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| this.camera.angles()[1])
        };
        let start = angle(&view, cx);
        click(cx, "tab-camera");
        click(cx, "spin");
        assert_eq!(view.update(cx, |this, _| this.turning), Some(Turn::Spin));
        cx.executor().advance_clock(Duration::from_millis(700));
        cx.run_until_parked();
        let turned = angle(&view, cx);
        assert!((turned - start).abs() > 5.0, "{start} to {turned}");
        // A press on the view stops it.
        let middle = view.update(cx, |this, _| this.cache.borrow().bounds.center());
        cx.simulate_mouse_down(middle, MouseButton::Left, Modifiers::none());
        cx.simulate_mouse_up(middle, MouseButton::Left, Modifiers::none());
        cx.run_until_parked();
        assert_eq!(view.update(cx, |this, _| this.turning), None);
        let still = angle(&view, cx);
        cx.executor().advance_clock(Duration::from_millis(500));
        cx.run_until_parked();
        assert_eq!(angle(&view, cx), still);

        click(cx, "tab-camera");
        click(cx, "rock");
        cx.executor().advance_clock(Duration::from_millis(1500));
        cx.run_until_parked();
        assert_ne!(angle(&view, cx), still);
        // The same again stops it.
        click(cx, "tab-camera");
        click(cx, "rock");
        assert_eq!(view.update(cx, |this, _| this.turning), None);

        click(cx, "tab-image");
        click(cx, "transparent");
        assert!(view.update(cx, |this, _| this.transparent));
        // Ray tracing is slow, so nothing is traced until it is asked for.
        assert!(!view.update(cx, |this, _| this.traced));
        assert!(!crate::movie::Options::default().traced);
        click(cx, "tab-image");
        click(cx, "traced");
        assert!(view.update(cx, |this, _| this.traced));
    }

    #[gpui::test]
    fn the_slab_moves_with_option_scroll_and_shortcuts_show(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        click(cx, "tab-camera");
        click(cx, "perspective");
        assert!(view.update(cx, |this, _| this.camera.perspective));
        // A double-click reframes, and keeps the lens.
        let middle = view.update(cx, |this, _| this.cache.borrow().bounds.center());
        cx.simulate_event(MouseDownEvent {
            position: middle,
            button: MouseButton::Left,
            modifiers: Modifiers::none(),
            click_count: 2,
            first_mouse: false,
        });
        cx.simulate_mouse_up(middle, MouseButton::Left, Modifiers::none());
        cx.run_until_parked();
        assert!(view.update(cx, |this, _| this.camera.perspective));

        click(cx, "tab-camera");
        click(cx, "slab");
        let slab = view
            .update(cx, |this, _| this.camera.slab)
            .expect("no slab");
        let scroll = |cx: &mut VisualTestContext, shift: bool| {
            cx.simulate_event(ScrollWheelEvent {
                position: middle,
                delta: ScrollDelta::Pixels(gpui::point(px(0.), px(40.))),
                modifiers: Modifiers {
                    alt: true,
                    shift,
                    ..Modifiers::none()
                },
                ..Default::default()
            });
            cx.run_until_parked();
        };
        let zoom = view.update(cx, |this, _| this.camera.zoom);
        scroll(cx, false);
        let moved = view.update(cx, |this, _| this.camera.slab).unwrap();
        assert!(
            moved[0] > slab[0] && moved[1] > slab[1],
            "{slab:?} to {moved:?}"
        );
        assert!(
            ((moved[0] - moved[1]) - (slab[0] - slab[1])).abs() < 1e-4,
            "it changed thickness"
        );
        assert_eq!(
            view.update(cx, |this, _| this.camera.zoom),
            zoom,
            "it zoomed as well"
        );
        scroll(cx, true);
        let thick = view.update(cx, |this, _| this.camera.slab).unwrap();
        assert!(thick[0] - thick[1] > moved[0] - moved[1]);
        click(cx, "tab-camera");
        click(cx, "slab");
        assert!(view.update(cx, |this, _| this.camera.slab.is_none()));

        view.update(cx, |this, cx| this.toggle_shortcuts(cx));
        cx.run_until_parked();
        click(cx, "shortcuts");
        assert!(view.update(cx, |this, _| !this.shortcuts));
    }

    #[gpui::test]
    fn a_selection_is_drawn_its_own_way(cx: &mut TestAppContext) {
        cx.update(crate::ui::input::bind_keys);
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        view.update(cx, |this, cx| {
            let s = formats::load_bytes(SAMPLES[5].file, SAMPLES[5].bytes.to_vec()).unwrap();
            this.show(s, cx);
        });
        cx.run_until_parked();
        let (tris, capsules) = view.update(cx, |this, _| {
            (this.base.tris.len(), this.base.capsules.len())
        });

        // A bad one says what is wrong, and adds nothing.
        assert!(!view.update(cx, |this, cx| this.add_layer("resn HIS and", cx)));
        assert!(view.update(cx, |this, _| this.select_error.is_some()
            && this.layers.is_empty()));

        // Typed into its field: the histidines by the irons, as sticks.
        let field = cx.debug_bounds("select-field").expect("no selection field");
        cx.simulate_click(field.center(), Modifiers::none());
        cx.simulate_input("byres (resn HIS and within 3 of elem Fe)");
        cx.simulate_keystrokes("enter");
        cx.run_until_parked();
        let (n, more) = view.update(cx, |this, _| {
            (
                this.layers
                    .first()
                    .map(|l| l.atoms.iter().filter(|a| **a).count()),
                this.base.capsules.len() > capsules,
            )
        });
        assert_eq!(n, Some(40), "four histidines of ten atoms");
        assert!(more, "the sticks are not drawn");
        assert!(view.update(cx, |this, cx| this.select_input.read(cx).text().is_empty()));

        // Hidden, chain A's ribbon goes.
        assert!(view.update(cx, |this, cx| this.add_layer("chain A", cx)));
        click(cx, "layer-1-Hide");
        let (hidden, now) = view.update(cx, |this, _| {
            (this.appearance.hidden.is_some(), this.base.tris.len())
        });
        assert!(hidden && now < tris * 4 / 5, "{now} of {tris}");
        click(cx, "layer-1-color");
        assert_eq!(
            view.update(cx, |this, _| this.layers[1].color),
            LAYER_COLORS[1]
        );
        click(cx, "layer-1-remove");
        assert!(view.update(cx, |this, _| this.appearance.hidden.is_none()
            && this.base.tris.len() == tris));

        // And it is kept in a session.
        let json = view.update(cx, |this, _| this.session_json()).unwrap();
        let (other, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        other
            .update(cx, |this, cx| this.open_session(&json, cx))
            .unwrap();
        let kept = other.update(cx, |this, _| {
            this.layers
                .iter()
                .map(|l| {
                    (
                        l.selection.text().to_string(),
                        l.atoms.iter().filter(|a| **a).count(),
                    )
                })
                .collect::<Vec<_>>()
        });
        assert_eq!(
            kept,
            vec![("byres (resn HIS and within 3 of elem Fe)".to_string(), 40)]
        );
    }

    #[gpui::test]
    fn the_plots_point_at_residues(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        view.update(cx, |this, cx| {
            let mut s = formats::load_bytes(SAMPLES[2].file, SAMPLES[2].bytes.to_vec()).unwrap();
            // A made-up error, small along the diagonal.
            let n = s.residues.iter().filter(|r| !r.outside_polymer()).count();
            let values = (0..n * n)
                .map(|k| ((k / n) as f32 - (k % n) as f32).abs())
                .collect();
            s.pae = Some(Arc::new(crate::protein::Pae {
                n,
                values,
                max: 30.0,
            }));
            this.show(s, cx);
        });
        // Tall enough for both plots to be in view in the panel.
        cx.simulate_resize(gpui::size(px(1400.), px(2600.)));
        cx.run_until_parked();
        // Both start closed, under rows of their own.
        assert!(cx.debug_bounds("ramachandran").is_none() && cx.debug_bounds("pae").is_none());
        for row in ["open-ramachandran", "open-pae"] {
            let at = cx
                .debug_bounds(row)
                .expect("no row to open a plot by")
                .center();
            cx.simulate_click(at, Modifiers::none());
            cx.run_until_parked();
        }
        let plot = cx
            .debug_bounds("ramachandran")
            .expect("no Ramachandran plot");
        // A dot to point at: residue 10's, where the plot puts it.
        let (r, phi, psi) = view.update(cx, |this, _| {
            let t = this.structure.as_ref().unwrap().torsions();
            let r = (5..40).find(|&r| t[r].is_some()).unwrap();
            (r, t[r].unwrap().0, t[r].unwrap().1)
        });
        let at = plot.origin
            + gpui::point(
                plot.size.width * ((phi + 180.0) / 360.0),
                plot.size.height * ((180.0 - psi) / 360.0),
            );
        cx.simulate_mouse_move(at, None, Modifiers::none());
        cx.run_until_parked();
        let pointed = view.update(cx, |this, _| this.rama_hover);
        assert!(pointed.is_some(), "no dot under the pointer");
        assert_eq!(view.update(cx, |this, _| this.seq_hover), pointed);
        cx.simulate_click(at, Modifiers::none());
        cx.run_until_parked();
        assert_eq!(view.update(cx, |this, _| this.selected_residue()), pointed);
        let _ = r;

        let pae = cx.debug_bounds("pae").expect("no PAE plot");
        let quarter = pae.origin + gpui::point(pae.size.width * 0.25, pae.size.height * 0.75);
        cx.simulate_mouse_move(quarter, None, Modifiers::none());
        cx.run_until_parked();
        let (i, j) = view
            .update(cx, |this, _| this.pae_hover)
            .expect("nothing pointed at in the PAE");
        assert!(i > j, "rows run down, columns across: {i} {j}");
    }

    #[gpui::test]
    fn the_film_dialog_sets_up_a_turning_trajectory(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        let melt = SAMPLES
            .iter()
            .find(|s| s.file.ends_with("lammpstrj.gz"))
            .unwrap();
        view.update(cx, |this, cx| {
            this.show(
                formats::load_bytes(melt.file, melt.bytes.to_vec()).unwrap(),
                cx,
            );
            this.open_movie(cx);
        });
        cx.run_until_parked();
        assert!(cx.debug_bounds("movie").is_some(), "the dialog is not up");
        let frames = view.update(cx, |this, _| this.frame_count());
        // Turning, every other frame, smoothed.
        view.update(cx, |this, cx| {
            this.change(cx, |o| {
                o.motion = crate::movie::Motion::Turn;
                o.every = 2;
                o.smooth = 1;
                o.trajectory = true;
                o.seconds = 0.0;
            })
        });
        cx.run_until_parked();
        let (count, _) = view.update(cx, |this, _| this.movie_length(&this.options()));
        let used = frames.div_ceil(2);
        assert_eq!(count, (used - 1) * 2 + 1);

        // Recording a path makes it the camera's motion.
        view.update(cx, |this, cx| {
            this.close_movie(cx);
            this.start_recording(cx);
        });
        assert!(view.update(cx, |this, _| !this.movie.open && this.recording()));
        view.update(cx, |this, cx| {
            this.camera.rotate(0.5, 0.0);
            this.stop_recording(cx);
        });
        cx.run_until_parked();
        let motion = view.update(cx, |this, _| this.options().motion);
        assert_eq!(motion, crate::movie::Motion::Recorded);
        assert!(
            view.update(cx, |this, _| this.movie.open && !this.recording()),
            "the dialog did not come back"
        );
    }

    #[gpui::test]
    fn a_trajectory_moves_the_structure_shown(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        let crambin = formats::load_bytes(SAMPLES[2].file, SAMPLES[2].bytes.to_vec()).unwrap();
        view.update(cx, |this, cx| this.show(crambin.clone(), cx));
        cx.run_until_parked();
        // Three frames, each shifted a little further along x.
        let frames: Vec<crate::protein::Snapshot> = (0..3)
            .map(|k| crate::protein::Snapshot {
                step: None,
                pos: crambin
                    .atoms
                    .iter()
                    .map(|a| [a.pos[0] + k as f32, a.pos[1], a.pos[2]])
                    .collect(),
                cell: None,
            })
            .collect();
        let traj = Structure::beads("run.dcd", crambin.atoms.len(), frames);
        view.update(cx, |this, cx| this.finish(Ok(traj), cx));
        cx.run_until_parked();
        let (count, name, first) = view.update(cx, |this, _| {
            let s = this.structure.as_ref().unwrap();
            (
                s.frame_count(),
                s.atoms[1].name().to_string(),
                s.atoms[0].pos,
            )
        });
        assert_eq!(count, 3);
        assert_eq!(
            name,
            crambin.atoms[1].name(),
            "the atoms are no longer the structure's"
        );
        assert_eq!(first, crambin.atoms[0].pos);
        view.update(cx, |this, cx| this.set_frame(2, cx));
        let moved = view.update(cx, |this, _| this.structure.as_ref().unwrap().atoms[0].pos);
        assert!((moved[0] - crambin.atoms[0].pos[0] - 2.0).abs() < 1e-5);

        // One that does not fit opens on its own, as beads.
        let other = Structure::beads(
            "other.xtc",
            10,
            vec![crate::protein::Snapshot {
                step: None,
                pos: vec![[0.0; 3]; 10],
                cell: None,
            }],
        );
        view.update(cx, |this, cx| this.finish(Ok(other), cx));
        assert_eq!(
            view.update(cx, |this, _| this.structure.as_ref().unwrap().atoms.len()),
            10
        );
    }

    #[test]
    fn a_trajectory_turning_is_written_as_a_looping_gif() {
        use image::AnimationDecoder;
        let melt = SAMPLES
            .iter()
            .find(|s| s.file.ends_with("lammpstrj.gz"))
            .unwrap();
        let s = formats::load_bytes(melt.file, melt.bytes.to_vec()).unwrap();
        assert!(s.frame_count() > 4);
        let options = crate::movie::Options {
            motion: crate::movie::Motion::Turn,
            last: 4,
            seconds: 0.0,
            size: crate::movie::Size::Custom,
            width: 120,
            height: 90,
            supersample: 2,
            format: crate::movie::Format::Gif,
            counter: true,
            scale_bar: true,
            caption: "melt".into(),
            ..Default::default()
        };
        let camera = Camera::frame(&s);
        let shots = crate::movie::plan(&options, &camera, 0, s.frame_count(), None);
        assert_eq!(shots.len(), 4);
        let job = crate::movie::Job {
            structure: s,
            build: Arc::new(|s: &Structure| {
                render::build_scene(
                    s,
                    Repr::Spacefill,
                    Coloring::Element,
                    1.0,
                    Style::Soft,
                    Palette::Spectrum,
                    true,
                    false,
                    &Appearance::default(),
                )
            }),
            style: Style::Soft,
            page: [1.0; 3],
            view: (120, 90),
            shots,
            options,
            legend: None,
            labels: Vec::new(),
        };
        let dir = std::env::temp_dir().join(format!("profold-gif-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let progress = crate::movie::Progress::default();
        let path = crate::movie::export(&job, &dir.join("melt.gif"), &progress).unwrap();
        assert_eq!(progress.done.load(std::sync::atomic::Ordering::Relaxed), 4);
        let bytes = std::fs::read(&path).unwrap();
        let decoded = image::codecs::gif::GifDecoder::new(std::io::Cursor::new(bytes))
            .unwrap()
            .into_frames()
            .collect_frames()
            .unwrap();
        assert_eq!(decoded.len(), 4);
        assert_eq!(decoded[0].buffer().dimensions(), (120, 90));
        assert_ne!(
            decoded[0].buffer().as_raw(),
            decoded[2].buffer().as_raw(),
            "every frame the same"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[gpui::test]
    fn undo_steps_back_through_how_it_was_shown(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        // Each change after a pause, so none is folded into the last.
        let pause = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| {
                this.last_change = std::time::Instant::now() - UNDO_PAUSE * 2
            })
        };
        let state = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| {
                (this.repr, this.coloring, this.measurements.len())
            })
        };
        let start = state(&view, cx);
        pause(&view, cx);
        view.update(cx, |this, cx| this.set_repr(Repr::Spacefill, cx));
        pause(&view, cx);
        view.update(cx, |this, cx| {
            this.coloring = Coloring::Rainbow;
            this.rebuild(cx);
        });
        pause(&view, cx);
        view.update(cx, |this, cx| {
            this.measurements.push(Measurement {
                kind: measure::Kind::Distance,
                which: 0,
                atoms: vec![0, 1],
            });
            this.compose(cx);
        });
        assert_eq!(state(&view, cx), (Repr::Spacefill, Coloring::Rainbow, 1));
        view.update(cx, |this, cx| this.undo(false, cx));
        assert_eq!(state(&view, cx), (Repr::Spacefill, Coloring::Rainbow, 0));
        view.update(cx, |this, cx| this.undo(false, cx));
        assert_eq!(state(&view, cx), (Repr::Spacefill, start.1, 0));
        view.update(cx, |this, cx| this.undo(false, cx));
        assert_eq!(state(&view, cx), start);
        // Nothing before the start.
        view.update(cx, |this, cx| this.undo(false, cx));
        assert_eq!(state(&view, cx), start);
        view.update(cx, |this, cx| this.undo(true, cx));
        view.update(cx, |this, cx| this.undo(true, cx));
        assert_eq!(state(&view, cx), (Repr::Spacefill, Coloring::Rainbow, 0));
        // A new change drops what could be redone.
        pause(&view, cx);
        view.update(cx, |this, cx| this.set_style(Style::Studio, cx));
        view.update(cx, |this, cx| this.undo(true, cx));
        assert_eq!(view.update(cx, |this, _| this.style), Style::Studio);
        // A slider's drag, many changes at once, is one step.
        pause(&view, cx);
        for size in [1.1, 1.3, 1.6] {
            view.update(cx, |this, cx| {
                this.appearance.set(
                    Species::Element(6),
                    Look {
                        size,
                        ..Look::default()
                    },
                );
                this.rebuild(cx);
            });
        }
        view.update(cx, |this, cx| this.undo(false, cx));
        assert_eq!(
            view.update(cx, |this, _| this.appearance.get(&Species::Element(6)).size),
            1.0
        );
    }

    #[gpui::test]
    fn a_viewpoint_saved_is_flown_back_to_and_recent_files_are_listed(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        view.update(cx, |this, cx| {
            this.camera.rotate(0.8, 0.3);
            this.camera.zoom = 2.0;
            cx.notify();
        });
        click(cx, "tab-camera");
        click(cx, "save-view");
        let saved = view.update(cx, |this, _| this.camera);
        view.update(cx, |this, cx| {
            this.camera.rotate(-1.5, 0.9);
            this.camera.zoom = 0.7;
            cx.notify();
        });
        click(cx, "tab-camera");
        click(cx, "view-0");
        cx.executor().advance_clock(Duration::from_secs(1));
        cx.run_until_parked();
        let now = view.update(cx, |this, _| this.camera);
        assert!((now.zoom - saved.zoom).abs() < 1e-4);
        assert!((0..3).all(
            |i| (0..3).all(|j| (now.orientation[i][j] - saved.orientation[i][j]).abs() < 1e-4)
        ));

        // What was opened is offered again.
        view.update(cx, |this, cx| {
            this.opening = Some("sample:1ubq.pdb.gz".into());
            let s = formats::load_bytes(SAMPLES[3].file, SAMPLES[3].bytes.to_vec()).unwrap();
            this.finish(Ok(s), cx);
        });
        assert_eq!(
            crate::settings::recent().first().map(String::as_str),
            Some("sample:1ubq.pdb.gz")
        );
        click(cx, "tab-samples");
        assert!(
            cx.debug_bounds("recent-0").is_some(),
            "no recent entry in the menu"
        );
    }

    #[gpui::test]
    fn menus_in_a_phone_shaped_window_stay_inside_it(cx: &mut TestAppContext) {
        let (_view, cx) = open(cx);
        theme::set_vertical(true);
        let width = px(380.);
        cx.simulate_resize(gpui::size(width, px(820.)));
        cx.run_until_parked();
        let bounds = |cx: &mut VisualTestContext, name: &str| {
            let selector: &'static str = Box::leak(name.to_string().into_boxed_str());
            cx.debug_bounds(selector)
                .unwrap_or_else(|| panic!("{name} is not on screen"))
        };
        // The card within the window, and each item within the card.
        let check = |cx: &mut VisualTestContext, card: &str, items: Vec<String>| {
            let c = bounds(cx, card);
            assert!(
                c.origin.x >= px(0.) && c.right() <= width,
                "{card} runs off the window: {c:?}"
            );
            for name in items {
                let b = bounds(cx, &name);
                assert!(
                    b.origin.x >= c.origin.x && b.right() <= c.right(),
                    "{name} runs out of {card}: {b:?} in {c:?}"
                );
            }
        };
        click(cx, "palette-button");
        let names = Palette::ALL
            .iter()
            .map(|p| format!("palette-{}", p.label()))
            .collect();
        check(cx, "palette-menu", names);
        click(cx, "palette-button");
        click(cx, "colour-menu");
        let names = vec![format!("colour-item-{}", Coloring::Element.label())];
        check(cx, "colour-menu-list", names);
        // At its own width, hung from its button, not spread across the row.
        let list = bounds(cx, "colour-menu-list");
        let button = bounds(cx, "colour-menu");
        assert_eq!(list.size.width, px(156.));
        assert!(
            list.bottom() <= button.top(),
            "the colour menu is not over its button: {list:?} {button:?}"
        );
        click(cx, "colour-menu");
        click(cx, "repr-menu");
        let names = vec!["repr-item-Spacefill".to_string()];
        check(cx, "repr-menu-list", names);
        theme::set_vertical(false);
    }

    #[gpui::test]
    fn a_phone_shaped_window_folds_the_rows_into_menus(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        theme::set_vertical(true);
        cx.simulate_resize(gpui::size(px(400.), px(820.)));
        cx.run_until_parked();
        assert!(
            view.update(cx, |this, _| this.compact),
            "no menus in a narrow window"
        );
        click(cx, "repr-menu");
        assert!(view.update(cx, |this, _| this.open == Some(Dropdown::Repr)));
        click(cx, "repr-item-Spacefill");
        assert!(view.update(cx, |this, _| this.repr == Repr::Spacefill
            && this.open.is_none()));

        cx.simulate_resize(gpui::size(px(1180.), px(760.)));
        cx.run_until_parked();
        assert!(
            !view.update(cx, |this, _| this.compact),
            "a wide window keeps its rows"
        );
        theme::set_vertical(false);
    }

    #[gpui::test]
    fn palette_menu_recolours_elements(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        click(cx, "element-colours-palette");
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
        let (orientation, zoom) =
            view.update(cx, |this, _| (this.camera.orientation, this.camera.zoom));
        assert_eq!(orientation, Camera::home(), "not back at the opening view");
        assert_eq!(zoom, 1.0);
    }

    #[gpui::test]
    fn standard_element_colours_ignore_the_palette(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        view.update(cx, |this, cx| {
            this.standard_colors = true;
            this.rebuild(cx);
        });
        // The elements' own colours, before the style and theme tone them,
        // which another test's dark mode would otherwise change under us.
        let elements = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| {
                let s = this.structure.as_ref().unwrap();
                render::element_colors(s, this.palette, this.standard_colors)
            })
        };
        let standard = elements(&view, cx);

        // The palette stays open to chains and ramps, but leaves elements be.
        click(cx, "palette-button");
        assert!(view.update(cx, |this, _| this.open == Some(Dropdown::Palette)));
        click(cx, "palette-Magma");
        assert_eq!(elements(&view, cx), standard);

        click(cx, "element-colours-palette");
        assert!(view.update(cx, |this, _| !this.standard_colors));
        assert_ne!(elements(&view, cx), standard);
        click(cx, "element-colours-standard");
        assert_eq!(elements(&view, cx), standard);
    }

    #[gpui::test]
    fn picking_a_colour_paints_the_element(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        click(cx, "dot-C");
        assert!(
            view.update(cx, |this, _| this.picking.is_some()),
            "the picker did not open"
        );
        click(cx, "swatch-1");
        let carbon = view.update(cx, |this, _| {
            this.appearance.get(&Species::Element(6)).color
        });
        assert_eq!(carbon, Some(SWATCHES[0]));
        assert!(
            colors(&view, cx).contains(&SWATCHES[0]) || {
                // The inked styles tone colours; the scene holds the toned one.
                view.update(cx, |this, _| this.style.inked())
            }
        );
    }

    #[gpui::test]
    fn picking_a_colour_paints_the_chain(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        view.update(cx, |this, cx| {
            let bytes = include_bytes!("../../samples/4hhb.pdb.gz").to_vec();
            this.show(formats::load_bytes("4hhb.pdb.gz", bytes).unwrap(), cx);
            this.coloring = Coloring::Rainbow;
        });
        cx.run_until_parked();
        click(cx, "chain-dot-1");
        assert!(view.update(cx, |this, _| this.picking == Some(Pick::Chain(1))));
        click(cx, "swatch-1");
        let (coloring, chain) = view.update(cx, |this, _| {
            let s = this.structure.as_ref().unwrap();
            (this.coloring, s.chains[1].color)
        });
        // A chain's colour shows only when colouring by chain.
        assert_eq!(coloring, Coloring::Chain);
        assert_eq!(chain, Some(SWATCHES[0]));

        view.update(cx, |this, cx| {
            this.picking = None;
            cx.notify();
        });
        cx.run_until_parked();
        click(cx, "reset-chains");
        let chain = view.update(cx, |this, _| {
            this.structure.as_ref().unwrap().chains[1].color
        });
        assert_eq!(chain, None);
    }

    #[test]
    fn the_palette_names_what_keeps_its_colour() {
        let names = |v: &[&str]| v.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        assert_eq!(kept_note(&[]), None);
        assert_eq!(kept_note(&names(&["C"])).unwrap(), "C keeps its own colour");
        assert_eq!(
            kept_note(&names(&["Chain B", "C", "O"])).unwrap(),
            "Chain B, C and O keep their own colours"
        );
        assert_eq!(
            kept_note(&names(&["A", "B", "C", "D", "E"])).unwrap(),
            "A, B, C and 2 more keep their own colours"
        );
    }

    #[gpui::test]
    fn a_palette_change_keeps_picked_colours_until_reset(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        view.update(cx, |this, cx| {
            this.set_color(Species::Element(6), Some(SWATCHES[0]), cx)
        });
        click(cx, "palette-button");
        click(cx, "palette-Viridis");
        let carbon = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| {
                this.appearance.get(&Species::Element(6)).color
            })
        };
        assert_eq!(carbon(&view, cx), Some(SWATCHES[0]));
        // Said straight away, with a way back.
        assert!(view.update(cx, |this, _| this.kept.is_some()));
        click(cx, "kept-reset");
        assert_eq!(carbon(&view, cx), None);
        assert!(view.update(cx, |this, _| this.kept.is_none()));

        // A palette change with nothing picked says nothing.
        click(cx, "palette-button");
        click(cx, "palette-Magma");
        assert!(view.update(cx, |this, _| this.kept.is_none()));

        // The menu's own Reset, and a row's own.
        view.update(cx, |this, cx| {
            this.set_color(Species::Element(6), Some(SWATCHES[0]), cx)
        });
        click(cx, "palette-button");
        click(cx, "reset-colours");
        assert_eq!(carbon(&view, cx), None);
        click(cx, "palette-button");
        view.update(cx, |this, cx| {
            this.set_color(Species::Element(6), Some(SWATCHES[0]), cx)
        });
        cx.run_until_parked();
        let row = view.update(cx, |this, _| {
            let s = this.structure.as_ref().unwrap();
            s.composition().iter().position(|(e, _)| e.0 == 6).unwrap()
        });
        let selector: &'static str = format!("species-unpick-{row}").leak();
        click(cx, selector);
        assert_eq!(carbon(&view, cx), None);
    }

    #[gpui::test]
    fn a_ligand_takes_a_colour_on_its_carbons_or_whole(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        view.update(cx, |this, cx| {
            let bytes = include_bytes!("../../samples/4hhb.pdb.gz").to_vec();
            this.show(formats::load_bytes("4hhb.pdb.gz", bytes).unwrap(), cx);
        });
        cx.run_until_parked();
        let haem = |view: &Entity<Viewer>, cx: &mut VisualTestContext| {
            view.update(cx, |this, _| this.ligand_paint("HEM"))
        };

        // A colour picked goes on the carbons first, and leaves the scheme be.
        click(cx, "ligand-dot-HEM");
        click(cx, "swatch-1");
        assert_eq!(haem(&view, cx), Some(LigandPaint::Carbons(SWATCHES[0])));
        assert_eq!(view.update(cx, |this, _| this.coloring), Coloring::Chain);
        view.update(cx, |this, cx| {
            this.picking = None;
            cx.notify();
        });
        cx.run_until_parked();

        click(cx, "ligand-HEM-Solid");
        assert_eq!(haem(&view, cx), Some(LigandPaint::Solid(SWATCHES[0])));
        click(cx, "ligand-HEM-Element");
        assert_eq!(haem(&view, cx), None);
        // Straight to Solid with no colour yet: gold.
        click(cx, "ligand-HEM-Solid");
        assert_eq!(haem(&view, cx), Some(LigandPaint::Solid(LIGAND_GOLD)));
        click(cx, "ligand-unpick-0");
        assert_eq!(haem(&view, cx), None);

        // Phosphate has no carbons, so it goes solid.
        assert!(cx.debug_bounds("ligand-PO4-Carbons").is_none());
        view.update(cx, |this, cx| {
            this.paint(Pick::Ligand("PO4".into()), Some(SWATCHES[0]), cx)
        });
        let phosphate = view.update(cx, |this, _| this.ligand_paint("PO4"));
        assert_eq!(phosphate, Some(LigandPaint::Solid(SWATCHES[0])));
    }

    fn hit(kind: Kind, id: &str, title: &str) -> Hit {
        Hit::for_test(kind, id, title)
    }

    #[gpui::test]
    fn results_list_opens_the_chosen_hit(cx: &mut TestAppContext) {
        cx.update(crate::ui::input::bind_keys);
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        cx.simulate_input("green fluorescent protein");
        let (detected, searching) = view.update(cx, |this, _| {
            (this.detected.clone(), this.searching.clone())
        });
        assert_eq!(detected, Input::Query("green fluorescent protein".into()));
        assert_eq!(
            searching.as_deref(),
            Some("green fluorescent protein"),
            "no search was started"
        );

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
        assert!(
            cx.debug_bounds("hit-2").is_some(),
            "the results list is not on screen"
        );

        // Down twice, up once: the AlphaFold model.
        cx.simulate_keystrokes("down down up");
        assert_eq!(view.update(cx, |this, _| this.highlight), 1);
        cx.simulate_keystrokes("escape");
        assert!(
            !view.update(cx, |this, _| this.dropdown),
            "escape left the list up"
        );
        cx.simulate_keystrokes("down");
        assert!(
            view.update(cx, |this, _| this.dropdown),
            "down did not bring the list back"
        );

        click(cx, "hit-0");
        let busy = view.update(cx, |this, _| {
            (this.dropdown, this.busy.is_some(), this.log.last().cloned())
        });
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
        assert!(
            cx.debug_bounds("source-link").is_some(),
            "no link to PubChem"
        );
        let page = view.update(cx, |this, _| this.structure.as_ref().unwrap().page());
        assert_eq!(
            page.as_deref(),
            Some("https://pubchem.ncbi.nlm.nih.gov/compound/2519")
        );
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
        let (busy, last) = view.update(cx, |this, _| {
            (this.busy.is_some(), this.log.last().cloned())
        });
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
        let (detected, searching) = view.update(cx, |this, _| {
            (this.detected.clone(), this.searching.clone())
        });
        assert!(matches!(detected, Input::Sequence(_)));
        assert!(searching.is_none(), "a fold tab searched");
    }

    #[gpui::test]
    fn a_viewer_never_offers_to_fold(cx: &mut TestAppContext) {
        let (view, cx) = cx.add_window_view(|window, cx| Viewer::new(Mode::View, window, cx));
        cx.simulate_input("1CRN");
        assert_eq!(
            view.update(cx, |this, _| this.detected.clone()),
            Input::Pdb("1CRN".into())
        );
        assert!(cx.debug_bounds("compare").is_none());
        assert!(cx.debug_bounds("model").is_none());
    }
}
