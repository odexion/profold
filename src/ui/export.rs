//! The film dialog: what moves, which frames, how large and on what, in
//! which format and with what written over it; a preview, the export with
//! its progress, and recording a camera path by hand to play back.

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use gpui::{div, prelude::*, px, relative, Context, MouseButton, SharedString, Task};

use super::{pill_keyed, spaced, swallow, swatch, Field, Pick, Viewer};
use crate::movie::{
    self, Axis, Background, Corner, Format, Legend, Motion, Options, Place, Sample, Size, TextSize,
    TimeUnit,
};
use crate::protein::Structure;
use crate::render::{self, Coloring, Repr, Scene, Species};
use crate::settings;
use crate::theme;

/// The values in the dialog that are typed.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Knob {
    Degrees,
    First,
    Last,
    Every,
    Fps,
    Seconds,
    Width,
    Height,
    Quality,
    Dt,
    Caption,
    LegendTitle,
}

struct Recording {
    started: Instant,
    samples: Vec<Sample>,
    _task: Task<()>,
}

/// What the dialog keeps.
#[derive(Default)]
pub(super) struct Movie {
    options: Option<Options>,
    pub(super) open: bool,
    /// The film as it plays in the dialog: the frame showing and which of
    /// the film's it is, and every frame drawn so far while they all fit in
    /// `PREVIEW_BYTES`, to play again without drawing them again.
    shown: Option<(usize, Arc<gpui::RenderImage>)>,
    kept: Vec<Arc<gpui::RenderImage>>,
    /// The width, in points, the preview was drawn to be shown at.
    drawn_for: f32,
    preview_task: Option<Task<()>>,
    /// The wait after a change before the preview is drawn again, so typing
    /// or clicking through choices does not start one for each.
    refresh_task: Option<Task<()>>,
    preview_cancel: Option<Arc<AtomicBool>>,
    previewing: bool,
    progress: Option<Arc<movie::Progress>>,
    task: Option<Task<()>>,
    recording: Option<Recording>,
    recorded: Option<Arc<Vec<Sample>>>,
}

impl Drop for Movie {
    fn drop(&mut self) {
        // A tab closed partway through: stop drawing for it.
        if let Some(p) = &self.progress {
            p.cancel.store(true, Ordering::Relaxed);
        }
        if let Some(c) = &self.preview_cancel {
            c.store(true, Ordering::Relaxed);
        }
    }
}

fn saved() -> Options {
    settings::file("movie.json")
        .and_then(|p| std::fs::read_to_string(p).ok())
        .and_then(|t| serde_json::from_str(&t).ok())
        .map(|v| Options::from_json(&v))
        .unwrap_or_default()
}

fn keep(o: &Options) {
    let Some(p) = settings::file("movie.json") else {
        return;
    };
    if let Some(dir) = p.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let _ = std::fs::write(p, o.to_json());
}

/// Seconds as minutes and seconds, to a tenth.
fn clock(s: f32) -> String {
    format!("{}:{:04.1}", (s / 60.0) as u32, s % 60.0)
}

/// The preview's width in the dialog, in points: as wide as the viewport
/// leaves room for, between these.
const PREVIEW_WIDTH: f32 = 296.0;
const PREVIEW_WIDEST: f32 = 640.0;

/// The dialog's width besides the preview: its settings, the gap between,
/// its padding, and a margin each side in the viewport.
const DIALOG_REST: f32 = 410.0 + 22.0 + 36.0 + 24.0;

/// What the preview's frames may hold between them, in bytes; a longer
/// film is drawn again each time round instead of kept.
const PREVIEW_BYTES: usize = 192 << 20;

/// A drawn frame as an image for the window.
fn render_image(p: movie::Picture) -> Arc<gpui::RenderImage> {
    let mut bgra = p.rgba;
    for px in bgra.chunks_exact_mut(4) {
        px.swap(0, 2);
    }
    let buffer =
        image::RgbaImage::from_raw(p.width as u32, p.height as u32, bgra).expect("sized to fit");
    Arc::new(gpui::RenderImage::new(smallvec::smallvec![
        image::Frame::new(buffer)
    ]))
}

/// How long settings must rest before the preview is drawn again.
const PREVIEW_PAUSE: Duration = Duration::from_millis(350);

impl Viewer {
    pub(super) fn options(&self) -> Options {
        self.movie.options.clone().unwrap_or_default()
    }

    fn set_options(&mut self, o: Options, cx: &mut Context<Self>) {
        let o = o.clamped();
        if self.movie.options.as_ref() != Some(&o) {
            keep(&o);
            self.movie.options = Some(o);
            // What was previewed is not what would be made now.
            self.refresh_preview(PREVIEW_PAUSE, cx);
        }
        cx.notify();
    }

    /// Draw the preview again once `after` has passed with no further change.
    fn refresh_preview(&mut self, after: Duration, cx: &mut Context<Self>) {
        if !self.movie.open || self.movie.progress.is_some() {
            return;
        }
        self.stop_preview();
        self.movie.previewing = true;
        self.movie.refresh_task = Some(cx.spawn(async move |this, cx| {
            if !after.is_zero() {
                cx.background_executor().timer(after).await;
            }
            let _ = this.update(cx, |this, cx| this.preview_movie(cx));
        }));
    }

    pub(super) fn change(&mut self, cx: &mut Context<Self>, f: impl FnOnce(&mut Options)) {
        let mut o = self.options();
        f(&mut o);
        self.set_options(o, cx);
    }

    /// Open the film dialog.
    pub fn open_movie(&mut self, cx: &mut Context<Self>) {
        if self.structure.is_none() {
            return;
        }
        self.open = None;
        if self.movie.options.is_none() {
            let mut o = saved();
            // A recording made this session is what was wanted.
            if self.movie.recorded.is_none() && o.motion == Motion::Recorded {
                o.motion = Motion::Still;
            }
            self.movie.options = Some(o);
        }
        // Without ffmpeg, a GIF is the film that can be made.
        if let Some(o) = self
            .movie
            .options
            .as_mut()
            .filter(|o| o.format.needs_ffmpeg() && movie::ffmpeg().is_none())
        {
            o.format = Format::Gif;
        }
        self.movie.open = true;
        self.refresh_preview(Duration::ZERO, cx);
        cx.notify();
    }

    pub(super) fn close_movie(&mut self, cx: &mut Context<Self>) {
        self.commit_edit(cx);
        self.movie.open = false;
        if matches!(self.picking, Some(Pick::FilmBackground | Pick::FilmInk)) {
            self.picking = None;
        }
        self.stop_preview();
        cx.notify();
    }

    fn stop_preview(&mut self) {
        if let Some(c) = self.movie.preview_cancel.take() {
            c.store(true, Ordering::Relaxed);
        }
        self.movie.preview_task = None;
        self.movie.refresh_task = None;
        self.movie.previewing = false;
    }

    /// Everything a film of what is on screen needs, drawn as the viewport
    /// draws it, but at full detail.
    fn movie_job(&self, o: &Options) -> Option<movie::Job> {
        let s = self.structure.as_ref()?;
        let (repr, coloring, style, palette, standard, cell) = (
            self.repr,
            self.coloring,
            self.style,
            self.palette,
            self.standard_colors,
            self.show_cell,
        );
        let mut appearance = self.appearance.clone();
        appearance.hidden = None;
        let layers: Vec<_> = self
            .layers
            .iter()
            .map(|l| (l.selection.clone(), l.style, l.color))
            .collect();
        let main = self.shown(0).is_some();
        let other = self.shown(1).map(|o| Structure::clone(o));
        // The volume's surface stands still through a film, built once.
        let iso = self.iso_drawn();
        let build = move |s: &Structure| -> Scene {
            let mut scene = Scene::default();
            if let Some(iso) = iso.as_ref().filter(|_| main) {
                scene.merge((**iso).clone());
            }
            if main {
                // Layers are found anew for each pose, as they are for each
                // frame in the viewport.
                let mut look = appearance.clone();
                let mut hidden = vec![false; s.atoms.len()];
                let found: Vec<_> = layers
                    .iter()
                    .map(|(sel, st, c)| {
                        let atoms = sel.atoms(s);
                        if *st == render::LayerStyle::Hidden {
                            for (h, a) in hidden.iter_mut().zip(&atoms) {
                                *h |= *a;
                            }
                        }
                        (atoms, *st, *c)
                    })
                    .collect();
                look.hidden = hidden.iter().any(|h| *h).then(|| Arc::new(hidden));
                scene.merge(render::build_scene(
                    s, repr, coloring, 1.0, style, palette, standard, cell, &look,
                ));
                for (atoms, st, c) in &found {
                    scene.merge(render::build_layer(
                        s, atoms, *st, *c, style, palette, standard, &look,
                    ));
                }
            }
            if let Some(other) = &other {
                scene.merge(render::build_scene(
                    other,
                    repr,
                    coloring,
                    1.0,
                    style,
                    palette,
                    standard,
                    cell,
                    &appearance,
                ));
            }
            scene
        };
        let recorded = self.movie.recorded.as_deref().map(|v| v.as_slice());
        let shots = movie::plan(o, &self.camera, s.frame, self.frame_count(), recorded);
        Some(movie::Job {
            structure: Structure::clone(s),
            build: Arc::new(build),
            style,
            page: theme::view(),
            view: self.cache.borrow().size,
            shots,
            options: o.clone(),
            legend: self.film_legend(),
            labels: self
                .labels
                .iter()
                .map(|l| (l.atom, l.text.clone()))
                .collect(),
        })
    }

    /// What the colours on screen mean, as the panel's key says it.
    fn film_legend(&self) -> Option<Legend> {
        let s = self.shown(0).or(self.shown(1))?;
        let palette = self.palette;
        let ramp =
            |f: &dyn Fn(f32) -> [f32; 3]| (0..24).map(|i| f(i as f32 / 23.0)).collect::<Vec<_>>();
        // No more than this many swatches; the rest are counted.
        const MOST: usize = 8;
        let swatches = |title: &str, all: Vec<([f32; 3], String)>| {
            let more = all.len().saturating_sub(MOST);
            Legend::Swatches {
                title: title.into(),
                items: all.into_iter().take(MOST).collect(),
                more,
            }
        };
        Some(match self.coloring {
            Coloring::Property => {
                let i = s.prop_color.map_or(0, |p| p.index);
                let (lo, hi) = s.prop_span(i)?;
                Legend::Ramp {
                    title: s.props[i].name.clone(),
                    colors: ramp(&|t| palette.ramp(t)),
                    low: crate::ui::app::figure(lo),
                    high: crate::ui::app::figure(hi),
                }
            }
            Coloring::Hydrophobic => Legend::Ramp {
                title: "Hydrophobicity".into(),
                colors: ramp(&|t| {
                    render::diverging(t * 2.0 - 1.0, render::HYDROPHILIC, render::HYDROPHOBIC)
                }),
                low: "hydrophilic".into(),
                high: "hydrophobic".into(),
            },
            Coloring::Charge => {
                let surface = self.repr == Repr::Surface;
                Legend::Ramp {
                    title: if surface {
                        "Electrostatic potential"
                    } else {
                        "Charge"
                    }
                    .into(),
                    colors: ramp(&|t| {
                        render::diverging(t * 2.0 - 1.0, render::NEGATIVE, render::POSITIVE)
                    }),
                    low: if surface { "−" } else { "acidic" }.into(),
                    high: if surface { "+" } else { "basic" }.into(),
                }
            }
            Coloring::BFactor => {
                let (lo, hi) = s.b_span()?;
                Legend::Ramp {
                    title: "B-factor".into(),
                    colors: ramp(&|t| palette.ramp(t)),
                    low: format!("{lo:.0} Å²"),
                    high: format!("{hi:.0} Å²"),
                }
            }
            Coloring::Rainbow if s.has_polymer() => Legend::Ramp {
                title: "Sequence".into(),
                colors: ramp(&|t| palette.ramp(t)),
                low: "N-terminus".into(),
                high: "C-terminus".into(),
            },
            Coloring::Structure => swatches(
                "Secondary structure",
                render::SS_COLORS
                    .iter()
                    .map(|(_, name, c)| (*c, name.to_string()))
                    .collect(),
            ),
            Coloring::Confidence if s.has_confidence => swatches(
                "pLDDT",
                render::PLDDT_BANDS
                    .iter()
                    .map(|(name, v)| (render::plddt_band_color(*v), name.to_string()))
                    .collect(),
            ),
            Coloring::Chain if s.has_polymer() => swatches(
                "Chains",
                s.chains
                    .iter()
                    .enumerate()
                    .map(|(i, c)| {
                        let name = if c.id.trim().is_empty() {
                            format!("Chain {}", i + 1)
                        } else {
                            format!("Chain {}", c.id)
                        };
                        (s.chain_color(i, palette), name)
                    })
                    .collect(),
            ),
            _ => {
                // Elements and atom types, as the atoms list colours them.
                let colors = render::element_colors(s, palette, self.standard_colors);
                let own = |species: Species, fallback: [f32; 3]| {
                    self.appearance.get(&species).color.unwrap_or(fallback)
                };
                let mut rows: Vec<([f32; 3], String)> = s
                    .composition()
                    .into_iter()
                    .filter(|(e, _)| e.is_known())
                    .map(|(e, _)| {
                        (
                            own(Species::Element(e.0), colors[e.0 as usize]),
                            e.symbol().to_string(),
                        )
                    })
                    .collect();
                let mut by_kind = vec![0usize; s.kinds.len()];
                for a in s.atoms.iter().filter(|a| !a.element.is_known()) {
                    if let Some(n) = by_kind.get_mut(a.kind as usize) {
                        *n += 1;
                    }
                }
                rows.extend(
                    by_kind
                        .iter()
                        .enumerate()
                        .filter(|(_, n)| **n > 0)
                        .map(|(k, _)| {
                            let name = s.kinds[k].clone();
                            (
                                own(Species::Kind(name.clone()), palette.chain(k)),
                                format!("Type {name}"),
                            )
                        }),
                );
                if rows.is_empty() {
                    return None;
                }
                swatches(
                    if by_kind.iter().any(|n| *n > 0) {
                        "Atoms"
                    } else {
                        "Elements"
                    },
                    rows,
                )
            }
        })
    }

    /// How many frames the film would have, and how long it would run.
    pub(super) fn movie_length(&self, o: &Options) -> (usize, f32) {
        let (now, count) = (
            self.structure.as_ref().map(|s| s.frame).unwrap_or(0),
            self.frame_count(),
        );
        let recorded = self.movie.recorded.as_deref().map(|v| v.as_slice());
        let n = movie::plan(o, &self.camera, now, count, recorded).len();
        (n, n as f32 / o.fps.max(1) as f32)
    }

    /// How wide the preview is shown, in points: what the viewport leaves
    /// room for beside the settings.
    fn preview_width(&self) -> f32 {
        let room = f32::from(self.cache.borrow().bounds.size.width);
        if room <= 0.0 {
            PREVIEW_WIDTH
        } else {
            (room - DIALOG_REST).clamp(PREVIEW_WIDTH, PREVIEW_WIDEST)
        }
    }

    /// The screen's pixels to the point, as the viewport is drawn.
    fn screen_scale(&self) -> f32 {
        let c = self.cache.borrow();
        let w = f32::from(c.bounds.size.width);
        if w > 0.0 {
            (c.size.0 as f32 / w).max(1.0)
        } else {
            2.0
        }
    }

    /// Play the film in the dialog, each frame drawn as it will be made at
    /// the size it is shown, at the film's frame rate as far as drawing
    /// keeps up.
    fn preview_movie(&mut self, cx: &mut Context<Self>) {
        let Some(job) = self.movie_job(&self.options()) else {
            return;
        };
        self.stop_preview();
        let cancel = Arc::new(AtomicBool::new(false));
        self.movie.preview_cancel = Some(cancel.clone());
        self.movie.previewing = true;
        self.movie.kept.clear();
        self.movie.drawn_for = self.preview_width();
        let fps = job.options.fps.max(1);
        let across = (self.movie.drawn_for * self.screen_scale()).round() as usize;
        let (pw, ph) = movie::preview_size(&job, across);
        let fits = job.shots.len() * pw * ph * 4 <= PREVIEW_BYTES;
        // Two frames drawn ahead at most, so drawing keeps pace with
        // playing rather than filling memory.
        // On a thread of its own, since it waits on the channel: waiting in
        // the executor's pool would hold up whatever else runs there.
        let (tx, rx) = std::sync::mpsc::sync_channel(2);
        std::thread::spawn(move || {
            movie::play(&job, across, !fits, &cancel, |i, p| {
                tx.send((i, render_image(p))).is_ok()
            });
        });
        self.movie.preview_task = Some(cx.spawn(async move |this, cx| {
            use std::sync::mpsc::TryRecvError;
            let tick = Duration::from_secs_f32(1.0 / fps as f32);
            loop {
                let next = rx.try_recv();
                let going = this.update(cx, |this, cx| {
                    if !this.movie.open {
                        return None;
                    }
                    let m = &mut this.movie;
                    match next {
                        Ok((i, image)) => {
                            m.previewing = false;
                            if fits {
                                m.kept.push(image.clone());
                            }
                            m.shown = Some((i, image));
                        }
                        // Not drawn yet: look again soon.
                        Err(TryRecvError::Empty) => return Some(false),
                        // Every frame drawn and kept: play them round again.
                        Err(TryRecvError::Disconnected) => {
                            if m.kept.len() < 2 {
                                return None;
                            }
                            let i = m.shown.as_ref().map_or(0, |(i, _)| (i + 1) % m.kept.len());
                            m.shown = Some((i, m.kept[i].clone()));
                        }
                    }
                    cx.notify();
                    Some(true)
                });
                match going {
                    Ok(Some(true)) => cx.background_executor().timer(tick).await,
                    Ok(Some(false)) => {
                        cx.background_executor()
                            .timer(Duration::from_millis(8))
                            .await
                    }
                    _ => return,
                }
            }
        }));
        cx.notify();
    }

    /// Ask where to, then make the film.
    fn choose_movie_path(&mut self, cx: &mut Context<Self>) {
        self.commit_edit(cx);
        let o = self.options();
        if o.format.needs_ffmpeg() && movie::ffmpeg().is_none() {
            self.error =
                Some(format!("{} needs ffmpeg, which was not found", o.format.label()).into());
            cx.notify();
            return;
        }
        let stem = self
            .structure
            .as_ref()
            .map(|s| s.id.to_lowercase())
            .unwrap_or_else(|| "profold".into());
        let suggested = match o.format {
            Format::Png => format!("{stem}-frames"),
            f => format!("{stem}.{}", f.extension()),
        };
        let directory = dirs::home_dir()
            .map(|h| h.join("Desktop"))
            .unwrap_or_else(|| PathBuf::from("."));
        let chosen = cx.prompt_for_new_path(&directory, Some(&suggested));
        self.export = Some(cx.spawn(async move |this, cx| {
            let Ok(Ok(Some(path))) = chosen.await else {
                return;
            };
            let _ = this.update(cx, |this, cx| this.export_movie(path, o, cx));
        }));
    }

    /// Make the film at `path`, with progress in the dialog.
    pub(super) fn export_movie(
        &mut self,
        path: PathBuf,
        o: Options,
        cx: &mut Context<Self>,
    ) -> bool {
        if self.movie.progress.is_some() {
            self.error = Some("a film is already being made".into());
            cx.notify();
            return false;
        }
        let Some(job) = self.movie_job(&o) else {
            return false;
        };
        self.stop_preview();
        let progress = Arc::new(movie::Progress::default());
        progress.total.store(job.shots.len(), Ordering::Relaxed);
        self.movie.progress = Some(progress.clone());
        let result = Arc::new(Mutex::new(None));
        {
            let (progress, result) = (progress.clone(), result.clone());
            cx.background_executor()
                .spawn(async move {
                    let r = movie::export(&job, &path, &progress);
                    *result.lock().expect("not poisoned") = Some(r);
                })
                .detach();
        }
        self.movie.task = Some(cx.spawn(async move |this, cx| loop {
            cx.background_executor()
                .timer(Duration::from_millis(120))
                .await;
            let done = result.lock().expect("not poisoned").take();
            let going = this.update(cx, |this, cx| {
                let Some(r) = done else {
                    cx.notify();
                    return true;
                };
                this.movie.progress = None;
                match r {
                    Ok(path) => {
                        this.note(format!("Saved {}", path.display()));
                        cx.reveal_path(&path);
                    }
                    Err(e) if e.to_string() == "cancelled" => {
                        this.note("Film cancelled".to_string())
                    }
                    Err(e) => this.error = Some(format!("{e:#}").into()),
                }
                cx.notify();
                false
            });
            if !matches!(going, Ok(true)) {
                return;
            }
        }));
        cx.notify();
        true
    }

    fn cancel_movie(&mut self, cx: &mut Context<Self>) {
        if let Some(p) = &self.movie.progress {
            p.cancel.store(true, Ordering::Relaxed);
        }
        cx.notify();
    }

    /// Start keeping the camera and frame shown, thirty times a second,
    /// for a film to follow later.
    pub(super) fn start_recording(&mut self, cx: &mut Context<Self>) {
        if self.structure.is_none() || self.movie.recording.is_some() {
            return;
        }
        self.commit_edit(cx);
        self.movie.open = false;
        self.stop_preview();
        let task = cx.spawn(async move |this, cx| loop {
            cx.background_executor()
                .timer(Duration::from_secs_f32(1.0 / movie::RECORD_HZ))
                .await;
            let going = this.update(cx, |this, cx| {
                let frame = this.structure.as_ref().map(|s| s.frame).unwrap_or(0);
                let camera = this.camera;
                let Some(rec) = this.movie.recording.as_mut() else {
                    return false;
                };
                let t = rec.started.elapsed().as_secs_f32();
                rec.samples.push(Sample { t, camera, frame });
                // Ten minutes is a long film.
                if t > 600.0 {
                    this.stop_recording(cx);
                    return false;
                }
                cx.notify();
                true
            });
            if !matches!(going, Ok(true)) {
                return;
            }
        });
        let frame = self.structure.as_ref().map(|s| s.frame).unwrap_or(0);
        self.movie.recording = Some(Recording {
            started: Instant::now(),
            samples: vec![Sample {
                t: 0.0,
                camera: self.camera,
                frame,
            }],
            _task: task,
        });
        self.note("Recording: move, turn and play the trajectory as the film should".to_string());
        cx.notify();
    }

    /// Stop recording, and open the dialog to make a film of it.
    pub(super) fn stop_recording(&mut self, cx: &mut Context<Self>) -> Option<f32> {
        let rec = self.movie.recording.take()?;
        let mut samples = rec.samples;
        let t = rec.started.elapsed().as_secs_f32();
        let frame = self.structure.as_ref().map(|s| s.frame).unwrap_or(0);
        samples.push(Sample {
            t,
            camera: self.camera,
            frame,
        });
        self.movie.recorded = Some(Arc::new(samples));
        let mut o = self.options();
        o.motion = Motion::Recorded;
        self.set_options(o, cx);
        self.open_movie(cx);
        Some(t)
    }

    pub(super) fn recording(&self) -> bool {
        self.movie.recording.is_some()
    }

    /// The red light in the top bar while recording, with its stop. Up
    /// there it covers nothing in the viewport.
    pub(super) fn recording_badge(&self, cx: &mut Context<Self>) -> Option<gpui::AnyElement> {
        let rec = self.movie.recording.as_ref()?;
        let t = rec.started.elapsed().as_secs_f32();
        Some(
            div()
                .id("recording")
                .debug_selector(|| "recording".into())
                .flex()
                .flex_none()
                .items_center()
                .gap(px(7.))
                .mr(px(14.))
                .px(px(10.))
                .py(px(2.))
                .rounded(px(12.))
                .bg(theme::sunken())
                .border_1()
                .border_color(theme::border())
                .cursor_pointer()
                .on_mouse_down(MouseButton::Left, swallow)
                .on_click(cx.listener(|this, _, _, cx| {
                    this.stop_recording(cx);
                }))
                .child(div().size(px(8.)).rounded_full().bg(gpui::rgb(0xe5484d)))
                .child(
                    div()
                        .font_family(theme::MONO)
                        .text_size(px(10.5))
                        .text_color(theme::text())
                        .child(clock(t)),
                )
                .child(
                    div()
                        .text_size(px(10.5))
                        .text_color(theme::muted())
                        .child("Stop recording"),
                )
                .into_any_element(),
        )
    }

    /// Take what was typed into one of the dialog's values.
    pub(super) fn commit_movie(&mut self, knob: Knob, text: &str, cx: &mut Context<Self>) {
        let n = || {
            text.trim()
                .trim_end_matches(['°', 's', '%'])
                .trim()
                .parse::<f32>()
                .ok()
                .filter(|v| v.is_finite())
        };
        let mut o = self.options();
        match knob {
            Knob::Degrees => o.degrees = n().unwrap_or(o.degrees),
            Knob::First => o.first = n().map(|v| v as usize).unwrap_or(o.first),
            Knob::Last => {
                let count = self.frame_count();
                o.last = n()
                    .map(|v| if v as usize >= count { 0 } else { v as usize })
                    .unwrap_or(o.last);
            }
            Knob::Every => o.every = n().map(|v| v as usize).unwrap_or(o.every),
            Knob::Fps => o.fps = n().map(|v| v as u32).unwrap_or(o.fps),
            Knob::Seconds => {
                let t = text.trim();
                o.seconds = if t.is_empty() || t.eq_ignore_ascii_case("auto") {
                    0.0
                } else {
                    n().unwrap_or(o.seconds)
                };
            }
            Knob::Width => o.width = n().map(|v| v as u32).unwrap_or(o.width),
            Knob::Height => o.height = n().map(|v| v as u32).unwrap_or(o.height),
            Knob::Quality => o.quality = n().map(|v| v as u32).unwrap_or(o.quality),
            Knob::Dt => o.dt = n().unwrap_or(o.dt),
            Knob::Caption => o.caption = text.trim().to_string(),
            Knob::LegendTitle => o.legend_title = text.trim().to_string(),
        }
        self.set_options(o, cx);
    }

    // -- the dialog -------------------------------------------------------

    /// A row of choices, the one in use filled.
    fn choices<T: Copy + PartialEq + 'static>(
        &self,
        id: &'static str,
        items: &[(T, &'static str, bool)],
        current: T,
        set: fn(&mut Options, T),
        cx: &mut Context<Self>,
    ) -> gpui::Div {
        div()
            .flex()
            .flex_wrap()
            .gap(px(2.))
            .children(items.iter().enumerate().map(|(i, &(v, name, enabled))| {
                let p = pill_keyed((id, i), name, v == current);
                if enabled {
                    p.cursor_pointer()
                        .on_click(cx.listener(move |this, _, _, cx| this.change(cx, |o| set(o, v))))
                        .into_any_element()
                } else {
                    p.opacity(0.4).into_any_element()
                }
            }))
    }

    /// A choice that is on or off.
    fn switch(
        &self,
        id: &'static str,
        name: &'static str,
        on: bool,
        flip: fn(&mut Options),
        cx: &mut Context<Self>,
    ) -> impl IntoElement {
        pill_keyed(
            (id, 0),
            if on {
                format!("✓ {name}")
            } else {
                name.to_string()
            },
            on,
        )
        .cursor_pointer()
        .on_click(cx.listener(move |this, _, _, cx| this.change(cx, flip)))
    }

    /// The colour in use, which opens the colour picker for `pick` when
    /// clicked, as an atom's colour does.
    fn color_button(
        &self,
        id: &'static str,
        pick: Pick,
        name: &str,
        default: [f32; 3],
        chosen: Option<[f32; 3]>,
        cx: &mut Context<Self>,
    ) -> gpui::AnyElement {
        let open = self.picking.as_ref() == Some(&pick);
        let shown = chosen.unwrap_or(default);
        let toggle = pick.clone();
        div()
            .relative()
            .child(
                div()
                    .id(id)
                    .debug_selector(move || id.into())
                    .flex()
                    .items_center()
                    .gap(px(6.))
                    .px(px(6.))
                    .py(px(2.))
                    .rounded(px(5.))
                    .border_1()
                    .border_color(if open { theme::ink() } else { theme::border() })
                    .cursor_pointer()
                    .hover(|s| s.bg(theme::sunken()))
                    .child(
                        div()
                            .size(px(12.))
                            .rounded(px(3.))
                            .border_1()
                            .border_color(theme::border())
                            .bg(swatch(shown)),
                    )
                    .child(
                        div()
                            .font_family(theme::MONO)
                            .text_size(px(10.))
                            .text_color(theme::muted())
                            .child(settings::hex(shown)),
                    )
                    .on_click(cx.listener(move |this, _, _, cx| {
                        this.picking = if this.picking.as_ref() == Some(&toggle) {
                            None
                        } else {
                            Some(toggle.clone())
                        };
                        cx.notify();
                    })),
            )
            .when(open, |d| {
                d.child(
                    gpui::deferred(self.color_picker(pick, name.to_string(), default, chosen, cx))
                        .with_priority(3),
                )
            })
            .into_any_element()
    }

    fn knob(
        &self,
        id: &'static str,
        knob: Knob,
        text: String,
        width: f32,
        cx: &mut Context<Self>,
    ) -> gpui::AnyElement {
        self.field(id, Field::Movie(knob), text, width, cx)
    }

    /// The film dialog, over the viewport, while it is open.
    pub(super) fn movie_dialog(&self, cx: &mut Context<Self>) -> Option<gpui::AnyElement> {
        if !self.movie.open {
            return None;
        }
        let o = self.options();
        let frames = self.frame_count();
        let has_traj = frames > 1;
        let traj = has_traj && o.trajectory && o.motion != Motion::Recorded;
        let ffmpeg = movie::ffmpeg().is_some();
        let (count, seconds) = self.movie_length(&o);
        let (w, h) = o.output_size(self.cache.borrow().size);

        let heading = |text: &str| {
            div()
                .pt(px(6.))
                .text_size(px(9.))
                .text_color(theme::faint())
                .child(spaced(&text.to_uppercase()))
        };
        let row = |name: &str, content: gpui::AnyElement| {
            div()
                .flex()
                .items_center()
                .gap(px(8.))
                .min_h(px(22.))
                .child(
                    div()
                        .w(px(78.))
                        .flex_none()
                        .text_size(px(10.5))
                        .text_color(theme::muted())
                        .child(name.to_string()),
                )
                .child(
                    div()
                        .flex_1()
                        .flex()
                        .flex_wrap()
                        .items_center()
                        .gap(px(4.))
                        .child(content),
                )
        };
        let note = |text: String| {
            div()
                .text_size(px(9.5))
                .text_color(theme::faint())
                .child(text)
        };
        let inline = |items: Vec<gpui::AnyElement>| {
            div()
                .flex()
                .flex_wrap()
                .items_center()
                .gap(px(4.))
                .children(items)
                .into_any_element()
        };
        let small = |t: &str| {
            div()
                .text_size(px(10.))
                .text_color(theme::faint())
                .child(t.to_string())
                .into_any_element()
        };

        // -- Motion
        let mut motions = vec![
            (Motion::Still, "Still", true),
            (Motion::Turn, "Turn", true),
            (Motion::Rock, "Rock", true),
        ];
        if self.movie.recorded.is_some() {
            motions.push((Motion::Recorded, "Recorded path", true));
        }
        let camera_row = self
            .choices("movie-motion", &motions, o.motion, |o, v| o.motion = v, cx)
            .into_any_element();
        let turning = matches!(o.motion, Motion::Turn | Motion::Rock);
        let degrees = inline(vec![
            self.knob(
                "movie-degrees",
                Knob::Degrees,
                format!("{}°", o.degrees),
                52.,
                cx,
            ),
            small(if o.motion == Motion::Rock {
                "either side"
            } else {
                "in all"
            }),
        ]);
        let axes = self
            .choices(
                "movie-axis",
                &[
                    (Axis::Vertical, "Upright", true),
                    (Axis::Horizontal, "Across", true),
                    (Axis::View, "Line of sight", true),
                ],
                o.axis,
                |o, v| o.axis = v,
                cx,
            )
            .into_any_element();
        let record = pill_keyed(("movie-record", 0), "● Record a camera path", false)
            .cursor_pointer()
            .debug_selector(|| "movie-record".into())
            .on_click(cx.listener(|this, _, _, cx| this.start_recording(cx)))
            .into_any_element();
        let recorded_note = self
            .movie
            .recorded
            .as_ref()
            .filter(|_| o.motion == Motion::Recorded)
            .map(|r| {
                let t = r.last().map(|s| s.t).unwrap_or(0.0);
                note(format!(
                    "{} recorded, played back as it was shown, trajectory and all",
                    clock(t)
                ))
            });

        // -- Frames
        let range = inline(vec![
            self.knob(
                "movie-first",
                Knob::First,
                o.first.min(frames).to_string(),
                44.,
                cx,
            ),
            small("to"),
            self.knob(
                "movie-last",
                Knob::Last,
                if o.last == 0 || o.last > frames {
                    frames.to_string()
                } else {
                    o.last.to_string()
                },
                44.,
                cx,
            ),
            small("every"),
            self.knob("movie-every", Knob::Every, o.every.to_string(), 32., cx),
        ]);
        let smooth = self
            .choices(
                "movie-smooth",
                &[
                    (0, "Off", true),
                    (1, "2×", true),
                    (3, "4×", true),
                    (7, "8×", true),
                ],
                o.smooth,
                |o, v| o.smooth = v,
                cx,
            )
            .into_any_element();
        let pingpong = self
            .switch(
                "movie-pingpong",
                "Ping-pong",
                o.pingpong,
                |o| o.pingpong = !o.pingpong,
                cx,
            )
            .into_any_element();
        let timing = inline(vec![
            self.knob("movie-fps", Knob::Fps, o.fps.to_string(), 36., cx),
            small("fps ·"),
            self.knob(
                "movie-seconds",
                Knob::Seconds,
                if o.seconds > 0.0 {
                    format!("{}s", o.seconds)
                } else {
                    "auto".into()
                },
                48.,
                cx,
            ),
            small("long"),
        ]);

        // -- Image
        let sizes: Vec<_> = Size::ALL.iter().map(|s| (*s, s.label(), true)).collect();
        let size = self
            .choices("movie-size", &sizes, o.size, |o, v| o.size = v, cx)
            .into_any_element();
        let custom = (o.size == Size::Custom).then(|| {
            inline(vec![
                self.knob("movie-width", Knob::Width, o.width.to_string(), 48., cx),
                small("×"),
                self.knob("movie-height", Knob::Height, o.height.to_string(), 48., cx),
            ])
        });
        let aa = self
            .choices(
                "movie-aa",
                &[(1, "Off", true), (2, "2×", true), (4, "4×", true)],
                o.supersample,
                |o, v| o.supersample = v,
                cx,
            )
            .into_any_element();
        // Traced, how long the film will be, as this machine last traced.
        let traced_note: String = match (o.traced, crate::trace::estimate(w as usize, h as usize)) {
            (true, Some(t)) => format!(
                "≈ {} a frame, {} for all {count}",
                crate::trace::roughly(t),
                crate::trace::roughly(t * count as u32)
            ),
            (true, None) => format!("seconds a frame, so minutes for all {count}"),
            (false, _) => "shadows and ambient occlusion; seconds a frame".into(),
        };
        let traced = inline(vec![
            self.switch(
                "movie-traced",
                "Ray traced",
                o.traced,
                |o| o.traced = !o.traced,
                cx,
            )
            .into_any_element(),
            small(&traced_note),
        ]);
        let quality = inline(vec![
            self.knob(
                "movie-quality",
                Knob::Quality,
                format!("{}", o.quality),
                36.,
                cx,
            ),
            small("of 100"),
        ]);
        let backgrounds = self
            .choices(
                "movie-bg",
                &[
                    (0u8, "As shown", true),
                    (1u8, "Colour", true),
                    (2u8, "Transparent", o.format.keeps_alpha()),
                ],
                match o.background {
                    Background::Theme => 0,
                    Background::Color(_) => 1,
                    Background::Transparent => 2,
                },
                |o, v| {
                    o.background = match v {
                        0 => Background::Theme,
                        1 => Background::Color(match o.background {
                            Background::Color(c) => c,
                            _ => theme::view(),
                        }),
                        _ => Background::Transparent,
                    }
                },
                cx,
            )
            .into_any_element();
        let hex = match o.background {
            Background::Color(c) => {
                let page = theme::view();
                Some(self.color_button(
                    "movie-bg-color",
                    Pick::FilmBackground,
                    "the background",
                    page,
                    (c != page).then_some(c),
                    cx,
                ))
            }
            _ => None,
        };

        // -- Format
        let formats: Vec<_> = Format::ALL
            .iter()
            .map(|f| (*f, f.label(), !f.needs_ffmpeg() || ffmpeg))
            .collect();
        let format = self
            .choices(
                "movie-format",
                &formats,
                o.format,
                |o, v| {
                    o.format = v;
                    if !v.keeps_alpha() && o.background == Background::Transparent {
                        o.background = Background::Theme;
                    }
                },
                cx,
            )
            .into_any_element();
        let format_note = match o.format {
            Format::Mp4 => "H.264, plays everywhere: slides, the web, phones",
            Format::Gif => "For chat and READMEs; large, 256 colours, hard-edged transparency",
            Format::Png => "A numbered PNG a frame in a folder, lossless, for editing",
            Format::WebM => "VP9 with soft transparency, for the web",
            Format::Mov => "ProRes 4444 with transparency, for Final Cut, Premiere, Keynote",
        };
        let ffmpeg_note =
            (!ffmpeg).then(|| note("MP4, WebM and ProRes need ffmpeg installed".to_string()));

        // -- Overlays
        let overlays = inline(vec![
            self.switch(
                "movie-counter",
                "Frame counter",
                o.counter,
                |o| o.counter = !o.counter,
                cx,
            )
            .into_any_element(),
            self.switch("movie-time", "Time", o.time, |o| o.time = !o.time, cx)
                .into_any_element(),
            self.switch(
                "movie-scale",
                "Scale bar",
                o.scale_bar,
                |o| o.scale_bar = !o.scale_bar,
                cx,
            )
            .into_any_element(),
            self.switch(
                "movie-legend",
                "Legend",
                o.legend,
                |o| o.legend = !o.legend,
                cx,
            )
            .into_any_element(),
        ]);
        let sizes_of = |id: &'static str,
                        current: TextSize,
                        set: fn(&mut Options, TextSize),
                        cx: &mut Context<Self>| {
            let items: Vec<_> = TextSize::ALL
                .iter()
                .map(|t| (*t, t.label(), true))
                .collect();
            self.choices(id, &items, current, set, cx)
                .into_any_element()
        };
        let film_legend = self.film_legend();
        let legend = o.legend.then(|| match &film_legend {
            Some(_) => {
                let corners: Vec<_> = Corner::ALL.iter().map(|c| (*c, c.label(), true)).collect();
                self.choices("movie-corner", &corners, o.corner, |o, v| o.corner = v, cx)
                    .into_any_element()
            }
            None => note("This colouring has no key to show".to_string()).into_any_element(),
        });
        // How the legend reads, once there is one to show.
        let legend_look = film_legend.as_ref().filter(|_| o.legend).map(|l| {
            let title = self.knob(
                "movie-legend-title",
                Knob::LegendTitle,
                if o.legend_title.is_empty() {
                    l.title().to_string()
                } else {
                    o.legend_title.clone()
                },
                160.,
                cx,
            );
            let backdrop = self
                .switch(
                    "movie-legend-box",
                    "Backdrop",
                    o.legend_box,
                    |o| o.legend_box = !o.legend_box,
                    cx,
                )
                .into_any_element();
            (
                inline(vec![title, backdrop]),
                sizes_of(
                    "movie-legend-size",
                    o.legend_size,
                    |o, v| o.legend_size = v,
                    cx,
                ),
            )
        });
        let steps = self
            .structure
            .as_ref()
            .is_some_and(|s| s.frames.iter().any(|f| f.step.is_some()));
        let time = (o.time).then(|| {
            if !steps {
                return note("This trajectory does not say which step each frame is".to_string())
                    .into_any_element();
            }
            inline(vec![
                small("timestep"),
                self.knob("movie-dt", Knob::Dt, format!("{}", o.dt), 52., cx),
                self.choices(
                    "movie-unit",
                    &[
                        (TimeUnit::Fs, "fs", true),
                        (TimeUnit::Ps, "ps", true),
                        (TimeUnit::Tau, "τ", true),
                    ],
                    o.unit,
                    |o, v| o.unit = v,
                    cx,
                )
                .into_any_element(),
            ])
        });
        let caption = self.knob(
            "movie-caption",
            Knob::Caption,
            if o.caption.is_empty() {
                "none".into()
            } else {
                o.caption.clone()
            },
            200.,
            cx,
        );
        let caption_look = (!o.caption.trim().is_empty()).then(|| {
            let p = o.caption_place;
            (
                inline(vec![
                    self.choices(
                        "movie-caption-edge",
                        &[(false, "Top", true), (true, "Bottom", true)],
                        p.bottom(),
                        |o, v| o.caption_place = Place::at(v, o.caption_place.across()),
                        cx,
                    )
                    .into_any_element(),
                    small("·"),
                    self.choices(
                        "movie-caption-across",
                        &[(0, "Left", true), (1, "Centre", true), (2, "Right", true)],
                        p.across(),
                        |o, v| o.caption_place = Place::at(o.caption_place.bottom(), v),
                        cx,
                    )
                    .into_any_element(),
                ]),
                sizes_of(
                    "movie-caption-size",
                    o.caption_size,
                    |o, v| o.caption_size = v,
                    cx,
                ),
            )
        });
        let writing = o.legend && film_legend.is_some()
            || !o.caption.trim().is_empty()
            || o.counter
            || o.time
            || o.scale_bar;
        let ink = writing.then(|| {
            let mut items = vec![self
                .choices(
                    "movie-ink",
                    &[(false, "Automatic", true), (true, "Colour", true)],
                    o.ink.is_some(),
                    |o, v| {
                        o.ink = v.then(|| o.ink.unwrap_or([0.85, 0.2, 0.2]));
                    },
                    cx,
                )
                .into_any_element()];
            if o.ink.is_some() {
                items.push(self.color_button(
                    "movie-ink-color",
                    Pick::FilmInk,
                    "the writing",
                    movie::auto_ink(o.background, theme::view()),
                    o.ink,
                    cx,
                ));
            }
            inline(items)
        });

        let left = div()
            .w(px(410.))
            .flex_none()
            .flex()
            .flex_col()
            .gap(px(4.))
            .child(heading("Motion"))
            .when(has_traj && o.motion != Motion::Recorded, |d| {
                d.child(row(
                    "Trajectory",
                    self.switch(
                        "movie-traj",
                        "Play the frames",
                        o.trajectory,
                        |o| o.trajectory = !o.trajectory,
                        cx,
                    )
                    .into_any_element(),
                ))
            })
            .child(row("Camera", camera_row))
            .when(turning, |d| {
                d.child(row("Angle", degrees)).child(row("About", axes))
            })
            .when(turning && o.motion == Motion::Turn, |d| {
                d.child(row(
                    "",
                    self.switch(
                        "movie-seamless",
                        "Loop seamlessly",
                        o.seamless,
                        |o| o.seamless = !o.seamless,
                        cx,
                    )
                    .into_any_element(),
                ))
            })
            .children(recorded_note)
            .child(row("", record))
            .when(traj || o.motion == Motion::Recorded && has_traj, |d| {
                d.child(heading("Frames"))
            })
            .children(traj.then(|| row("Range", range)))
            .children(
                (traj || o.motion == Motion::Recorded && has_traj)
                    .then(|| row("Smoothing", smooth)),
            )
            .children(traj.then(|| row("", pingpong)))
            .child(row("Timing", timing))
            .child(heading("Image"))
            .child(row("Size", size))
            .children(custom.map(|c| row("", c)))
            .when(!o.traced, |d| d.child(row("Antialias", aa)))
            .child(row("", traced))
            .when(o.format.lossy(), |d| d.child(row("Quality", quality)))
            .child(row("Background", backgrounds))
            .children(hex.map(|c| row("", c)))
            .child(heading("Format"))
            .child(row("", format))
            .child(note(format_note.to_string()))
            .children(ffmpeg_note)
            .child(heading("Overlays"))
            .child(row("Show", overlays))
            .children(time.map(|t| row("", t)))
            .children(legend.map(|l| row("Legend", l)))
            .children(legend_look.map(|(title, size)| {
                div()
                    .flex()
                    .flex_col()
                    .gap(px(4.))
                    .child(row("Legend title", title))
                    .child(row("Legend size", size))
            }))
            .child(row("Caption", caption))
            .children(caption_look.map(|(place, size)| {
                div()
                    .flex()
                    .flex_col()
                    .gap(px(4.))
                    .child(row("Placed", place))
                    .child(row("Caption size", size))
            }))
            .children(ink.map(|i| row("Text colour", i)));

        // -- Preview and the buttons
        let ratio = h as f32 / w.max(1) as f32;
        let wide = self.preview_width();
        let ph = (wide * ratio).min(wide * 1.2);
        let pw = ph / ratio;
        let shown = self.movie.shown.clone();
        let screen = div()
            .id("movie-preview")
            .debug_selector(|| "movie-preview".into())
            .w(px(pw))
            .h(px(ph))
            .rounded(px(6.))
            .overflow_hidden()
            .border_1()
            .border_color(theme::border())
            .bg(theme::sunken())
            .flex()
            .items_center()
            .justify_center()
            .when_some(shown.clone(), |d, (_, image)| {
                // What is showing is out of date while the next is drawn.
                let stale = self.movie.previewing;
                d.child(
                    gpui::img(gpui::ImageSource::Render(image))
                        .size_full()
                        .when(stale, |i| i.opacity(0.5)),
                )
            })
            .when(shown.is_none(), |d| {
                d.child(div().text_size(px(10.5)).text_color(theme::faint()).child(
                    if self.movie.previewing {
                        "Drawing the preview…"
                    } else {
                        "Nothing to preview"
                    },
                ))
            });
        let at = shown.as_ref().map(|(shot, _)| {
            let updating = if self.movie.previewing {
                " · updating…"
            } else {
                ""
            };
            note(format!(
                "{} / {} · {}{updating}",
                shot + 1,
                count,
                clock(*shot as f32 / o.fps.max(1) as f32)
            ))
        });
        let summary = format!(
            "{} × {} · {} frames · {} at {} fps",
            w,
            h,
            count,
            clock(seconds),
            o.fps
        );
        let busy = self.movie.progress.clone();
        let primary = |id: &'static str, text: &'static str| {
            div()
                .id(id)
                .debug_selector(move || id.into())
                .px(px(14.))
                .py(px(5.))
                .rounded(px(6.))
                .bg(theme::ink())
                .text_color(theme::on_ink())
                .text_size(px(11.5))
                .cursor_pointer()
                .child(text)
        };
        let plain = |id: &'static str, text: &'static str| {
            div()
                .id(id)
                .debug_selector(move || id.into())
                .px(px(12.))
                .py(px(5.))
                .rounded(px(6.))
                .border_1()
                .border_color(theme::border())
                .text_color(theme::muted())
                .text_size(px(11.5))
                .cursor_pointer()
                .hover(|s| s.bg(theme::sunken()).text_color(theme::text()))
                .child(text)
        };
        let actions: gpui::AnyElement = match &busy {
            Some(p) => {
                let (done, total) = (
                    p.done.load(Ordering::Relaxed),
                    p.total.load(Ordering::Relaxed).max(1),
                );
                let f = done as f32 / total as f32;
                div()
                    .flex()
                    .flex_col()
                    .gap(px(6.))
                    .child(
                        div()
                            .w_full()
                            .h(px(5.))
                            .rounded(px(3.))
                            .bg(theme::sunken())
                            .child(
                                div()
                                    .h_full()
                                    .w(relative(f))
                                    .rounded(px(3.))
                                    .bg(theme::ink()),
                            ),
                    )
                    .child(
                        div()
                            .flex()
                            .justify_between()
                            .items_center()
                            .child(note(format!("Frame {done} of {total}")))
                            .child(
                                plain("movie-cancel", "Cancel")
                                    .on_click(cx.listener(|this, _, _, cx| this.cancel_movie(cx))),
                            ),
                    )
                    .into_any_element()
            }
            None => div()
                .flex()
                .gap(px(6.))
                .justify_end()
                .child(
                    plain("movie-close", "Close")
                        .on_click(cx.listener(|this, _, _, cx| this.close_movie(cx))),
                )
                .child(
                    primary("movie-export", "Export…")
                        .on_click(cx.listener(|this, _, _, cx| this.choose_movie_path(cx))),
                )
                .into_any_element(),
        };
        let right = div()
            .w(px(wide))
            .flex_none()
            .flex()
            .flex_col()
            .gap(px(8.))
            .pt(px(6.))
            .child(screen)
            .children(at)
            .child(
                div()
                    .text_size(px(10.5))
                    .text_color(theme::text())
                    .child(SharedString::from(summary)),
            )
            .when(o.format == Format::Gif && o.fps > 50, |d| {
                d.child(note("GIFs play at most 50 fps".to_string()))
            })
            .child(div().flex_1())
            .child(actions);

        Some(
            div()
                .id("movie")
                .debug_selector(|| "movie".into())
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
                .on_scroll_wheel(|_, _, cx| cx.stop_propagation())
                .child(
                    div()
                        .id("movie-card")
                        .max_h(relative(0.94))
                        .overflow_y_scroll()
                        .p(px(18.))
                        .rounded(px(10.))
                        .bg(theme::surface())
                        .border_1()
                        .border_color(theme::border())
                        .shadow_lg()
                        .flex()
                        .flex_col()
                        .gap(px(4.))
                        .child(
                            div()
                                .text_size(px(13.))
                                .text_color(theme::ink())
                                .child("Export a film"),
                        )
                        .child(div().flex().gap(px(22.)).child(left).child(right)),
                )
                .into_any_element(),
        )
    }

    /// `movie <file> [key=value …]`: make a film without the dialog.
    pub(super) fn movie_command(
        &mut self,
        rest: &str,
        cx: &mut Context<Self>,
    ) -> anyhow::Result<String> {
        let words = movie::words(rest);
        // Alone, it opens the dialog.
        let Some((file, args)) = words.split_first() else {
            self.open_movie(cx);
            return Ok("the film dialog is open".into());
        };
        let path = PathBuf::from(file.trim_matches('"'));
        let o = movie::parse_args(Options::default(), &path, args)?;
        if o.motion == Motion::Recorded && self.movie.recorded.is_none() {
            anyhow::bail!("nothing has been recorded; use record, then record stop");
        }
        let (n, secs) = self.movie_length(&o);
        if !self.export_movie(path.clone(), o, cx) {
            anyhow::bail!("could not start the film");
        }
        Ok(format!(
            "making {} ({n} frames, {})",
            path.display(),
            clock(secs)
        ))
    }
}
