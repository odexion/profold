//! The film dialog: what moves, which frames, how large and on what, in
//! which format and with what written over it; a preview, the export with
//! its progress, and recording a camera path by hand to play back.

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use gpui::{div, prelude::*, px, relative, Context, MouseButton, SharedString, Task};

use super::{pill_keyed, spaced, swallow, Field, Viewer};
use crate::movie::{self, Axis, Background, Format, Motion, Options, Sample, Size, TimeUnit};
use crate::protein::Structure;
use crate::render::{self, Scene};
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
    Hex,
    Dt,
    Caption,
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
    /// Small frames of the film and how many of its frames each stands for,
    /// and the one showing.
    preview: Vec<Arc<gpui::RenderImage>>,
    preview_stride: usize,
    preview_at: usize,
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
    let Some(p) = settings::file("movie.json") else { return };
    if let Some(dir) = p.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let _ = std::fs::write(p, o.to_json());
}

/// Seconds as minutes and seconds, to a tenth.
fn clock(s: f32) -> String {
    format!("{}:{:04.1}", (s / 60.0) as u32, s % 60.0)
}

const PREVIEW_WIDTH: f32 = 296.0;

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
        self.movie.open = true;
        self.refresh_preview(Duration::ZERO, cx);
        cx.notify();
    }

    pub(super) fn close_movie(&mut self, cx: &mut Context<Self>) {
        self.commit_edit(cx);
        self.movie.open = false;
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
        let (repr, coloring, style, palette, standard, cell) =
            (self.repr, self.coloring, self.style, self.palette, self.standard_colors, self.show_cell);
        let mut appearance = self.appearance.clone();
        appearance.hidden = None;
        let layers: Vec<_> = self.layers.iter().map(|l| (l.selection.clone(), l.style, l.color)).collect();
        let main = self.shown(0).is_some();
        let other = self.shown(1).map(|o| Structure::clone(o));
        let build = move |s: &Structure| -> Scene {
            let mut scene = Scene::default();
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
                scene.merge(render::build_scene(s, repr, coloring, 1.0, style, palette, standard, cell, &look));
                for (atoms, st, c) in &found {
                    scene.merge(render::build_layer(s, atoms, *st, *c, style, palette, standard, &look));
                }
            }
            if let Some(other) = &other {
                scene.merge(render::build_scene(other, repr, coloring, 1.0, style, palette, standard, cell, &appearance));
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
        })
    }

    /// How many frames the film would have, and how long it would run.
    pub(super) fn movie_length(&self, o: &Options) -> (usize, f32) {
        let (now, count) = (self.structure.as_ref().map(|s| s.frame).unwrap_or(0), self.frame_count());
        let recorded = self.movie.recorded.as_deref().map(|v| v.as_slice());
        let n = movie::plan(o, &self.camera, now, count, recorded).len();
        (n, n as f32 / o.fps.max(1) as f32)
    }

    /// Draw a few small frames of the film and play them in the dialog.
    fn preview_movie(&mut self, cx: &mut Context<Self>) {
        let Some(job) = self.movie_job(&self.options()) else { return };
        self.stop_preview();
        let cancel = Arc::new(AtomicBool::new(false));
        self.movie.preview_cancel = Some(cancel.clone());
        self.movie.previewing = true;
        let fps = job.options.fps.max(1);
        let made = cx.background_executor().spawn(async move {
            let (pictures, stride) = movie::preview(&job, (PREVIEW_WIDTH * 2.0) as usize, 150, &cancel);
            let images: Vec<Arc<gpui::RenderImage>> = pictures
                .into_iter()
                .map(|p| {
                    let mut bgra = p.rgba;
                    for px in bgra.chunks_exact_mut(4) {
                        px.swap(0, 2);
                    }
                    let buffer = image::RgbaImage::from_raw(p.width as u32, p.height as u32, bgra).expect("sized to fit");
                    Arc::new(gpui::RenderImage::new(smallvec::smallvec![image::Frame::new(buffer)]))
                })
                .collect();
            (images, stride)
        });
        self.movie.preview_task = Some(cx.spawn(async move |this, cx| {
            let (images, stride) = made.await;
            let delay = Duration::from_secs_f32(stride as f32 / fps as f32);
            let ok = this.update(cx, |this, cx| {
                this.movie.previewing = false;
                this.movie.preview = images;
                this.movie.preview_stride = stride;
                this.movie.preview_at = 0;
                cx.notify();
                this.movie.preview.len() > 1
            });
            if !matches!(ok, Ok(true)) {
                return;
            }
            loop {
                cx.background_executor().timer(delay).await;
                let going = this.update(cx, |this, cx| {
                    let n = this.movie.preview.len();
                    if n == 0 || !this.movie.open {
                        return false;
                    }
                    this.movie.preview_at = (this.movie.preview_at + 1) % n;
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

    /// Ask where to, then make the film.
    fn choose_movie_path(&mut self, cx: &mut Context<Self>) {
        self.commit_edit(cx);
        let o = self.options();
        if o.format.needs_ffmpeg() && movie::ffmpeg().is_none() {
            self.error = Some(format!("{} needs ffmpeg, which was not found", o.format.label()).into());
            cx.notify();
            return;
        }
        let stem = self.structure.as_ref().map(|s| s.id.to_lowercase()).unwrap_or_else(|| "profold".into());
        let suggested = match o.format {
            Format::Png => format!("{stem}-frames"),
            f => format!("{stem}.{}", f.extension()),
        };
        let directory = dirs::home_dir().map(|h| h.join("Desktop")).unwrap_or_else(|| PathBuf::from("."));
        let chosen = cx.prompt_for_new_path(&directory, Some(&suggested));
        self.export = Some(cx.spawn(async move |this, cx| {
            let Ok(Ok(Some(path))) = chosen.await else { return };
            let _ = this.update(cx, |this, cx| this.export_movie(path, o, cx));
        }));
    }

    /// Make the film at `path`, with progress in the dialog.
    pub(super) fn export_movie(&mut self, path: PathBuf, o: Options, cx: &mut Context<Self>) -> bool {
        if self.movie.progress.is_some() {
            self.error = Some("a film is already being made".into());
            cx.notify();
            return false;
        }
        let Some(job) = self.movie_job(&o) else { return false };
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
            cx.background_executor().timer(Duration::from_millis(120)).await;
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
                    Err(e) if e.to_string() == "cancelled" => this.note("Film cancelled".to_string()),
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
            cx.background_executor().timer(Duration::from_secs_f32(1.0 / movie::RECORD_HZ)).await;
            let going = this.update(cx, |this, cx| {
                let frame = this.structure.as_ref().map(|s| s.frame).unwrap_or(0);
                let camera = this.camera;
                let Some(rec) = this.movie.recording.as_mut() else { return false };
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
            samples: vec![Sample { t: 0.0, camera: self.camera, frame }],
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
        samples.push(Sample { t, camera: self.camera, frame });
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

    /// The red light over the viewport while recording, with its stop.
    pub(super) fn recording_badge(&self, cx: &mut Context<Self>) -> Option<gpui::AnyElement> {
        let rec = self.movie.recording.as_ref()?;
        let t = rec.started.elapsed().as_secs_f32();
        Some(
            div()
                .absolute()
                .top(px(14.))
                .left_0()
                .right_0()
                .flex()
                .justify_center()
                .child(
                    div()
                        .id("recording")
                        .debug_selector(|| "recording".into())
                        .flex()
                        .items_center()
                        .gap(px(8.))
                        .px(px(12.))
                        .py(px(5.))
                        .rounded(px(14.))
                        .bg(theme::surface())
                        .border_1()
                        .border_color(theme::border())
                        .shadow_md()
                        .cursor_pointer()
                        .on_mouse_down(MouseButton::Left, swallow)
                        .on_click(cx.listener(|this, _, _, cx| {
                            this.stop_recording(cx);
                        }))
                        .child(div().size(px(9.)).rounded_full().bg(gpui::rgb(0xe5484d)))
                        .child(div().font_family(theme::MONO).text_size(px(11.)).text_color(theme::text()).child(clock(t)))
                        .child(div().text_size(px(11.)).text_color(theme::muted()).child("Stop recording")),
                )
                .into_any_element(),
        )
    }

    /// Take what was typed into one of the dialog's values.
    pub(super) fn commit_movie(&mut self, knob: Knob, text: &str, cx: &mut Context<Self>) {
        let n = || text.trim().trim_end_matches(['°', 's', '%']).trim().parse::<f32>().ok().filter(|v| v.is_finite());
        let mut o = self.options();
        match knob {
            Knob::Degrees => o.degrees = n().unwrap_or(o.degrees),
            Knob::First => o.first = n().map(|v| v as usize).unwrap_or(o.first),
            Knob::Last => {
                let count = self.frame_count();
                o.last = n().map(|v| if v as usize >= count { 0 } else { v as usize }).unwrap_or(o.last);
            }
            Knob::Every => o.every = n().map(|v| v as usize).unwrap_or(o.every),
            Knob::Fps => o.fps = n().map(|v| v as u32).unwrap_or(o.fps),
            Knob::Seconds => {
                let t = text.trim();
                o.seconds = if t.is_empty() || t.eq_ignore_ascii_case("auto") { 0.0 } else { n().unwrap_or(o.seconds) };
            }
            Knob::Width => o.width = n().map(|v| v as u32).unwrap_or(o.width),
            Knob::Height => o.height = n().map(|v| v as u32).unwrap_or(o.height),
            Knob::Quality => o.quality = n().map(|v| v as u32).unwrap_or(o.quality),
            Knob::Hex => {
                if let Some(c) = settings::parse_hex(text) {
                    o.background = Background::Color(c);
                }
            }
            Knob::Dt => o.dt = n().unwrap_or(o.dt),
            Knob::Caption => o.caption = text.trim().to_string(),
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
        div().flex().flex_wrap().gap(px(2.)).children(items.iter().enumerate().map(|(i, &(v, name, enabled))| {
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
    fn switch(&self, id: &'static str, name: &'static str, on: bool, flip: fn(&mut Options), cx: &mut Context<Self>) -> impl IntoElement {
        pill_keyed((id, 0), if on { format!("✓ {name}") } else { name.to_string() }, on)
            .cursor_pointer()
            .on_click(cx.listener(move |this, _, _, cx| this.change(cx, flip)))
    }

    fn knob(&self, id: &'static str, knob: Knob, text: String, width: f32, cx: &mut Context<Self>) -> gpui::AnyElement {
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

        let heading = |text: &str| div().pt(px(6.)).text_size(px(9.)).text_color(theme::faint()).child(spaced(&text.to_uppercase()));
        let row = |name: &str, content: gpui::AnyElement| {
            div()
                .flex()
                .items_center()
                .gap(px(8.))
                .min_h(px(22.))
                .child(div().w(px(78.)).flex_none().text_size(px(10.5)).text_color(theme::muted()).child(name.to_string()))
                .child(div().flex_1().flex().flex_wrap().items_center().gap(px(4.)).child(content))
        };
        let note = |text: String| div().text_size(px(9.5)).text_color(theme::faint()).child(text);
        let inline = |items: Vec<gpui::AnyElement>| div().flex().flex_wrap().items_center().gap(px(4.)).children(items).into_any_element();
        let small = |t: &str| div().text_size(px(10.)).text_color(theme::faint()).child(t.to_string()).into_any_element();

        // -- Motion
        let mut motions = vec![(Motion::Still, "Still", true), (Motion::Turn, "Turn", true), (Motion::Rock, "Rock", true)];
        if self.movie.recorded.is_some() {
            motions.push((Motion::Recorded, "Recorded path", true));
        }
        let camera_row = self.choices("movie-motion", &motions, o.motion, |o, v| o.motion = v, cx).into_any_element();
        let turning = matches!(o.motion, Motion::Turn | Motion::Rock);
        let degrees = inline(vec![
            self.knob("movie-degrees", Knob::Degrees, format!("{}°", o.degrees), 52., cx),
            small(if o.motion == Motion::Rock { "either side" } else { "in all" }),
        ]);
        let axes = self
            .choices(
                "movie-axis",
                &[(Axis::Vertical, "Upright", true), (Axis::Horizontal, "Across", true), (Axis::View, "Line of sight", true)],
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
        let recorded_note = self.movie.recorded.as_ref().filter(|_| o.motion == Motion::Recorded).map(|r| {
            let t = r.last().map(|s| s.t).unwrap_or(0.0);
            note(format!("{} recorded, played back as it was shown, trajectory and all", clock(t)))
        });

        // -- Frames
        let range = inline(vec![
            self.knob("movie-first", Knob::First, o.first.min(frames).to_string(), 44., cx),
            small("to"),
            self.knob("movie-last", Knob::Last, if o.last == 0 || o.last > frames { frames.to_string() } else { o.last.to_string() }, 44., cx),
            small("every"),
            self.knob("movie-every", Knob::Every, o.every.to_string(), 32., cx),
        ]);
        let smooth = self
            .choices(
                "movie-smooth",
                &[(0, "Off", true), (1, "2×", true), (3, "4×", true), (7, "8×", true)],
                o.smooth,
                |o, v| o.smooth = v,
                cx,
            )
            .into_any_element();
        let pingpong = self.switch("movie-pingpong", "Ping-pong", o.pingpong, |o| o.pingpong = !o.pingpong, cx).into_any_element();
        let timing = inline(vec![
            self.knob("movie-fps", Knob::Fps, o.fps.to_string(), 36., cx),
            small("fps ·"),
            self.knob("movie-seconds", Knob::Seconds, if o.seconds > 0.0 { format!("{}s", o.seconds) } else { "auto".into() }, 48., cx),
            small("long"),
        ]);

        // -- Image
        let sizes: Vec<_> = Size::ALL.iter().map(|s| (*s, s.label(), true)).collect();
        let size = self.choices("movie-size", &sizes, o.size, |o, v| o.size = v, cx).into_any_element();
        let custom = (o.size == Size::Custom).then(|| {
            inline(vec![
                self.knob("movie-width", Knob::Width, o.width.to_string(), 48., cx),
                small("×"),
                self.knob("movie-height", Knob::Height, o.height.to_string(), 48., cx),
            ])
        });
        let aa = self
            .choices("movie-aa", &[(1, "Off", true), (2, "2×", true), (4, "4×", true)], o.supersample, |o, v| o.supersample = v, cx)
            .into_any_element();
        let quality = inline(vec![self.knob("movie-quality", Knob::Quality, format!("{}", o.quality), 36., cx), small("of 100")]);
        let color = match o.background {
            Background::Color(c) => c,
            _ => theme::view(),
        };
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
        let hex = matches!(o.background, Background::Color(_)).then(|| {
            let b = color.map(|v| (v.clamp(0.0, 1.0) * 255.0).round() as u8);
            inline(vec![
                div()
                    .size(px(14.))
                    .rounded(px(3.))
                    .border_1()
                    .border_color(theme::border())
                    .bg(gpui::Rgba { r: color[0], g: color[1], b: color[2], a: 1.0 })
                    .into_any_element(),
                self.knob("movie-hex", Knob::Hex, format!("#{:02x}{:02x}{:02x}", b[0], b[1], b[2]), 64., cx),
            ])
        });

        // -- Format
        let formats: Vec<_> = Format::ALL.iter().map(|f| (*f, f.label(), !f.needs_ffmpeg() || ffmpeg)).collect();
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
        let ffmpeg_note = (!ffmpeg).then(|| note("WebM and ProRes need ffmpeg installed".to_string()));

        // -- Overlays
        let overlays = inline(vec![
            self.switch("movie-counter", "Frame counter", o.counter, |o| o.counter = !o.counter, cx).into_any_element(),
            self.switch("movie-time", "Time", o.time, |o| o.time = !o.time, cx).into_any_element(),
            self.switch("movie-scale", "Scale bar", o.scale_bar, |o| o.scale_bar = !o.scale_bar, cx).into_any_element(),
        ]);
        let steps = self.structure.as_ref().is_some_and(|s| s.frames.iter().any(|f| f.step.is_some()));
        let time = (o.time).then(|| {
            if !steps {
                return note("This trajectory does not say which step each frame is".to_string()).into_any_element();
            }
            inline(vec![
                small("timestep"),
                self.knob("movie-dt", Knob::Dt, format!("{}", o.dt), 52., cx),
                self.choices(
                    "movie-unit",
                    &[(TimeUnit::Fs, "fs", true), (TimeUnit::Ps, "ps", true), (TimeUnit::Tau, "τ", true)],
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
            if o.caption.is_empty() { "none".into() } else { o.caption.clone() },
            200.,
            cx,
        );

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
                    self.switch("movie-traj", "Play the frames", o.trajectory, |o| o.trajectory = !o.trajectory, cx).into_any_element(),
                ))
            })
            .child(row("Camera", camera_row))
            .when(turning, |d| d.child(row("Angle", degrees)).child(row("About", axes)))
            .when(turning && o.motion == Motion::Turn, |d| {
                d.child(row(
                    "",
                    self.switch("movie-seamless", "Loop seamlessly", o.seamless, |o| o.seamless = !o.seamless, cx).into_any_element(),
                ))
            })
            .children(recorded_note)
            .child(row("", record))
            .when(traj || o.motion == Motion::Recorded && has_traj, |d| d.child(heading("Frames")))
            .children(traj.then(|| row("Range", range)))
            .children((traj || o.motion == Motion::Recorded && has_traj).then(|| row("Smoothing", smooth)))
            .children(traj.then(|| row("", pingpong)))
            .child(row("Timing", timing))
            .child(heading("Image"))
            .child(row("Size", size))
            .children(custom.map(|c| row("", c)))
            .child(row("Antialias", aa))
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
            .child(row("Caption", caption));

        // -- Preview and the buttons
        let ratio = h as f32 / w.max(1) as f32;
        let ph = (PREVIEW_WIDTH * ratio).min(PREVIEW_WIDTH * 1.2);
        let pw = ph / ratio;
        let shown = self.movie.preview.get(self.movie.preview_at).cloned();
        let screen = div()
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
            .when_some(shown, |d, image| {
                // What is showing is out of date while the next is drawn.
                let stale = self.movie.previewing;
                d.child(gpui::img(gpui::ImageSource::Render(image)).size_full().when(stale, |i| i.opacity(0.5)))
            })
            .when(self.movie.preview.is_empty(), |d| {
                d.child(
                    div()
                        .text_size(px(10.5))
                        .text_color(theme::faint())
                        .child(if self.movie.previewing { "Drawing the preview…" } else { "Nothing to preview" }),
                )
            });
        let at = (!self.movie.preview.is_empty()).then(|| {
            let shot = self.movie.preview_at * self.movie.preview_stride.max(1);
            let updating = if self.movie.previewing { " · updating…" } else { "" };
            note(format!("{} / {} · {}{updating}", shot + 1, count, clock(shot as f32 / o.fps.max(1) as f32)))
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
                let (done, total) = (p.done.load(Ordering::Relaxed), p.total.load(Ordering::Relaxed).max(1));
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
                            .child(div().h_full().w(relative(f)).rounded(px(3.)).bg(theme::ink())),
                    )
                    .child(
                        div()
                            .flex()
                            .justify_between()
                            .items_center()
                            .child(note(format!("Frame {done} of {total}")))
                            .child(plain("movie-cancel", "Cancel").on_click(cx.listener(|this, _, _, cx| this.cancel_movie(cx)))),
                    )
                    .into_any_element()
            }
            None => div()
                .flex()
                .gap(px(6.))
                .justify_end()
                .child(plain("movie-close", "Close").on_click(cx.listener(|this, _, _, cx| this.close_movie(cx))))
                .child(primary("movie-export", "Export…").on_click(cx.listener(|this, _, _, cx| this.choose_movie_path(cx))))
                .into_any_element(),
        };
        let right = div()
            .w(px(PREVIEW_WIDTH))
            .flex_none()
            .flex()
            .flex_col()
            .gap(px(8.))
            .pt(px(6.))
            .child(screen)
            .children(at)
            .child(div().text_size(px(10.5)).text_color(theme::text()).child(SharedString::from(summary)))
            .when(o.format == Format::Gif && o.fps > 50, |d| d.child(note("GIFs play at most 50 fps".to_string())))
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
                .bg(gpui::Hsla { a: 0.35, ..theme::ink() })
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
                        .child(div().text_size(px(13.)).text_color(theme::ink()).child("Export a film"))
                        .child(div().flex().gap(px(22.)).child(left).child(right)),
                )
                .into_any_element(),
        )
    }

    /// `movie <file> [key=value …]`: make a film without the dialog.
    pub(super) fn movie_command(&mut self, rest: &str, cx: &mut Context<Self>) -> anyhow::Result<String> {
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
        Ok(format!("making {} ({n} frames, {})", path.display(), clock(secs)))
    }
}
