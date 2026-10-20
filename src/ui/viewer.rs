//! The 3D viewport.
//!
//! A GPUI `canvas` that draws the molecule offscreen, on the graphics card or
//! failing that the CPU, and uploads the result as an image. The frame is
//! cached against the camera and viewport size, so idle redraws cost nothing.
//!
//! Frames are drawn off the main thread. While one is being drawn the last
//! one stays up, and the newest view is drawn once it is done, so the window
//! keeps answering however long a frame takes.

use std::cell::RefCell;
use std::collections::VecDeque;
use std::rc::Rc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

use gpui::{canvas, px, Corners, IntoElement, RenderImage, Styled};
use image::{Frame, RgbaImage};
use smallvec::smallvec;

use crate::render::{self, Camera, Scene, Style};
use crate::theme;

#[derive(Clone, PartialEq, Eq)]
struct Key {
    width: usize,
    height: usize,
    /// Camera parameters, compared bitwise.
    cam: [u32; 15],
    /// The scene's own stamp. Its address will not do: a freed scene's
    /// memory is reused, and a recoloured scene has the same number of
    /// primitives, so the old frame would be shown again.
    scene: u64,
    style: Style,
}

fn key(scene: &Arc<Scene>, cam: &Camera, width: usize, height: usize, style: Style) -> Key {
    Key {
        width,
        height,
        cam: [
            cam.orientation[0][0].to_bits(),
            cam.orientation[0][1].to_bits(),
            cam.orientation[0][2].to_bits(),
            cam.orientation[1][0].to_bits(),
            cam.orientation[1][1].to_bits(),
            cam.orientation[1][2].to_bits(),
            cam.orientation[2][0].to_bits(),
            cam.orientation[2][1].to_bits(),
            cam.orientation[2][2].to_bits(),
            cam.zoom.to_bits(),
            cam.pan[0].to_bits(),
            cam.pan[1].to_bits(),
            cam.perspective as u32,
            cam.slab.map_or(0, |s| s[0].to_bits()),
            cam.slab.map_or(0, |s| s[1].to_bits()),
        ],
        scene: scene.version,
        style,
    }
}

#[derive(Default)]
pub struct Cache {
    key: Option<Key>,
    image: Option<Arc<RenderImage>>,
    /// The view of the frame being drawn, if one is.
    drawing: Option<Key>,
    /// The pixel size last painted, which is what an exported image matches.
    pub size: (usize, usize),
    /// Where the viewport was last painted in the window.
    pub bounds: gpui::Bounds<gpui::Pixels>,
    /// When each frame of the last second went up.
    shown: VecDeque<Instant>,
    /// Whether a repaint is due for when the frame rate falls idle.
    settling: bool,
    /// The frame rate over `shown`, while there were two frames to time.
    rate: f32,
    /// The view the frame up now was ray traced at, if it was.
    traced: Option<Key>,
    /// A trace of the still view, under way.
    tracing: Option<Tracing>,
}

struct Tracing {
    key: Key,
    /// Set to give it up, as soon as the view moves.
    cancel: Arc<AtomicBool>,
    /// Rows traced, out of `rows`.
    done: Arc<AtomicUsize>,
    rows: usize,
}

/// How long the view must hold still before it is traced, so a pause in a
/// drag does not start one.
const STILL: Duration = Duration::from_millis(250);

/// How far back the frame rate looks, and how long the view must be still
/// before it is called idle.
const WINDOW: Duration = Duration::from_secs(1);

impl Cache {
    /// Frames put up per second, as last measured while the view moved,
    /// and whether it is moving still. Frames are drawn only when the view
    /// changes, so a still view keeps the rate it last had.
    /// How far along the trace of the still view is, while one is.
    pub fn tracing(&self) -> Option<f32> {
        let t = self.tracing.as_ref()?;
        Some(t.done.load(Ordering::Relaxed) as f32 / t.rows.max(1) as f32)
    }

    pub fn fps(&self) -> (f32, bool) {
        let live = self.shown.back().is_some_and(|t| t.elapsed() <= WINDOW);
        (self.rate, live)
    }

    /// Note a new frame, and see that the rate reads idle once the view
    /// stops: nothing else would paint the window again to say so.
    fn shown(this: &Rc<RefCell<Cache>>, window: &mut gpui::Window, cx: &mut gpui::App) {
        let mut cache = this.borrow_mut();
        let now = Instant::now();
        cache.shown.push_back(now);
        while cache.shown.front().is_some_and(|t| now - *t > WINDOW) {
            cache.shown.pop_front();
        }
        if let (Some(first), true) = (cache.shown.front(), cache.shown.len() > 1) {
            let span = (now - *first).as_secs_f32();
            if span > 0.0 {
                cache.rate = (cache.shown.len() - 1) as f32 / span;
            }
        }
        if std::mem::replace(&mut cache.settling, true) {
            return;
        }
        let this = this.clone();
        window
            .spawn(cx, async move |cx| {
                loop {
                    let last = this.borrow().shown.back().copied();
                    let Some(wait) = last.and_then(|t| (t + WINDOW).checked_duration_since(Instant::now()))
                    else {
                        break;
                    };
                    cx.background_executor().timer(wait + Duration::from_millis(10)).await;
                }
                this.borrow_mut().settling = false;
                let _ = cx.update(|window, _| window.refresh());
            })
            .detach();
    }
}

/// `interactive` renders at 1x instead of the display's scale factor, which
/// keeps dragging smooth on a retina screen when the frame is drawn on the
/// CPU. The graphics card keeps up at full resolution.
pub fn viewport(
    scene: Arc<Scene>,
    cam: Camera,
    cache: Rc<RefCell<Cache>>,
    interactive: bool,
    style: Style,
    traced: bool,
) -> impl IntoElement {
    canvas(
        |_, _, _| {},
        move |bounds, _, window, cx| {
            let scale = if interactive && !render::accelerated() {
                1.0
            } else {
                window.scale_factor()
            };
            let width = (f32::from(bounds.size.width) * scale).round().max(1.0) as usize;
            let height = (f32::from(bounds.size.height) * scale).round().max(1.0) as usize;

            let _timer = Phase::Paint.timer();
            let mut current = cache.borrow_mut();
            current.size = (width, height);
            current.bounds = bounds;
            let k = key(&scene, &cam, width, height, style);
            let resized = current
                .key
                .as_ref()
                .is_none_or(|old| (old.width, old.height) != (width, height));
            if resized {
                // Nothing of the right size to show meanwhile: draw it now.
                let image = draw(&scene, &cam, width, height, style);
                if let Some(old) = current.image.replace(image) {
                    let _ = window.drop_image(old);
                }
                current.key = Some(k.clone());
                drop(current);
                Cache::shown(&cache, window, cx);
                current = cache.borrow_mut();
            } else if current.key.as_ref() != Some(&k) && current.drawing.is_none() {
                current.drawing = Some(k.clone());
                let (scene, cache) = (scene.clone(), cache.clone());
                let drawn = cx
                    .background_executor()
                    .spawn(async move { draw(&scene, &cam, width, height, style) });
                window
                    .spawn(cx, async move |cx| {
                        let image = drawn.await;
                        let _ = cx.update(|window, cx| {
                            let this = cache.clone();
                            let mut cache = cache.borrow_mut();
                            // A resize may have drawn past it in the meantime.
                            let Some(k) = cache.drawing.take() else { return };
                            if cache.key.as_ref().is_some_and(|old| {
                                (old.width, old.height) != (k.width, k.height)
                            }) {
                                let _ = window.drop_image(image);
                            } else {
                                if let Some(old) = cache.image.replace(image) {
                                    let _ = window.drop_image(old);
                                }
                                cache.key = Some(k);
                                drop(cache);
                                Cache::shown(&this, window, cx);
                            }
                            // Paint again, which starts on the newest view if
                            // it has moved on since.
                            window.refresh();
                        });
                    })
                    .detach();
            }

            // Ray traced once the view holds still: the quick frame stays up
            // until the trace replaces it, and any move gives the trace up.
            if let Some(t) = current.tracing.as_ref().filter(|t| !traced || t.key != k) {
                t.cancel.store(true, Ordering::Relaxed);
                current.tracing = None;
            }
            if !traced && current.traced.take().is_some_and(|t| t == k) {
                // The traced frame is up for this very view; draw it plainly.
                let image = draw(&scene, &cam, width, height, style);
                if let Some(old) = current.image.replace(image) {
                    let _ = window.drop_image(old);
                }
            }
            let quick_is_up = current.key.as_ref() == Some(&k) && current.drawing.is_none();
            if traced && quick_is_up && current.traced.as_ref() != Some(&k) && current.tracing.is_none() {
                let cancel = Arc::new(AtomicBool::new(false));
                let done = Arc::new(AtomicUsize::new(0));
                current.tracing = Some(Tracing {
                    key: k.clone(),
                    cancel: cancel.clone(),
                    done: done.clone(),
                    rows: crate::trace::rows(height, width),
                });
                let (scene, cache) = (scene.clone(), cache.clone());
                window
                    .spawn(cx, async move |cx| {
                        cx.background_executor().timer(STILL).await;
                        if cancel.load(Ordering::Relaxed) {
                            return;
                        }
                        let finished = Arc::new(AtomicBool::new(false));
                        let drawn = cx.background_executor().spawn({
                            let (cancel, done, finished) = (cancel.clone(), done.clone(), finished.clone());
                            async move {
                                let control = crate::trace::Control { cancel: Some(&cancel), done: Some(&done) };
                                let frame = crate::trace::render(&scene, &cam, width, height, theme::view(), style, &control);
                                finished.store(true, Ordering::Relaxed);
                                frame.map(image_of)
                            }
                        });
                        // Painted again now and then, so the progress shows.
                        while !finished.load(Ordering::Relaxed) {
                            let _ = cx.update(|window, _| window.refresh());
                            cx.background_executor().timer(Duration::from_millis(150)).await;
                        }
                        let image = drawn.await;
                        let _ = cx.update(|window, _| {
                            let mut cache = cache.borrow_mut();
                            let current = cache.tracing.as_ref().is_some_and(|t| Arc::ptr_eq(&t.cancel, &cancel));
                            match image {
                                Some(image) if current && !cancel.load(Ordering::Relaxed) => {
                                    if let Some(old) = cache.image.replace(image) {
                                        let _ = window.drop_image(old);
                                    }
                                    cache.traced = cache.tracing.take().map(|t| t.key);
                                }
                                Some(image) => {
                                    let _ = window.drop_image(image);
                                }
                                None => {}
                            }
                            if current {
                                cache.tracing = None;
                            }
                            window.refresh();
                        });
                    })
                    .detach();
            }

            if let Some(image) = current.image.clone() {
                let _ = window.paint_image(bounds, Corners::all(px(0.)), image, 0, false);
            }
        },
    )
    .size_full()
}

fn draw(scene: &Scene, cam: &Camera, width: usize, height: usize, style: Style) -> Arc<RenderImage> {
    let frame = {
        let _timer = Phase::Draw.timer();
        render::render(scene, cam, width, height, theme::view(), style)
    };
    image_of(frame)
}

fn image_of(frame: render::Frame) -> Arc<RenderImage> {
    let buffer = RgbaImage::from_raw(frame.width as u32, frame.height as u32, frame.bgra)
        .expect("frame buffer is width * height * 4");
    Arc::new(RenderImage::new(smallvec![Frame::new(buffer)]))
}

/// `PROFOLD_SPIN=1` turns the camera every frame and prints where the time
/// goes, for measuring the viewer as it is used rather than in isolation.
pub fn spinning() -> bool {
    static SPIN: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *SPIN.get_or_init(|| std::env::var_os("PROFOLD_SPIN").is_some())
}

#[derive(Clone, Copy)]
pub enum Phase {
    /// Building the window's element tree.
    Tree,
    /// Painting the viewport, the draw included.
    Paint,
    /// Drawing the frame itself.
    Draw,
}

/// Time spent in each phase, in nanoseconds, and frames shown.
static SPENT: [std::sync::atomic::AtomicU64; 3] = [const { std::sync::atomic::AtomicU64::new(0) }; 3];
static DRAWN: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);

thread_local! {
    static TRACE: std::cell::RefCell<(u32, Option<std::time::Instant>)> = Default::default();
}

pub struct Timer(Phase, std::time::Instant);

impl Phase {
    pub fn timer(self) -> Option<Timer> {
        spinning().then(|| Timer(self, std::time::Instant::now()))
    }
}

impl Drop for Timer {
    fn drop(&mut self) {
        use std::sync::atomic::Ordering::Relaxed;
        let spent = self.1.elapsed();
        SPENT[self.0 as usize].fetch_add(spent.as_nanos() as u64, Relaxed);
        match self.0 {
            Phase::Draw => {
                DRAWN.fetch_add(1, Relaxed);
                return;
            }
            Phase::Tree => return,
            Phase::Paint => {}
        }
        TRACE.with_borrow_mut(|(frames, since)| {
            *frames += 1;
            let start = *since.get_or_insert_with(std::time::Instant::now);
            if *frames == 120 {
                let n = *frames;
                let wall = start.elapsed();
                let drawn = DRAWN.swap(0, Relaxed).max(1);
                let per = |i: usize, n: u32| {
                    std::time::Duration::from_nanos(SPENT[i].swap(0, Relaxed) / n as u64)
                };
                eprintln!(
                    "{:.1} fps  {:.1} drawn/s  tree {:.2?}  paint {:.2?}  draw {:.2?}",
                    n as f64 / wall.as_secs_f64(),
                    drawn as f64 / wall.as_secs_f64(),
                    per(0, n),
                    per(1, n),
                    per(2, drawn),
                );
                *frames = 0;
                *since = Some(std::time::Instant::now());
            }
        });
    }
}
