//! Animations written out as films: a trajectory played through, the camera
//! turning or rocking, a path flown by hand and recorded, or all of these at
//! once, at any size and quality, on any background, with a frame counter,
//! simulated time, a scale bar and a caption laid over it.
//!
//! The work is split in three. `plan` turns the options into shots, each a
//! camera and a place along the trajectory. `Job` draws a shot: it poses the
//! structure there, blending between frames where asked, draws it larger
//! than wanted and scales it down to smooth the edges, puts it on its
//! background and writes the overlays. A `Sink` takes the drawn frames into
//! a file: H.264 in MP4 and GIF are made here, numbered PNGs too, and WebM
//! and ProRes, which keep transparency in a film, are handed to ffmpeg when
//! it is installed.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;

use anyhow::{anyhow, bail, Context as _, Result};
use serde_json::{json, Value};

use crate::protein::{Cell, Structure};
use crate::render::{self, Camera, Scene, Style};

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------

/// How the camera moves while the film runs.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Motion {
    Still,
    /// Round through `degrees`, a full turn by default.
    Turn,
    /// Back and forth, `degrees` either side of where it started.
    Rock,
    /// Along the path recorded in the viewport.
    Recorded,
}

/// The screen axis the camera turns or rocks about.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Axis {
    /// The screen's upright: the molecule turns as on a turntable.
    Vertical,
    /// Across the screen: it tumbles towards the viewer.
    Horizontal,
    /// The line of sight: it rolls in the plane of the screen.
    View,
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Background {
    /// The page as the app draws it, light or dark.
    Theme,
    Color([f32; 3]),
    /// Nothing at all, where the format keeps transparency.
    Transparent,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Format {
    /// H.264, through ffmpeg.
    Mp4,
    Gif,
    /// Numbered PNGs in a folder, lossless and with transparency.
    Png,
    /// VP9 with transparency, through ffmpeg.
    WebM,
    /// ProRes 4444 with transparency, through ffmpeg.
    Mov,
}

impl Format {
    pub const ALL: [Format; 5] = [
        Format::Mp4,
        Format::Gif,
        Format::Png,
        Format::WebM,
        Format::Mov,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Format::Mp4 => "MP4",
            Format::Gif => "GIF",
            Format::Png => "PNG frames",
            Format::WebM => "WebM",
            Format::Mov => "ProRes",
        }
    }

    pub fn extension(self) -> &'static str {
        match self {
            Format::Mp4 => "mp4",
            Format::Gif => "gif",
            Format::Png => "",
            Format::WebM => "webm",
            Format::Mov => "mov",
        }
    }

    pub fn from_extension(ext: &str) -> Option<Format> {
        Some(match ext.to_ascii_lowercase().as_str() {
            "mp4" | "m4v" => Format::Mp4,
            "gif" => Format::Gif,
            "png" | "" => Format::Png,
            "webm" => Format::WebM,
            "mov" => Format::Mov,
            _ => return None,
        })
    }

    /// Whether a transparent background survives in it. A GIF's pixels are
    /// either there or not, so soft edges go hard.
    pub fn keeps_alpha(self) -> bool {
        !matches!(self, Format::Mp4)
    }

    /// Made with ffmpeg rather than here.
    pub fn needs_ffmpeg(self) -> bool {
        matches!(self, Format::Mp4 | Format::WebM | Format::Mov)
    }

    /// Whether a quality setting means anything to it.
    pub fn lossy(self) -> bool {
        matches!(self, Format::Mp4 | Format::WebM)
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Size {
    /// The viewport's own size.
    Window,
    Hd720,
    Hd1080,
    Uhd,
    /// 1080 a side.
    Square,
    Custom,
}

impl Size {
    pub const ALL: [Size; 6] = [
        Size::Window,
        Size::Hd720,
        Size::Hd1080,
        Size::Uhd,
        Size::Square,
        Size::Custom,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Size::Window => "Window",
            Size::Hd720 => "720p",
            Size::Hd1080 => "1080p",
            Size::Uhd => "4K",
            Size::Square => "Square",
            Size::Custom => "Custom",
        }
    }
}

/// A corner of the picture, for the legend.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Corner {
    TopLeft,
    TopRight,
    BottomLeft,
    BottomRight,
}

impl Corner {
    pub const ALL: [Corner; 4] = [
        Corner::BottomRight,
        Corner::BottomLeft,
        Corner::TopRight,
        Corner::TopLeft,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Corner::TopLeft => "Top left",
            Corner::TopRight => "Top right",
            Corner::BottomLeft => "Bottom left",
            Corner::BottomRight => "Bottom right",
        }
    }

    fn index(self) -> usize {
        self as usize
    }

    fn right(self) -> bool {
        matches!(self, Corner::TopRight | Corner::BottomRight)
    }

    fn bottom(self) -> bool {
        matches!(self, Corner::BottomLeft | Corner::BottomRight)
    }
}

/// What the colours mean, as the panel's key shows it, for writing into a
/// film.
#[derive(Clone, Debug, PartialEq)]
pub enum Legend {
    /// A continuous scale, with what each end means.
    Ramp {
        title: String,
        colors: Vec<[f32; 3]>,
        low: String,
        high: String,
    },
    /// A colour each, and how many more were left out.
    Swatches {
        title: String,
        items: Vec<([f32; 3], String)>,
        more: usize,
    },
}

/// The unit a simulation's timestep is given in.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TimeUnit {
    Fs,
    Ps,
    /// Reduced (Lennard-Jones) time.
    Tau,
}

/// Everything a film is made with. Kept between launches.
#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    /// Play the trajectory's frames through, rather than holding the one
    /// shown.
    pub trajectory: bool,
    pub motion: Motion,
    /// Through how far a turn goes, or how far either side a rock swings.
    pub degrees: f32,
    pub axis: Axis,
    /// The camera ends where it began, so the film loops without a jump.
    pub seamless: bool,
    /// The trajectory's frames used, from 1 as the transport counts them:
    /// first, last (0 for the end) and every how many.
    pub first: usize,
    pub last: usize,
    pub every: usize,
    pub fps: u32,
    /// How long the film runs, in seconds; 0 for as long as the frames
    /// take at `fps`.
    pub seconds: f32,
    /// Frames made between each pair of the trajectory's, by blending the
    /// two, so a sparse dump plays smoothly.
    pub smooth: u32,
    /// Forward then back again.
    pub pingpong: bool,
    pub size: Size,
    pub width: u32,
    pub height: u32,
    /// Each frame is drawn this many times larger, then scaled down.
    pub supersample: u32,
    /// Each frame is ray traced, with shadows and ambient occlusion, and
    /// smoothed by the tracer's own oversampling rather than `supersample`.
    pub traced: bool,
    /// 1–100, for formats that lose detail to stay small.
    pub quality: u32,
    pub background: Background,
    pub format: Format,
    pub counter: bool,
    pub time: bool,
    /// The simulation's timestep, for reading steps as time.
    pub dt: f32,
    pub unit: TimeUnit,
    pub scale_bar: bool,
    /// The colour key, in a corner of its own.
    pub legend: bool,
    pub corner: Corner,
    pub caption: String,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            trajectory: true,
            motion: Motion::Still,
            degrees: 360.0,
            axis: Axis::Vertical,
            seamless: true,
            first: 1,
            last: 0,
            every: 1,
            fps: 30,
            seconds: 0.0,
            smooth: 0,
            pingpong: false,
            size: Size::Hd1080,
            width: 1920,
            height: 1080,
            supersample: 2,
            traced: false,
            quality: 75,
            background: Background::Theme,
            format: Format::Mp4,
            counter: false,
            time: false,
            dt: 1.0,
            unit: TimeUnit::Fs,
            scale_bar: false,
            legend: false,
            corner: Corner::BottomRight,
            caption: String::new(),
        }
    }
}

fn name<T: Copy + PartialEq>(table: &[(T, &'static str)], v: T) -> &'static str {
    table
        .iter()
        .find(|(t, _)| *t == v)
        .map(|(_, n)| *n)
        .unwrap_or("")
}

fn named<T: Copy>(table: &[(T, &'static str)], s: &str) -> Option<T> {
    table
        .iter()
        .find(|(_, n)| n.eq_ignore_ascii_case(s))
        .map(|(t, _)| *t)
}

const MOTIONS: [(Motion, &str); 4] = [
    (Motion::Still, "still"),
    (Motion::Turn, "turn"),
    (Motion::Rock, "rock"),
    (Motion::Recorded, "recorded"),
];
const AXES: [(Axis, &str); 3] = [
    (Axis::Vertical, "y"),
    (Axis::Horizontal, "x"),
    (Axis::View, "z"),
];
const FORMATS: [(Format, &str); 5] = [
    (Format::Mp4, "mp4"),
    (Format::Gif, "gif"),
    (Format::Png, "png"),
    (Format::WebM, "webm"),
    (Format::Mov, "mov"),
];
const SIZES: [(Size, &str); 6] = [
    (Size::Window, "window"),
    (Size::Hd720, "720p"),
    (Size::Hd1080, "1080p"),
    (Size::Uhd, "4k"),
    (Size::Square, "square"),
    (Size::Custom, "custom"),
];
const CORNERS: [(Corner, &str); 4] = [
    (Corner::TopLeft, "top-left"),
    (Corner::TopRight, "top-right"),
    (Corner::BottomLeft, "bottom-left"),
    (Corner::BottomRight, "bottom-right"),
];
const UNITS: [(TimeUnit, &str); 3] = [
    (TimeUnit::Fs, "fs"),
    (TimeUnit::Ps, "ps"),
    (TimeUnit::Tau, "tau"),
];

impl Options {
    pub fn to_json(&self) -> String {
        let background = match self.background {
            Background::Theme => json!("theme"),
            Background::Transparent => json!("transparent"),
            Background::Color(c) => json!(hex(c)),
        };
        serde_json::to_string_pretty(&json!({
            "trajectory": self.trajectory,
            "motion": name(&MOTIONS, self.motion),
            "degrees": self.degrees,
            "axis": name(&AXES, self.axis),
            "seamless": self.seamless,
            "every": self.every,
            "fps": self.fps,
            "seconds": self.seconds,
            "smooth": self.smooth,
            "pingpong": self.pingpong,
            "size": name(&SIZES, self.size),
            "width": self.width,
            "height": self.height,
            "supersample": self.supersample,
            "traced": self.traced,
            "quality": self.quality,
            "background": background,
            "format": name(&FORMATS, self.format),
            "counter": self.counter,
            "time": self.time,
            "dt": self.dt,
            "unit": name(&UNITS, self.unit),
            "scale_bar": self.scale_bar,
            "legend": self.legend,
            "corner": name(&CORNERS, self.corner),
            "caption": self.caption,
        }))
        .unwrap_or_default()
    }

    /// What was kept, over the defaults. The frame range is not kept: it
    /// belongs to the trajectory it was set for.
    pub fn from_json(v: &Value) -> Options {
        let mut o = Options::default();
        let b = |k: &str, d: bool| v.get(k).and_then(Value::as_bool).unwrap_or(d);
        let f = |k: &str, d: f32| {
            v.get(k)
                .and_then(Value::as_f64)
                .map(|x| x as f32)
                .unwrap_or(d)
        };
        let u = |k: &str, d: u32| {
            v.get(k)
                .and_then(Value::as_u64)
                .map(|x| x as u32)
                .unwrap_or(d)
        };
        let s = |k: &str| v.get(k).and_then(Value::as_str).unwrap_or("");
        o.trajectory = b("trajectory", o.trajectory);
        // A recording is gone at the next launch.
        o.motion = named(&MOTIONS, s("motion"))
            .filter(|m| *m != Motion::Recorded)
            .unwrap_or(o.motion);
        o.degrees = f("degrees", o.degrees);
        o.axis = named(&AXES, s("axis")).unwrap_or(o.axis);
        o.seamless = b("seamless", o.seamless);
        o.every = u("every", o.every as u32).max(1) as usize;
        o.fps = u("fps", o.fps);
        o.seconds = f("seconds", o.seconds);
        o.smooth = u("smooth", o.smooth);
        o.pingpong = b("pingpong", o.pingpong);
        o.size = named(&SIZES, s("size")).unwrap_or(o.size);
        o.width = u("width", o.width);
        o.height = u("height", o.height);
        o.supersample = u("supersample", o.supersample);
        o.traced = b("traced", o.traced);
        o.quality = u("quality", o.quality);
        o.background = match s("background") {
            "theme" => Background::Theme,
            "transparent" => Background::Transparent,
            other => crate::settings::parse_hex(other)
                .map(Background::Color)
                .unwrap_or(o.background),
        };
        o.format = named(&FORMATS, s("format")).unwrap_or(o.format);
        o.counter = b("counter", o.counter);
        o.time = b("time", o.time);
        o.dt = f("dt", o.dt);
        o.unit = named(&UNITS, s("unit")).unwrap_or(o.unit);
        o.scale_bar = b("scale_bar", o.scale_bar);
        o.legend = b("legend", o.legend);
        o.corner = named(&CORNERS, s("corner")).unwrap_or(o.corner);
        o.caption = s("caption").to_string();
        o.clamped()
    }

    /// The same with every number brought into a range that makes sense.
    pub fn clamped(mut self) -> Options {
        self.degrees = self.degrees.clamp(1.0, 3600.0);
        self.every = self.every.max(1);
        self.first = self.first.max(1);
        self.fps = self.fps.clamp(1, 120);
        self.seconds = self.seconds.clamp(0.0, 3600.0);
        self.smooth = self.smooth.min(30);
        self.width = self.width.clamp(16, 7680);
        self.height = self.height.clamp(16, 7680);
        self.supersample = match self.supersample {
            0 | 1 => 1,
            2 | 3 => 2,
            _ => 4,
        };
        self.quality = self.quality.clamp(1, 100);
        if !self.dt.is_finite() || self.dt <= 0.0 {
            self.dt = 1.0;
        }
        // A film with no transparency in it is drawn on the page.
        if !self.format.keeps_alpha() && self.background == Background::Transparent {
            self.background = Background::Theme;
        }
        self
    }

    /// The film's size in pixels, for a viewport `view` across. Both even,
    /// as H.264 needs.
    pub fn output_size(&self, view: (usize, usize)) -> (usize, usize) {
        let (w, h) = match self.size {
            Size::Window => (view.0.max(16), view.1.max(16)),
            Size::Hd720 => (1280, 720),
            Size::Hd1080 => (1920, 1080),
            Size::Uhd => (3840, 2160),
            Size::Square => (1080, 1080),
            Size::Custom => (self.width as usize, self.height as usize),
        };
        (w & !1, h & !1)
    }

    /// The trajectory's frames used, as indices from 0.
    pub fn keys(&self, count: usize, now: usize) -> Vec<usize> {
        if !self.trajectory || count < 2 {
            return vec![now.min(count.saturating_sub(1))];
        }
        let first = (self.first.max(1) - 1).min(count - 1);
        let last = if self.last == 0 {
            count - 1
        } else {
            (self.last - 1).min(count - 1)
        };
        let (first, last) = if last < first {
            (last, first)
        } else {
            (first, last)
        };
        (first..=last).step_by(self.every.max(1)).collect()
    }
}

fn hex(c: [f32; 3]) -> String {
    let b = c.map(|v| (v.clamp(0.0, 1.0) * 255.0).round() as u8);
    format!("#{:02x}{:02x}{:02x}", b[0], b[1], b[2])
}

// ---------------------------------------------------------------------------
// Recording
// ---------------------------------------------------------------------------

/// One moment of a recording: when, the camera then, and the trajectory's
/// frame shown.
#[derive(Clone, Copy, Debug)]
pub struct Sample {
    pub t: f32,
    pub camera: Camera,
    pub frame: usize,
}

/// How often the viewport is sampled while recording.
pub const RECORD_HZ: f32 = 30.0;

// ---------------------------------------------------------------------------
// Planning
// ---------------------------------------------------------------------------

/// One frame of the film: the camera, and where along the trajectory,
/// `t` of the way from frame `from` to frame `to`.
#[derive(Clone, Copy, Debug)]
pub struct Shot {
    pub camera: Camera,
    pub from: usize,
    pub to: usize,
    pub t: f32,
}

impl Shot {
    /// The trajectory frame nearest the shot.
    pub fn nearest(&self) -> usize {
        if self.t < 0.5 {
            self.from
        } else {
            self.to
        }
    }
}

/// Seconds a film that only moves the camera runs, when not told.
const CAMERA_SECONDS: f32 = 6.0;

/// Turn `base` through `angle` radians about `axis`.
pub fn turned(base: &Camera, axis: Axis, angle: f32) -> Camera {
    let mut cam = *base;
    match axis {
        Axis::Vertical => cam.rotate(angle, 0.0),
        Axis::Horizontal => cam.rotate(0.0, angle),
        Axis::View => {
            let (s, c) = angle.sin_cos();
            let roll = [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]];
            cam.orientation = render::orthonormal(render::mul(roll, cam.orientation));
        }
    }
    cam
}

/// Somewhere between two cameras.
fn between(a: &Camera, b: &Camera, t: f32) -> Camera {
    let mix = |x: f32, y: f32| x + (y - x) * t;
    let mut o = a.orientation;
    for (r, row) in o.iter_mut().enumerate() {
        for (k, v) in row.iter_mut().enumerate() {
            *v = mix(a.orientation[r][k], b.orientation[r][k]);
        }
    }
    Camera {
        orientation: render::orthonormal(o),
        zoom: mix(a.zoom, b.zoom),
        pan: [mix(a.pan[0], b.pan[0]), mix(a.pan[1], b.pan[1])],
        center: std::array::from_fn(|k| mix(a.center[k], b.center[k])),
        radius: mix(a.radius, b.radius),
        perspective: if t < 0.5 {
            a.perspective
        } else {
            b.perspective
        },
        slab: if t < 0.5 { a.slab } else { b.slab },
    }
}

/// The film's shots, from the options, the camera as it is, the trajectory
/// frame shown and how many there are, and a recording if there is one.
pub fn plan(
    o: &Options,
    base: &Camera,
    now: usize,
    count: usize,
    recording: Option<&[Sample]>,
) -> Vec<Shot> {
    let fps = o.fps.max(1) as f32;
    if let (Motion::Recorded, Some(rec)) = (o.motion, recording.filter(|r| !r.is_empty())) {
        return plan_recorded(o, rec, fps);
    }
    let keys = o.keys(count, now);
    let moving = keys.len() > 1;
    let natural = if moving {
        (keys.len() - 1) * (o.smooth as usize + 1) + 1
    } else {
        (CAMERA_SECONDS * fps) as usize
    };
    // Ping-pong goes there and back, the ends once each.
    let total = if o.seconds > 0.0 {
        (o.seconds * fps).round().max(1.0) as usize
    } else if o.pingpong {
        (2 * natural).saturating_sub(2).max(1)
    } else {
        natural
    };
    let span = |i: usize| -> f32 {
        if total < 2 {
            0.0
        } else if o.pingpong {
            let x = i as f32 / total as f32;
            1.0 - (1.0 - 2.0 * x).abs()
        } else {
            i as f32 / (total - 1) as f32
        }
    };
    (0..total)
        .map(|i| {
            // Along the trajectory.
            let u = span(i) * (keys.len() - 1) as f32;
            let a = (u.floor() as usize).min(keys.len() - 1);
            let b = (a + 1).min(keys.len() - 1);
            let frac = u - a as f32;
            let (from, to, t) = if o.smooth > 0 {
                (keys[a], keys[b], frac)
            } else {
                let k = keys[if frac < 0.5 { a } else { b }];
                (k, k, 0.0)
            };
            // And the camera, which keeps going round however the frames go.
            let c = if total < 2 {
                0.0
            } else if o.seamless || o.motion == Motion::Rock {
                i as f32 / total as f32
            } else {
                i as f32 / (total - 1) as f32
            };
            let camera = match o.motion {
                Motion::Turn => turned(base, o.axis, o.degrees.to_radians() * c),
                Motion::Rock => turned(
                    base,
                    o.axis,
                    o.degrees.to_radians() * (std::f32::consts::TAU * c).sin(),
                ),
                Motion::Still | Motion::Recorded => *base,
            };
            Shot {
                camera,
                from,
                to,
                t,
            }
        })
        .collect()
}

fn plan_recorded(o: &Options, rec: &[Sample], fps: f32) -> Vec<Shot> {
    let start = rec[0].t;
    let length = (rec[rec.len() - 1].t - start).max(0.0);
    let total = ((length * fps).round() as usize + 1).max(1);
    let mut j = 0;
    (0..total)
        .map(|i| {
            let t = start + i as f32 / fps;
            while j + 1 < rec.len() - 1 && rec[j + 1].t <= t {
                j += 1;
            }
            let (a, b) = (&rec[j], &rec[(j + 1).min(rec.len() - 1)]);
            let gap = b.t - a.t;
            let f = if gap > 1e-6 {
                ((t - a.t) / gap).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let camera = between(&a.camera, &b.camera, f);
            let (from, to, t) = if o.smooth > 0 && a.frame != b.frame {
                (a.frame, b.frame, f)
            } else {
                let k = if f < 0.5 { a.frame } else { b.frame };
                (k, k, 0.0)
            };
            Shot {
                camera,
                from,
                to,
                t,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Posing
// ---------------------------------------------------------------------------

fn inverse(m: [[f32; 3]; 3]) -> Option<[[f32; 3]; 3]> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-12 {
        return None;
    }
    let d = 1.0 / det;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * d,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * d,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * d,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * d,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * d,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * d,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * d,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * d,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * d,
        ],
    ])
}

/// How far `b` is from `a`, the short way round the cell where there is
/// one: an atom that left one side and came back in at the other has moved
/// a little, not across the box.
fn displacement(a: [f32; 3], b: [f32; 3], cell: Option<&(Cell, [[f32; 3]; 3])>) -> [f32; 3] {
    let mut d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    if let Some((c, inv)) = cell {
        // Cartesian = vectors^T · fractional, so fractional = inv(vectors^T) · d,
        // with `inv` that inverse.
        let f: [f32; 3] =
            std::array::from_fn(|i| inv[i][0] * d[0] + inv[i][1] * d[1] + inv[i][2] * d[2]);
        let r = f.map(|x| x.round());
        if r != [0.0; 3] {
            for k in 0..3 {
                d[k] -= r[0] * c.vectors[0][k] + r[1] * c.vectors[1][k] + r[2] * c.vectors[2][k];
            }
        }
    }
    d
}

/// Put `s` where a shot is along its trajectory: at a frame, or blended
/// between two.
pub fn pose(s: &mut Structure, shot: &Shot) {
    if s.frame_count() < 2 {
        return;
    }
    if shot.from == shot.to || shot.t <= 1e-4 {
        s.set_frame(shot.from);
        return;
    }
    if shot.t >= 1.0 - 1e-4 {
        s.set_frame(shot.to);
        return;
    }
    let frames = s.frames.clone();
    let (Some(a), Some(b)) = (frames.get(shot.from), frames.get(shot.to)) else {
        return;
    };
    s.set_frame(shot.from);
    if a.pos.len() != s.atoms.len() || b.pos.len() != s.atoms.len() {
        return;
    }
    let t = shot.t;
    let cell = a.cell.and_then(|c| {
        let vt: [[f32; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| c.vectors[j][i]));
        Some((c, inverse(vt)?))
    });
    for (atom, (pa, pb)) in s.atoms.iter_mut().zip(a.pos.iter().zip(&b.pos)) {
        let d = displacement(*pa, *pb, cell.as_ref());
        atom.pos = [pa[0] + d[0] * t, pa[1] + d[1] * t, pa[2] + d[2] * t];
    }
    if let (Some(ca), Some(cb)) = (a.cell, b.cell) {
        let mix = |x: f32, y: f32| x + (y - x) * t;
        s.cell = Some(Cell {
            origin: std::array::from_fn(|k| mix(ca.origin[k], cb.origin[k])),
            vectors: std::array::from_fn(|i| {
                std::array::from_fn(|k| mix(ca.vectors[i][k], cb.vectors[i][k]))
            }),
        });
    }
}

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

/// What a frame of the film is drawn from, besides its shot.
pub struct Job {
    /// The structure as loaded, posed afresh for each shot.
    pub structure: Structure,
    /// How the app draws a structure posed, with all its representations,
    /// layers and colours.
    pub build: Arc<dyn Fn(&Structure) -> Scene + Send + Sync>,
    pub style: Style,
    /// The page the app draws on, for `Background::Theme`.
    pub page: [f32; 3],
    /// The viewport's size, which the camera's pan is measured against.
    pub view: (usize, usize),
    pub shots: Vec<Shot>,
    pub options: Options,
    /// The colour key, drawn when the options ask for it.
    pub legend: Option<Legend>,
    /// Text pinned to atoms, by atom index, following them frame by frame.
    pub labels: Vec<(u32, String)>,
}

/// A drawn frame, as straight (not premultiplied) RGBA.
pub struct Picture {
    pub width: usize,
    pub height: usize,
    pub rgba: Vec<u8>,
}

/// Neither side of a frame is drawn larger than this, however much it is
/// supersampled: graphics cards stop at 8192 or 16384.
const MAX_SIDE: usize = 8192;

/// What a job keeps from one frame to the next, so a scene is built only
/// when the structure has moved.
pub struct Drawer<'a> {
    job: &'a Job,
    posed: Structure,
    last: Option<(usize, usize, u32)>,
    scene: Option<Scene>,
    font: Option<ab_glyph::FontRef<'static>>,
}

static FONT: &[u8] = include_bytes!("../assets/fonts/Inter-Regular.ttf");

impl<'a> Drawer<'a> {
    pub fn new(job: &'a Job) -> Self {
        Drawer {
            job,
            posed: job.structure.clone(),
            last: None,
            scene: None,
            font: ab_glyph::FontRef::try_from_slice(FONT).ok(),
        }
    }

    /// Shot `i`, `width` by `height`, drawn `ss` times over.
    pub fn draw(&mut self, i: usize, width: usize, height: usize, ss: usize) -> Picture {
        let job = self.job;
        let o = &job.options;
        let shot = job.shots[i];
        let key = (shot.from, shot.to, (shot.t * 10_000.0) as u32);
        if self.last != Some(key) || self.scene.is_none() {
            pose(&mut self.posed, &shot);
            self.scene = Some((job.build)(&self.posed).stamp());
            self.last = Some(key);
        }
        let scene = self.scene.as_ref().expect("built above");

        if o.traced {
            let mut cam = shot.camera;
            let k = width.min(height) as f32 / job.view.0.min(job.view.1).max(1) as f32;
            cam.pan = [cam.pan[0] * k, cam.pan[1] * k];
            let traced =
                crate::trace::trace(scene, &cam, width, height, job.style, &Default::default())
                    .expect("not cancelled");
            let rgba = match o.background {
                Background::Transparent => traced.transparent(),
                Background::Theme | Background::Color(_) => {
                    let page = match o.background {
                        Background::Color(c) => c,
                        _ => job.page,
                    };
                    let mut rgba = traced.develop(page).bgra;
                    for p in rgba.chunks_exact_mut(4) {
                        p.swap(0, 2);
                    }
                    rgba
                }
            };
            let mut pic = Picture {
                width,
                height,
                rgba,
            };
            self.pin_labels(&mut pic, &cam);
            self.overlay(&mut pic, &shot, &cam, 1);
            return pic;
        }

        let mut ss = ss.max(1);
        while ss > 1 && (width * ss > MAX_SIDE || height * ss > MAX_SIDE) {
            ss /= 2;
        }
        let (rw, rh) = (width * ss, height * ss);
        // The pan is in the viewport's pixels; the picture is framed on its
        // shorter side, as the viewport is.
        let mut cam = shot.camera;
        let k = rw.min(rh) as f32 / job.view.0.min(job.view.1).max(1) as f32;
        cam.pan = [cam.pan[0] * k, cam.pan[1] * k];

        let big = match o.background {
            Background::Transparent => render::render_transparent(scene, &cam, rw, rh, job.style),
            Background::Theme | Background::Color(_) => {
                let page = match o.background {
                    Background::Color(c) => c,
                    _ => job.page,
                };
                let mut rgba = render::render(scene, &cam, rw, rh, page, job.style).bgra;
                for p in rgba.chunks_exact_mut(4) {
                    p.swap(0, 2);
                    p[3] = 255;
                }
                rgba
            }
        };
        let mut pic = Picture {
            width,
            height,
            rgba: shrink(&big, rw, rh, ss),
        };
        let mut small = cam;
        small.pan = [cam.pan[0] / ss as f32, cam.pan[1] / ss as f32];
        self.pin_labels(&mut pic, &small);
        self.overlay(&mut pic, &shot, &cam, ss);
        pic
    }

    /// The labels pinned to atoms, where the atoms are in this shot.
    fn pin_labels(&self, pic: &mut Picture, cam: &Camera) {
        let Some(font) = &self.font else { return };
        let pinned: Vec<([f32; 3], &str)> = self
            .job
            .labels
            .iter()
            .filter_map(|(a, t)| Some((self.posed.atoms.get(*a as usize)?.pos, t.as_str())))
            .collect();
        let size = (pic.width.min(pic.height) as f32 * 0.022).max(10.0);
        draw_labels(pic, font, cam, &pinned, size, self.ink(), self.job.page);
    }

    /// The text, the scale bar and the legend, over the frame. Each keeps
    /// to a corner, and what shares a corner is stacked inwards from it.
    fn overlay(&self, pic: &mut Picture, shot: &Shot, cam: &Camera, ss: usize) {
        let o = &self.job.options;
        let Some(font) = &self.font else { return };
        let (w, h) = (pic.width as f32, pic.height as f32);
        let short = w.min(h);
        let margin = (short * 0.035).max(6.0);
        let small = (short * 0.032).max(10.0);
        let gap = small * 0.9;
        let ink = self.ink();
        // How far in from its edge each corner is filled.
        let mut used = [0.0f32; 4];

        if !o.caption.trim().is_empty() {
            let size = (short * 0.05).max(12.0);
            text(
                pic,
                font,
                o.caption.trim(),
                margin,
                margin + size * 0.8,
                size,
                ink,
            );
            used[Corner::TopLeft.index()] += size + gap;
        }

        let mut info = Vec::new();
        let count = self.job.structure.frame_count();
        if o.counter && count > 1 {
            info.push(format!("frame {} / {}", shot.nearest() + 1, count));
        }
        if o.time {
            if let Some(step) = self.step(shot) {
                info.push(format!("step {}", thousands(step.round() as i64)));
                info.push(time_label(step * o.dt as f64, o.unit));
            }
        }
        if !info.is_empty() {
            text(
                pic,
                font,
                &info.join("   ·   "),
                margin,
                h - margin,
                small,
                ink,
            );
            used[Corner::BottomLeft.index()] += small + gap;
        }

        if o.scale_bar {
            // Pixels to the ångström at the depth of the centre, in the
            // picture as finally drawn.
            let per = cam_scale(cam, pic.width * ss, pic.height * ss) / ss as f32;
            if per > 0.0 {
                let (len, label) = nice_length(w * 0.18 / per);
                let bar = len * per;
                let thick = (short * 0.007).max(2.0);
                let (x1, y) = (w - margin, h - margin);
                rect(pic, x1 - bar, y - thick, bar, thick, ink);
                let tw = text_width(font, &label, small);
                text(
                    pic,
                    font,
                    &label,
                    x1 - bar / 2.0 - tw / 2.0,
                    y - thick - small * 0.45,
                    small,
                    ink,
                );
                used[Corner::BottomRight.index()] += thick + small * 1.3 + gap;
            }
        }

        if let (true, Some(legend)) = (o.legend, &self.job.legend) {
            let c = o.corner;
            let (lw, lh) = legend_size(font, legend, small, short);
            let x = if c.right() { w - margin - lw } else { margin };
            let y = if c.bottom() {
                h - margin - used[c.index()] - lh
            } else {
                margin + used[c.index()]
            };
            draw_legend(pic, font, legend, x, y, small, short, ink);
        }
    }

    /// The simulation step a shot is at, blended as its positions are.
    fn step(&self, shot: &Shot) -> Option<f64> {
        let frames = &self.job.structure.frames;
        let a = frames.get(shot.from)?.step? as f64;
        let b = frames
            .get(shot.to)
            .and_then(|f| f.step)
            .map(|s| s as f64)
            .unwrap_or(a);
        Some(a + (b - a) * shot.t as f64)
    }

    /// What colour writing stands out from the background in.
    fn ink(&self) -> [f32; 3] {
        let ground = match self.job.options.background {
            Background::Color(c) => c,
            Background::Theme | Background::Transparent => self.job.page,
        };
        let luma = 0.2126 * ground[0] + 0.7152 * ground[1] + 0.0722 * ground[2];
        if luma > 0.5 {
            [0.10, 0.10, 0.12]
        } else {
            [0.94, 0.94, 0.95]
        }
    }
}

/// Pixels to the ångström at the centre's depth, as the renderers frame a
/// picture `width` by `height`.
fn cam_scale(cam: &Camera, width: usize, height: usize) -> f32 {
    (width.min(height) as f32) / (2.2 * cam.radius.max(1e-6)) * cam.zoom
}

/// A round length for a scale bar no longer than `most` ångström, and how
/// it reads.
fn nice_length(most: f32) -> (f32, String) {
    let mut best = 1.0;
    let mut p = 0.01f32;
    while p < 1e7 {
        for m in [1.0, 2.0, 5.0] {
            if m * p <= most {
                best = m * p;
            }
        }
        p *= 10.0;
    }
    let label = if best >= 10.0 {
        format!("{} nm", trim(best / 10.0))
    } else {
        format!("{} Å", trim(best))
    };
    (best, label)
}

fn trim(x: f32) -> String {
    let s = format!("{x:.2}");
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

fn thousands(n: i64) -> String {
    let digits = n.unsigned_abs().to_string();
    let mut out = String::new();
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    if n < 0 {
        format!("-{out}")
    } else {
        out
    }
}

/// Simulated time, in the largest unit that keeps it above one.
fn time_label(t: f64, unit: TimeUnit) -> String {
    let fs = match unit {
        TimeUnit::Tau => return format!("t = {} τ", trim(t as f32)),
        TimeUnit::Fs => t,
        TimeUnit::Ps => t * 1e3,
    };
    let (v, u) = if fs.abs() >= 1e9 {
        (fs / 1e9, "µs")
    } else if fs.abs() >= 1e6 {
        (fs / 1e6, "ns")
    } else if fs.abs() >= 1e3 {
        (fs / 1e3, "ps")
    } else {
        (fs, "fs")
    };
    format!("t = {v:.2} {u}")
}

/// The legend's sizes, from the size of its writing: the title's, the
/// labels', a ramp's length and thickness, and a row of swatches' height.
fn legend_metrics(small: f32, short: f32) -> (f32, f32, f32, f32, f32) {
    (
        small * 0.95,
        small * 0.9,
        (short * 0.26).max(80.0),
        (small * 0.55).max(4.0),
        small * 1.45,
    )
}

fn legend_size(
    font: &ab_glyph::FontRef<'static>,
    legend: &Legend,
    small: f32,
    short: f32,
) -> (f32, f32) {
    let (title, label, length, thick, row) = legend_metrics(small, short);
    match legend {
        Legend::Ramp {
            title: t,
            low,
            high,
            ..
        } => {
            let ends = text_width(font, low, label) + text_width(font, high, label) + label;
            (
                length.max(ends).max(text_width(font, t, title)),
                title * 1.5 + thick + label * 1.5,
            )
        }
        Legend::Swatches {
            title: t,
            items,
            more,
        } => {
            let dot = label * 1.3;
            let widest = items
                .iter()
                .map(|(_, n)| text_width(font, n, label))
                .fold(0.0, f32::max);
            let rows = items.len() + usize::from(*more > 0);
            (
                text_width(font, t, title).max(dot + widest),
                title * 1.5 + rows as f32 * row,
            )
        }
    }
}

/// The legend with its top left corner at `x`, `y`.
#[allow(clippy::too_many_arguments)]
fn draw_legend(
    pic: &mut Picture,
    font: &ab_glyph::FontRef<'static>,
    legend: &Legend,
    x: f32,
    y: f32,
    small: f32,
    short: f32,
    ink: [f32; 3],
) {
    let (title_size, label, _, thick, row) = legend_metrics(small, short);
    let faint = ink.map(|v| v * 0.75 + 0.125);
    match legend {
        Legend::Ramp {
            title,
            colors,
            low,
            high,
        } => {
            text(pic, font, title, x, y + title_size * 0.85, title_size, ink);
            let (w, _) = legend_size(font, legend, small, short);
            let top = y + title_size * 1.5;
            let n = colors.len().max(1);
            // Each colour its share of the bar, blended into the next.
            for px in 0..w.ceil() as usize {
                let t = (px as f32 + 0.5) / w * (n - 1) as f32;
                let (i, f) = (t.floor() as usize, t.fract());
                let a = colors[i.min(n - 1)];
                let b = colors[(i + 1).min(n - 1)];
                let c = [
                    a[0] + (b[0] - a[0]) * f,
                    a[1] + (b[1] - a[1]) * f,
                    a[2] + (b[2] - a[2]) * f,
                ];
                rect(pic, x + px as f32, top, 1.0, thick, c);
            }
            let base = top + thick + label * 1.15;
            text(pic, font, low, x, base, label, faint);
            text(
                pic,
                font,
                high,
                x + w - text_width(font, high, label),
                base,
                label,
                faint,
            );
        }
        Legend::Swatches { title, items, more } => {
            text(pic, font, title, x, y + title_size * 0.85, title_size, ink);
            let mut top = y + title_size * 1.5;
            let r = label * 0.36;
            for (color, name) in items {
                let mid = top + row * 0.45;
                disc(pic, x + r, mid, r, *color);
                text(
                    pic,
                    font,
                    name,
                    x + label * 1.3,
                    mid + label * 0.35,
                    label,
                    ink,
                );
                top += row;
            }
            if *more > 0 {
                text(
                    pic,
                    font,
                    &format!("and {more} more"),
                    x + label * 1.3,
                    top + row * 0.45 + label * 0.35,
                    label,
                    faint,
                );
            }
        }
    }
}

/// A filled circle, its edge smoothed.
fn disc(pic: &mut Picture, cx: f32, cy: f32, r: f32, color: [f32; 3]) {
    let (x0, x1) = ((cx - r - 1.0).floor() as i64, (cx + r + 1.0).ceil() as i64);
    let (y0, y1) = ((cy - r - 1.0).floor() as i64, (cy + r + 1.0).ceil() as i64);
    for py in y0..y1 {
        for px in x0..x1 {
            let d = ((px as f32 + 0.5 - cx).powi(2) + (py as f32 + 0.5 - cy).powi(2)).sqrt();
            blend(pic, px, py, (r - d + 0.5).clamp(0.0, 1.0), color);
        }
    }
}

/// Scale an RGBA picture down by `k`, each pixel the average of the block
/// it covers. Colour is averaged weighted by how opaque it is, so a soft
/// edge on nothing does not darken.
fn shrink(big: &[u8], bw: usize, bh: usize, k: usize) -> Vec<u8> {
    if k == 1 {
        return big.to_vec();
    }
    let (w, h) = (bw / k, bh / k);
    let mut out = vec![0u8; w * h * 4];
    use rayon::prelude::*;
    out.par_chunks_mut(w * 4).enumerate().for_each(|(y, row)| {
        for x in 0..w {
            let mut acc = [0u32; 4];
            for dy in 0..k {
                let line = ((y * k + dy) * bw + x * k) * 4;
                for dx in 0..k {
                    let p = &big[line + dx * 4..line + dx * 4 + 4];
                    let a = p[3] as u32;
                    acc[0] += p[0] as u32 * a;
                    acc[1] += p[1] as u32 * a;
                    acc[2] += p[2] as u32 * a;
                    acc[3] += a;
                }
            }
            let o = &mut row[x * 4..x * 4 + 4];
            if acc[3] > 0 {
                for c in 0..3 {
                    o[c] = ((acc[c] + acc[3] / 2) / acc[3]).min(255) as u8;
                }
            }
            o[3] = ((acc[3] + (k * k) as u32 / 2) / (k * k) as u32) as u8;
        }
    });
    out
}

/// Lay `color` over a pixel with coverage `a`.
fn blend(pic: &mut Picture, x: i64, y: i64, a: f32, color: [f32; 3]) {
    if x < 0 || y < 0 || x >= pic.width as i64 || y >= pic.height as i64 || a <= 0.0 {
        return;
    }
    let o = (y as usize * pic.width + x as usize) * 4;
    let p = &mut pic.rgba[o..o + 4];
    let a = a.min(1.0);
    let da = p[3] as f32 / 255.0;
    let oa = a + da * (1.0 - a);
    for c in 0..3 {
        let d = p[c] as f32 / 255.0;
        let v = if oa > 0.0 {
            (color[c] * a + d * da * (1.0 - a)) / oa
        } else {
            0.0
        };
        p[c] = (v * 255.0 + 0.5) as u8;
    }
    p[3] = (oa * 255.0 + 0.5) as u8;
}

fn rect(pic: &mut Picture, x: f32, y: f32, w: f32, h: f32, color: [f32; 3]) {
    let (x0, x1, y0, y1) = (
        x.floor() as i64,
        (x + w).ceil() as i64,
        y.floor() as i64,
        (y + h).ceil() as i64,
    );
    for py in y0..y1 {
        let cy = ((py as f32 + 1.0).min(y + h) - (py as f32).max(y)).clamp(0.0, 1.0);
        for px in x0..x1 {
            let cx = ((px as f32 + 1.0).min(x + w) - (px as f32).max(x)).clamp(0.0, 1.0);
            blend(pic, px, py, cx * cy, color);
        }
    }
}

fn text_width(font: &ab_glyph::FontRef<'static>, s: &str, size: f32) -> f32 {
    use ab_glyph::{Font, ScaleFont};
    let f = font.as_scaled(size);
    let mut w = 0.0;
    let mut prev = None;
    for ch in s.chars() {
        let g = f.glyph_id(ch);
        if let Some(p) = prev {
            w += f.kern(p, g);
        }
        w += f.h_advance(g);
        prev = Some(g);
    }
    w
}

/// Write `s` with its baseline at `y`, from `x`.
fn text(
    pic: &mut Picture,
    font: &ab_glyph::FontRef<'static>,
    s: &str,
    x: f32,
    y: f32,
    size: f32,
    color: [f32; 3],
) {
    use ab_glyph::{point, Font, ScaleFont};
    let f = font.as_scaled(size);
    let mut caret = x;
    let mut prev = None;
    for ch in s.chars() {
        let id = f.glyph_id(ch);
        if let Some(p) = prev {
            caret += f.kern(p, id);
        }
        let glyph = id.with_scale_and_position(size, point(caret, y));
        caret += f.h_advance(id);
        prev = Some(id);
        if let Some(outlined) = font.outline_glyph(glyph) {
            let b = outlined.px_bounds();
            outlined.draw(|gx, gy, c| {
                blend(
                    pic,
                    b.min.x as i64 + gx as i64,
                    b.min.y as i64 + gy as i64,
                    c,
                    color,
                );
            });
        }
    }
}

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

/// Text pinned to points of the structure, as `cam` sees them in `pic`: a
/// little to the right of each point, on a patch of the page so it reads
/// over whatever is behind.
pub fn draw_labels(
    pic: &mut Picture,
    font: &ab_glyph::FontRef<'static>,
    cam: &Camera,
    labels: &[([f32; 3], &str)],
    size: f32,
    ink: [f32; 3],
    page: [f32; 3],
) {
    let (w, h) = (pic.width as f32, pic.height as f32);
    let m = cam.basis();
    let scale = w.min(h) / (2.2 * cam.radius) * cam.zoom;
    for (p, label) in labels {
        let q = [
            p[0] - cam.center[0],
            p[1] - cam.center[1],
            p[2] - cam.center[2],
        ];
        let d = |r: [f32; 3]| r[0] * q[0] + r[1] * q[1] + r[2] * q[2];
        let (v, _) = cam.warp([d(m[0]), d(m[1]), d(m[2])]);
        let x = w * 0.5 + cam.pan[0] + v[0] * scale + size * 0.4;
        let y = h * 0.5 + cam.pan[1] - v[1] * scale;
        let tw = text_width(font, label, size);
        let pad = size * 0.25;
        for py in (y - size * 0.6 - pad) as i64..(y + size * 0.45 + pad) as i64 {
            for px in (x - pad) as i64..(x + tw + pad) as i64 {
                blend(pic, px, py, 0.85, page);
            }
        }
        text(pic, font, label, x, y + size * 0.35, size, ink);
    }
}

/// The font the film's text is set in, for anything else drawing text.
pub fn font() -> Option<ab_glyph::FontRef<'static>> {
    ab_glyph::FontRef::try_from_slice(FONT).ok()
}

/// Takes drawn frames into a file.
trait Sink {
    fn push(&mut self, pic: &Picture) -> Result<()>;
    fn finish(self: Box<Self>) -> Result<()>;
}

/// Where ffmpeg is, if it is anywhere: on the path, or where package
/// managers put it, since an app opened from the Finder or the Start menu
/// is not given the shell's path.
pub fn ffmpeg() -> Option<PathBuf> {
    static FOUND: std::sync::OnceLock<Option<PathBuf>> = std::sync::OnceLock::new();
    FOUND
        .get_or_init(|| {
            which::which("ffmpeg").ok().or_else(|| {
                [
                    "/opt/homebrew/bin/ffmpeg",
                    "/usr/local/bin/ffmpeg",
                    "/usr/bin/ffmpeg",
                    "/snap/bin/ffmpeg",
                ]
                .into_iter()
                .map(PathBuf::from)
                .find(|p| p.is_file())
            })
        })
        .clone()
}

/// The H.264 encoder ffmpeg has, best first.
fn h264_encoder() -> Option<&'static str> {
    static FOUND: std::sync::OnceLock<Option<&'static str>> = std::sync::OnceLock::new();
    *FOUND.get_or_init(|| {
        let out = std::process::Command::new(ffmpeg()?)
            .args(["-hide_banner", "-encoders"])
            .output()
            .ok()?;
        let list = String::from_utf8_lossy(&out.stdout).to_string();
        [
            "libx264",
            "h264_videotoolbox",
            "h264_mf",
            "h264_nvenc",
            "h264_qsv",
            "h264_vaapi",
        ]
        .into_iter()
        .find(|e| list.split_whitespace().any(|w| w == *e))
    })
}

fn sink(o: &Options, path: &Path, width: usize, height: usize) -> Result<Box<dyn Sink>> {
    Ok(match o.format {
        Format::Gif => Box::new(Gif::new(path, o.fps)?),
        Format::Png => Box::new(Pngs::new(path)?),
        Format::Mp4 | Format::WebM | Format::Mov => Box::new(Ffmpeg::new(o, path, width, height)?),
    })
}

struct Gif {
    encoder: image::codecs::gif::GifEncoder<std::io::BufWriter<std::fs::File>>,
    delay: u32,
}

impl Gif {
    fn new(path: &Path, fps: u32) -> Result<Self> {
        let file = std::io::BufWriter::new(
            std::fs::File::create(path)
                .with_context(|| format!("could not write {}", path.display()))?,
        );
        let mut encoder = image::codecs::gif::GifEncoder::new_with_speed(file, 10);
        encoder.set_repeat(image::codecs::gif::Repeat::Infinite)?;
        // GIF counts in hundredths of a second, and most viewers draw
        // anything quicker than two of them far slower instead.
        let delay = (1000 / fps.max(1)).max(20);
        Ok(Gif { encoder, delay })
    }
}

impl Sink for Gif {
    fn push(&mut self, pic: &Picture) -> Result<()> {
        // A GIF pixel is there or not.
        let mut rgba = pic.rgba.clone();
        for p in rgba.chunks_exact_mut(4) {
            p[3] = if p[3] >= 128 { 255 } else { 0 };
        }
        let buffer = image::RgbaImage::from_raw(pic.width as u32, pic.height as u32, rgba)
            .expect("sized to fit");
        self.encoder.encode_frame(image::Frame::from_parts(
            buffer,
            0,
            0,
            image::Delay::from_numer_denom_ms(self.delay, 1),
        ))?;
        Ok(())
    }

    fn finish(self: Box<Self>) -> Result<()> {
        // The trailer is written as the encoder goes.
        drop(self.encoder);
        Ok(())
    }
}

struct Pngs {
    dir: PathBuf,
    stem: String,
    next: usize,
}

impl Pngs {
    /// Frames go in a folder made at `path`, each named after it.
    fn new(path: &Path) -> Result<Self> {
        let dir = path.with_extension("");
        std::fs::create_dir_all(&dir)
            .with_context(|| format!("could not make {}", dir.display()))?;
        let stem = dir
            .file_name()
            .map(|n| n.to_string_lossy().to_string())
            .unwrap_or_else(|| "frame".into());
        Ok(Pngs { dir, stem, next: 1 })
    }
}

impl Sink for Pngs {
    fn push(&mut self, pic: &Picture) -> Result<()> {
        let path = self.dir.join(format!("{}_{:05}.png", self.stem, self.next));
        image::save_buffer(
            &path,
            &pic.rgba,
            pic.width as u32,
            pic.height as u32,
            image::ColorType::Rgba8,
        )?;
        self.next += 1;
        Ok(())
    }

    fn finish(self: Box<Self>) -> Result<()> {
        Ok(())
    }
}

struct Ffmpeg {
    child: Option<std::process::Child>,
    stdin: Option<std::process::ChildStdin>,
}

impl Ffmpeg {
    fn new(o: &Options, path: &Path, width: usize, height: usize) -> Result<Self> {
        let exe = ffmpeg()
            .ok_or_else(|| anyhow!("{} needs ffmpeg, which was not found", o.format.label()))?;
        let mut cmd = std::process::Command::new(exe);
        cmd.args([
            "-hide_banner",
            "-loglevel",
            "error",
            "-y",
            "-f",
            "rawvideo",
            "-pix_fmt",
            "rgba",
        ])
        .args([
            "-s",
            &format!("{width}x{height}"),
            "-r",
            &o.fps.to_string(),
            "-i",
            "-",
        ]);
        match o.format {
            Format::Mp4 => {
                // x264 where ffmpeg has it, as most builds do; else the
                // system's own encoder.
                let q = o.quality as f32;
                match h264_encoder() {
                    Some("libx264") => {
                        let crf = (36.0 - q * 0.2).round().clamp(12.0, 40.0) as u32;
                        cmd.args([
                            "-c:v",
                            "libx264",
                            "-preset",
                            "slow",
                            "-crf",
                            &crf.to_string(),
                        ]);
                    }
                    Some("h264_videotoolbox") => {
                        cmd.args([
                            "-c:v",
                            "h264_videotoolbox",
                            "-q:v",
                            &((q * 0.9 + 10.0) as u32).to_string(),
                        ]);
                    }
                    Some(other) => {
                        let bps = (0.02 + 0.22 * (q / 100.0).powi(2))
                            * (width * height) as f32
                            * o.fps as f32;
                        cmd.args(["-c:v", other, "-b:v", &format!("{}", bps as u64)]);
                    }
                    None => {
                        bail!("this ffmpeg has no H.264 encoder; choose WebM, GIF or PNG frames")
                    }
                }
                cmd.args(["-pix_fmt", "yuv420p", "-movflags", "+faststart"]);
            }
            Format::WebM => {
                let crf = (52.0 - o.quality as f32 * 0.37).round().clamp(4.0, 63.0) as u32;
                cmd.args([
                    "-c:v",
                    "libvpx-vp9",
                    "-pix_fmt",
                    "yuva420p",
                    "-b:v",
                    "0",
                    "-crf",
                    &crf.to_string(),
                ])
                .args(["-row-mt", "1", "-auto-alt-ref", "0"]);
            }
            _ => {
                cmd.args([
                    "-c:v",
                    "prores_ks",
                    "-profile:v",
                    "4444",
                    "-pix_fmt",
                    "yuva444p10le",
                    "-vendor",
                    "apl0",
                ]);
            }
        }
        cmd.arg(path)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::piped());
        #[cfg(windows)]
        {
            use std::os::windows::process::CommandExt;
            const CREATE_NO_WINDOW: u32 = 0x0800_0000;
            cmd.creation_flags(CREATE_NO_WINDOW);
        }
        let mut child = cmd.spawn().context("could not start ffmpeg")?;
        let stdin = child.stdin.take();
        Ok(Ffmpeg {
            child: Some(child),
            stdin,
        })
    }
}

impl Sink for Ffmpeg {
    fn push(&mut self, pic: &Picture) -> Result<()> {
        use std::io::Write;
        let Some(stdin) = self.stdin.as_mut() else {
            bail!("ffmpeg has stopped")
        };
        if stdin.write_all(&pic.rgba).is_err() {
            bail!("ffmpeg stopped taking frames");
        }
        Ok(())
    }

    fn finish(mut self: Box<Self>) -> Result<()> {
        drop(self.stdin.take());
        let out = self.child.take().expect("taken once").wait_with_output()?;
        if !out.status.success() {
            let err = String::from_utf8_lossy(&out.stderr);
            bail!(
                "ffmpeg failed: {}",
                err.lines().last().unwrap_or("no reason given")
            );
        }
        Ok(())
    }
}

impl Drop for Ffmpeg {
    fn drop(&mut self) {
        // Cancelled partway: stop it rather than leave it waiting.
        if let Some(child) = self.child.as_mut() {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}

/// How far an export has got, for the progress bar, and the flag that
/// stops it.
#[derive(Default)]
pub struct Progress {
    pub done: AtomicUsize,
    pub total: AtomicUsize,
    pub cancel: AtomicBool,
}

/// Draw every shot and write the film to `path`. Stops early, and removes
/// what was written, when cancelled.
pub fn export(job: &Job, path: &Path, progress: &Progress) -> Result<PathBuf> {
    let o = &job.options;
    if o.format.needs_ffmpeg() && ffmpeg().is_none() {
        bail!(
            "{} needs ffmpeg; install it, or choose MP4, GIF or PNG frames",
            o.format.label()
        );
    }
    let path = match o.format {
        Format::Png => path.with_extension(""),
        f => path.with_extension(f.extension()),
    };
    let (width, height) = o.output_size(job.view);
    progress.total.store(job.shots.len(), Ordering::Relaxed);
    let mut out = sink(o, &path, width, height)?;
    let mut drawer = Drawer::new(job);
    for i in 0..job.shots.len() {
        if progress.cancel.load(Ordering::Relaxed) {
            drop(out);
            if o.format == Format::Png {
                let _ = std::fs::remove_dir_all(&path);
            } else {
                let _ = std::fs::remove_file(&path);
            }
            bail!("cancelled");
        }
        let pic = drawer.draw(i, width, height, o.supersample as usize);
        out.push(&pic)?;
        progress.done.store(i + 1, Ordering::Relaxed);
    }
    out.finish()?;
    Ok(path)
}

/// A few small frames of the film, for seeing it before it is made: no
/// more than `most`, taken evenly through it, `across` pixels wide at most.
/// With each, how many of the film's frames it stands for.
pub fn preview(
    job: &Job,
    across: usize,
    most: usize,
    cancel: &AtomicBool,
) -> (Vec<Picture>, usize) {
    let (w, h) = job.options.output_size(job.view);
    let k = (across as f32 / w.max(1) as f32).min(1.0);
    let (pw, ph) = (
        ((w as f32 * k) as usize).max(2),
        ((h as f32 * k) as usize).max(2),
    );
    let stride = job.shots.len().div_ceil(most.max(1)).max(1);
    let mut drawer = Drawer::new(job);
    let mut out = Vec::new();
    for i in (0..job.shots.len()).step_by(stride) {
        if cancel.load(Ordering::Relaxed) {
            break;
        }
        out.push(drawer.draw(i, pw, ph, 1));
    }
    (out, stride)
}

// ---------------------------------------------------------------------------
// The console's way of saying it
// ---------------------------------------------------------------------------

/// Options from `key=value` words, over `base`. The format comes from the
/// file's extension. Quoted values may have spaces in them.
pub fn parse_args(base: Options, path: &Path, words: &[String]) -> Result<Options> {
    let mut o = base;
    let ext = path
        .extension()
        .map(|e| e.to_string_lossy().to_string())
        .unwrap_or_default();
    o.format = Format::from_extension(&ext).ok_or_else(|| {
        anyhow!("a film is .mp4, .gif, .webm, .mov, or a folder of PNGs, not .{ext}")
    })?;
    let on = |v: &str| -> Result<bool> {
        match v {
            "on" | "yes" | "true" | "1" | "" => Ok(true),
            "off" | "no" | "false" | "0" => Ok(false),
            _ => bail!("\"{v}\" is neither on nor off"),
        }
    };
    let num = |k: &str, v: &str| -> Result<f32> {
        v.parse::<f32>()
            .map_err(|_| anyhow!("{k} takes a number, not \"{v}\""))
    };
    for word in words {
        let (k, v) = word.split_once('=').unwrap_or((word.as_str(), ""));
        let v = v.trim_matches('"');
        match k.to_ascii_lowercase().as_str() {
            "frames" => {
                // first:last:every, any part left out.
                let parts: Vec<&str> = v.split(':').collect();
                let at = |i: usize| {
                    parts
                        .get(i)
                        .filter(|s| !s.is_empty())
                        .map(|s| s.parse::<usize>())
                };
                o.trajectory = true;
                if let Some(x) = at(0) {
                    o.first = x.context("frames is first:last:every")?;
                }
                if let Some(x) = at(1) {
                    o.last = x.context("frames is first:last:every")?;
                }
                if let Some(x) = at(2) {
                    o.every = x.context("frames is first:last:every")?;
                }
            }
            "trajectory" => o.trajectory = on(v)?,
            "turn" | "rock" => {
                o.motion = if k == "turn" {
                    Motion::Turn
                } else {
                    Motion::Rock
                };
                if !v.is_empty() {
                    o.degrees = num(k, v)?;
                }
            }
            "still" => o.motion = Motion::Still,
            "recorded" => o.motion = Motion::Recorded,
            "axis" => {
                o.axis = named(&AXES, v)
                    .or(match v {
                        "vertical" => Some(Axis::Vertical),
                        "horizontal" => Some(Axis::Horizontal),
                        "view" | "roll" => Some(Axis::View),
                        _ => None,
                    })
                    .ok_or_else(|| anyhow!("axis is x, y or z"))?
            }
            "loop" | "seamless" => o.seamless = on(v)?,
            "fps" => o.fps = num(k, v)? as u32,
            "seconds" | "length" => o.seconds = num(k, v)?,
            "smooth" => o.smooth = num(k, v)? as u32,
            "pingpong" => o.pingpong = on(v)?,
            "size" => {
                if let Some(s) = named(&SIZES, v) {
                    o.size = s;
                } else if let Some((w, h)) = v.split_once(['x', '×']) {
                    o.size = Size::Custom;
                    o.width = num(k, w)? as u32;
                    o.height = num(k, h)? as u32;
                } else {
                    bail!("size is window, 720p, 1080p, 4k, square, or WIDTHxHEIGHT");
                }
            }
            "aa" | "supersample" => o.supersample = num(k, v)? as u32,
            "trace" | "traced" | "raytrace" => o.traced = on(v)?,
            "quality" => o.quality = num(k, v)? as u32,
            "bg" | "background" => {
                o.background =
                    match v {
                        "theme" | "page" => Background::Theme,
                        "transparent" | "none" => Background::Transparent,
                        other => Background::Color(crate::settings::parse_hex(other).ok_or_else(
                            || anyhow!("bg is theme, transparent or a #rrggbb colour"),
                        )?),
                    }
            }
            "counter" => o.counter = on(v)?,
            "time" => o.time = on(v)?,
            "dt" => {
                o.time = true;
                let (n, unit) = match v.find(|c: char| c.is_ascii_alphabetic() || c == 'τ') {
                    Some(i) => (&v[..i], &v[i..]),
                    None => (v, ""),
                };
                o.dt = num(k, n)?;
                if !unit.is_empty() {
                    o.unit = named(&UNITS, unit)
                        .or((unit == "τ").then_some(TimeUnit::Tau))
                        .ok_or_else(|| anyhow!("dt's unit is fs, ps or tau"))?;
                }
            }
            "scalebar" | "scale" => o.scale_bar = on(v)?,
            "legend" => {
                o.legend = true;
                if let Some(c) = named(&CORNERS, v) {
                    o.corner = c;
                } else {
                    o.legend = on(v).map_err(|_| {
                        anyhow!(
                            "legend is on, off, top-left, top-right, bottom-left or bottom-right"
                        )
                    })?;
                }
            }
            "caption" | "title" => o.caption = v.to_string(),
            _ => bail!("{k} is not something a film takes"),
        }
    }
    Ok(o.clamped())
}

/// Split a line into words, keeping what is in quotes together.
pub fn words(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut quoted = false;
    for c in line.chars() {
        match c {
            '"' => {
                quoted = !quoted;
                cur.push(c);
            }
            c if c.is_whitespace() && !quoted => {
                if !cur.is_empty() {
                    out.push(std::mem::take(&mut cur));
                }
            }
            c => cur.push(c),
        }
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protein::Snapshot;

    fn camera() -> Camera {
        Camera {
            orientation: Camera::home(),
            zoom: 1.0,
            pan: [0.0; 2],
            center: [0.0; 3],
            radius: 10.0,
            perspective: false,
            slab: None,
        }
    }

    #[test]
    fn a_trajectory_is_one_shot_a_frame() {
        let o = Options {
            first: 3,
            last: 9,
            every: 2,
            ..Options::default()
        };
        let shots = plan(&o, &camera(), 0, 20, None);
        let frames: Vec<usize> = shots.iter().map(|s| s.from).collect();
        assert_eq!(frames, vec![2, 4, 6, 8]);
    }

    #[test]
    fn smoothing_blends_between_the_frames_used() {
        let o = Options {
            every: 2,
            last: 5,
            smooth: 1,
            ..Options::default()
        };
        let shots = plan(&o, &camera(), 0, 10, None);
        // Frames 0, 2, 4 with one blend between each pair.
        assert_eq!(shots.len(), 5);
        assert_eq!((shots[1].from, shots[1].to), (0, 2));
        assert!((shots[1].t - 0.5).abs() < 1e-4);
        assert_eq!(shots[4].nearest(), 4);
    }

    #[test]
    fn ping_pong_comes_back_without_doubling_the_ends() {
        let o = Options {
            last: 4,
            pingpong: true,
            ..Options::default()
        };
        let frames: Vec<usize> = plan(&o, &camera(), 0, 4, None)
            .iter()
            .map(|s| s.nearest())
            .collect();
        assert_eq!(frames, vec![0, 1, 2, 3, 2, 1]);
    }

    #[test]
    fn a_seamless_turn_stops_short_of_where_it_began() {
        let o = Options {
            trajectory: false,
            motion: Motion::Turn,
            seconds: 1.0,
            fps: 4,
            ..Options::default()
        };
        let base = camera();
        let shots = plan(&o, &base, 0, 1, None);
        assert_eq!(shots.len(), 4);
        // The next shot after the last would be the first again.
        let next = turned(&base, Axis::Vertical, std::f32::consts::TAU);
        for r in 0..3 {
            for k in 0..3 {
                assert!((next.orientation[r][k] - base.orientation[r][k]).abs() < 1e-4);
            }
        }
        assert!((shots[1].camera.orientation[0][0] - base.orientation[0][0]).abs() > 1e-3);
    }

    #[test]
    fn a_recording_plays_back_at_the_films_rate() {
        let mut far = camera();
        far.zoom = 3.0;
        let rec = [
            Sample {
                t: 0.0,
                camera: camera(),
                frame: 0,
            },
            Sample {
                t: 1.0,
                camera: far,
                frame: 4,
            },
        ];
        let o = Options {
            motion: Motion::Recorded,
            fps: 10,
            smooth: 1,
            ..Options::default()
        };
        let shots = plan(&o, &camera(), 0, 5, Some(&rec));
        assert_eq!(shots.len(), 11);
        assert!((shots[5].camera.zoom - 2.0).abs() < 1e-3);
        assert_eq!((shots[5].from, shots[5].to), (0, 4));
    }

    #[test]
    fn blending_takes_the_short_way_across_the_cell() {
        let cell = Cell::from_bounds([0.0; 3], [10.0; 3]);
        let frames = vec![
            Snapshot {
                step: Some(0),
                pos: vec![[9.5, 5.0, 5.0]],
                cell: Some(cell),
            },
            Snapshot {
                step: Some(100),
                pos: vec![[0.5, 5.0, 5.0]],
                cell: Some(cell),
            },
        ];
        let mut s = Structure::beads("box", 1, frames);
        s.set_frame(0);
        pose(
            &mut s,
            &Shot {
                camera: camera(),
                from: 0,
                to: 1,
                t: 0.5,
            },
        );
        // Through the wall at 10, not back across the box through 5.
        assert!(
            (s.atoms[0].pos[0] - 10.0).abs() < 1e-4,
            "{:?}",
            s.atoms[0].pos
        );
    }

    #[test]
    fn shrinking_averages_by_coverage() {
        // One opaque red pixel and three empty ones.
        let big = [[255, 0, 0, 255], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]].concat();
        let small = shrink(&big, 2, 2, 2);
        assert_eq!(&small[..3], &[255, 0, 0]);
        assert_eq!(small[3], 64);
    }

    #[test]
    fn scale_bars_are_round_lengths() {
        assert_eq!(nice_length(37.0).1, "2 nm");
        assert_eq!(nice_length(6.0).1, "5 Å");
        assert_eq!(nice_length(0.3).1, "0.2 Å");
    }

    #[test]
    fn times_read_in_the_largest_unit() {
        assert_eq!(time_label(2_500_000.0, TimeUnit::Fs), "t = 2.50 ns");
        assert_eq!(time_label(1.5, TimeUnit::Ps), "t = 1.50 ps");
        assert_eq!(time_label(12.5, TimeUnit::Tau), "t = 12.5 τ");
    }

    #[test]
    fn the_console_says_it_in_words() {
        let w = words(
            r#"frames=1:500:2 fps=24 turn=360 bg=transparent size=1280x720 caption="Melting LJ crystal" dt=0.005tau legend=bottom-left"#,
        );
        let o = parse_args(Options::default(), Path::new("out.webm"), &w).unwrap();
        assert_eq!((o.first, o.last, o.every), (1, 500, 2));
        assert_eq!(o.fps, 24);
        assert_eq!(o.motion, Motion::Turn);
        assert_eq!(o.background, Background::Transparent);
        assert_eq!((o.size, o.width, o.height), (Size::Custom, 1280, 720));
        assert_eq!(o.caption, "Melting LJ crystal");
        assert_eq!(o.unit, TimeUnit::Tau);
        assert_eq!(o.format, Format::WebM);
        assert!(o.legend);
        assert_eq!(o.corner, Corner::BottomLeft);
        assert!(parse_args(Options::default(), Path::new("out.avi"), &[]).is_err());
    }

    #[test]
    fn options_survive_a_round_trip() {
        let o = Options {
            motion: Motion::Rock,
            degrees: 30.0,
            background: Background::Color([1.0, 0.0, 0.0]),
            caption: "hi".into(),
            format: Format::Gif,
            ..Options::default()
        };
        let back = Options::from_json(&serde_json::from_str(&o.to_json()).unwrap());
        assert_eq!(back, o);
    }

    #[test]
    fn an_mp4_is_written_through_ffmpeg() {
        let Some(exe) = ffmpeg() else { return };
        let dir = std::env::temp_dir().join(format!("profold-movie-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.mp4");
        let o = Options {
            fps: 10,
            ..Options::default()
        };
        let mut sink = Ffmpeg::new(&o, &path, 64, 48).unwrap();
        for i in 0..12u8 {
            let rgba: Vec<u8> = (0..64 * 48)
                .flat_map(|p| [i * 20, (p % 64) as u8 * 4, 128, 255])
                .collect();
            sink.push(&Picture {
                width: 64,
                height: 48,
                rgba,
            })
            .unwrap();
        }
        Box::new(sink).finish().unwrap();
        // Decoded back, it has every frame.
        let out = std::process::Command::new(exe)
            .args(["-hide_banner", "-i"])
            .arg(&path)
            .args(["-f", "null", "-"])
            .output()
            .unwrap();
        let log = String::from_utf8_lossy(&out.stderr);
        assert!(log.contains("h264"), "{log}");
        assert!(
            log.contains("frame=   12") || log.contains("frame=12"),
            "{log}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_traced_film_frame_is_drawn_at_its_size_and_differs_from_a_plain_one() {
        let melt = crate::samples::SAMPLES
            .iter()
            .find(|s| s.file.ends_with("lammpstrj.gz"))
            .unwrap();
        let s = crate::formats::load_bytes(melt.file, melt.bytes.to_vec()).unwrap();
        let camera = Camera::frame(&s);
        let frame = |traced: bool| {
            let options = Options {
                traced,
                trajectory: false,
                ..Options::default()
            };
            let shots = plan(&options, &camera, 0, s.frame_count(), None);
            let job = Job {
                structure: s.clone(),
                build: Arc::new(|s: &Structure| {
                    render::build_scene(
                        s,
                        render::Repr::Spacefill,
                        render::Coloring::Element,
                        1.0,
                        Style::Studio,
                        render::Palette::Spectrum,
                        true,
                        true,
                        &render::Appearance::default(),
                    )
                }),
                style: Style::Studio,
                page: [0.97, 0.97, 0.96],
                view: (320, 240),
                shots,
                options,
                legend: None,
                labels: Vec::new(),
            };
            Drawer::new(&job).draw(0, 160, 120, 2)
        };
        let (plain, traced) = (frame(false), frame(true));
        assert_eq!(
            (traced.width, traced.height, traced.rgba.len()),
            (160, 120, 160 * 120 * 4)
        );
        assert_ne!(plain.rgba, traced.rgba);
        let o = parse_args(Options::default(), Path::new("out.mp4"), &words("trace=on")).unwrap();
        assert!(o.traced);
    }

    #[test]
    fn a_label_is_drawn_where_its_atom_is_in_each_frame() {
        let melt = crate::samples::SAMPLES
            .iter()
            .find(|s| s.file.ends_with("lammpstrj.gz"))
            .unwrap();
        let s = crate::formats::load_bytes(melt.file, melt.bytes.to_vec()).unwrap();
        let camera = Camera::frame(&s);
        let frame = |labels: Vec<(u32, String)>| {
            let options = Options {
                trajectory: false,
                ..Options::default()
            };
            let shots = plan(&options, &camera, 0, s.frame_count(), None);
            let job = Job {
                structure: s.clone(),
                build: Arc::new(|s: &Structure| {
                    render::build_scene(
                        s,
                        render::Repr::Spacefill,
                        render::Coloring::Element,
                        1.0,
                        Style::Flat,
                        render::Palette::Spectrum,
                        true,
                        true,
                        &render::Appearance::default(),
                    )
                }),
                style: Style::Flat,
                page: [0.97, 0.97, 0.96],
                view: (320, 240),
                shots,
                options,
                legend: None,
                labels,
            };
            Drawer::new(&job).draw(0, 320, 240, 1)
        };
        let (plain, labelled) = (frame(Vec::new()), frame(vec![(0, "first atom".into())]));
        let changed = plain
            .rgba
            .chunks_exact(4)
            .zip(labelled.rgba.chunks_exact(4))
            .filter(|(a, b)| a != b)
            .count();
        assert!(changed > 100, "{changed} pixels changed");
    }

    #[test]
    fn a_legend_sits_in_its_corner() {
        let font = ab_glyph::FontRef::try_from_slice(FONT).unwrap();
        let legend = Legend::Ramp {
            title: "Hydrophobicity".into(),
            colors: vec![[0.2, 0.5, 0.8], [1.0, 1.0, 1.0], [0.9, 0.5, 0.1]],
            low: "hydrophilic".into(),
            high: "hydrophobic".into(),
        };
        let mut pic = Picture {
            width: 400,
            height: 300,
            rgba: vec![0; 400 * 300 * 4],
        };
        let (lw, lh) = legend_size(&font, &legend, 12.0, 300.0);
        draw_legend(
            &mut pic,
            &font,
            &legend,
            400.0 - 10.0 - lw,
            300.0 - 10.0 - lh,
            12.0,
            300.0,
            [0.0; 3],
        );
        let alpha = |x: usize, y: usize| pic.rgba[(y * 400 + x) * 4 + 3];
        let ink_in = |x0: usize, x1: usize, y0: usize, y1: usize| {
            (y0..y1).any(|y| (x0..x1).any(|x| alpha(x, y) > 0))
        };
        assert!(ink_in(200, 400, 150, 300), "nothing in the bottom right");
        assert!(
            !ink_in(0, 200, 0, 300) && !ink_in(200, 400, 0, 150),
            "drawn outside its corner"
        );
    }
}

/// Films of the melting crystal in every format, for looking at:
/// `PROFOLD_FILMS=<dir> cargo test --release film_samples -- --ignored`.
#[cfg(test)]
#[test]
#[ignore]
fn film_samples() {
    let Ok(dir) = std::env::var("PROFOLD_FILMS") else {
        return;
    };
    let dir = PathBuf::from(dir);
    std::fs::create_dir_all(&dir).unwrap();
    let melt = crate::samples::SAMPLES
        .iter()
        .find(|s| s.file.ends_with("lammpstrj.gz"))
        .unwrap();
    let s = crate::formats::load_bytes(melt.file, melt.bytes.to_vec()).unwrap();
    let camera = Camera::frame(&s);
    for (name, format, background) in [
        ("melt.mp4", Format::Mp4, Background::Theme),
        (
            "melt.gif",
            Format::Gif,
            Background::Color([0.08, 0.09, 0.12]),
        ),
        ("melt.webm", Format::WebM, Background::Transparent),
        ("melt.mov", Format::Mov, Background::Transparent),
        ("melt-frames", Format::Png, Background::Transparent),
    ] {
        let options = Options {
            motion: Motion::Rock,
            degrees: 25.0,
            smooth: 3,
            size: Size::Hd720,
            format,
            background,
            counter: true,
            time: true,
            dt: 0.005,
            unit: TimeUnit::Tau,
            scale_bar: true,
            legend: true,
            corner: Corner::TopRight,
            caption: "Lennard-Jones crystal melting".into(),
            ..Options::default()
        };
        let shots = plan(&options, &camera, 0, s.frame_count(), None);
        let job = Job {
            structure: s.clone(),
            build: Arc::new(|s: &Structure| {
                render::build_scene(
                    s,
                    render::Repr::Spacefill,
                    render::Coloring::Element,
                    1.0,
                    Style::Illustrative,
                    render::Palette::Spectrum,
                    true,
                    true,
                    &render::Appearance::default(),
                )
            }),
            style: Style::Illustrative,
            page: [0.97, 0.97, 0.96],
            view: (1280, 720),
            shots,
            options,
            legend: Some(Legend::Swatches {
                title: "Atoms".into(),
                items: vec![
                    ([0.13, 0.45, 0.85], "Type 1".into()),
                    ([0.90, 0.35, 0.15], "Type 2".into()),
                ],
                more: 0,
            }),
            labels: Vec::new(),
        };
        let t = std::time::Instant::now();
        let progress = Progress::default();
        let path = export(&job, &dir.join(name), &progress).unwrap();
        eprintln!(
            "{} · {} frames · {:.1}s",
            path.display(),
            job.shots.len(),
            t.elapsed().as_secs_f32()
        );
    }
}
