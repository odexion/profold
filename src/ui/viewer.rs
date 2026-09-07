//! The 3D viewport.
//!
//! A GPUI `canvas` that rasterizes the molecule on the CPU and uploads the
//! result as an image. The frame is cached against the camera and viewport
//! size, so idle redraws cost nothing.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use gpui::{canvas, px, Corners, IntoElement, RenderImage, Styled};
use image::{Frame, RgbaImage};
use smallvec::smallvec;

use crate::render::{self, Camera, Scene, Style};
use crate::theme;

#[derive(PartialEq, Eq)]
struct Key {
    width: usize,
    height: usize,
    /// Camera parameters, compared bitwise.
    cam: [u32; 5],
    scene: usize,
    len: usize,
    style: Style,
}

fn key(scene: &Arc<Scene>, cam: &Camera, width: usize, height: usize, style: Style) -> Key {
    Key {
        width,
        height,
        cam: [
            cam.yaw.to_bits(),
            cam.pitch.to_bits(),
            cam.zoom.to_bits(),
            cam.pan[0].to_bits(),
            cam.pan[1].to_bits(),
        ],
        scene: Arc::as_ptr(scene) as usize,
        len: scene.len(),
        style,
    }
}

#[derive(Default)]
pub struct Cache {
    key: Option<Key>,
    image: Option<Arc<RenderImage>>,
}

/// `interactive` renders at 1x instead of the display's scale factor, which
/// keeps dragging smooth on a retina screen.
pub fn viewport(
    scene: Arc<Scene>,
    cam: Camera,
    cache: Rc<RefCell<Cache>>,
    interactive: bool,
    style: Style,
) -> impl IntoElement {
    canvas(
        |_, _, _| {},
        move |bounds, _, window, _cx| {
            let scale = if interactive {
                1.0
            } else {
                window.scale_factor()
            };
            let width = (f32::from(bounds.size.width) * scale).round().max(1.0) as usize;
            let height = (f32::from(bounds.size.height) * scale).round().max(1.0) as usize;

            let mut cache = cache.borrow_mut();
            let k = key(&scene, &cam, width, height, style);
            if cache.key.as_ref() != Some(&k) {
                // A drawn figure sits on plain paper; the soft style keeps
                // its gradient.
                let (top, bottom) = match style {
                    Style::Illustrative => (theme::VIEW_FLAT, theme::VIEW_FLAT),
                    Style::Soft => (theme::VIEW_TOP, theme::VIEW_BOTTOM),
                };
                let frame = render::render(&scene, &cam, width, height, top, bottom, style);
                let buffer =
                    RgbaImage::from_raw(frame.width as u32, frame.height as u32, frame.bgra)
                    .expect("frame buffer is width * height * 4");
                let image = Arc::new(RenderImage::new(smallvec![Frame::new(buffer)]));
                if let Some(old) = cache.image.take() {
                    let _ = window.drop_image(old);
                }
                cache.image = Some(image);
                cache.key = Some(k);
            }

            if let Some(image) = cache.image.clone() {
                let _ = window.paint_image(bounds, Corners::all(px(0.)), image, 0, false);
            }
        },
    )
    .size_full()
}
