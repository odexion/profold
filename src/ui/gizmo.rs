//! The orientation axes in the bottom right corner of the viewport.
//!
//! The structure's own x, y and z, turning with it, drawn in ink like the
//! rest of the viewer. An axis that comes out of the screen towards the eye
//! is a solid line ending in a dot; one that goes back into it is a faint
//! dashed line ending in a ring. Clicking an axis's letter turns the view to
//! look straight down it at the structure, from its positive end, or from
//! its negative end with the command key (control off the Mac).
//!
//! Drawn over the viewport only: it is never part of a saved picture.

use std::rc::Rc;

use gpui::{
    canvas, div, px, App, InteractiveElement, IntoElement, MouseButton, ParentElement, Styled,
    Window,
};

use crate::theme;

/// The side of the square the axes are drawn in.
pub const SIZE: f32 = 56.0;

/// How far an axis lying flat in the screen reaches from the middle.
const REACH: f32 = 17.0;

/// One end of an axis as drawn: which axis, which end, where from the
/// middle in pixels (down the screen is +y), and how far it points out of
/// the screen, from -1 straight in to 1 straight out at the eye.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct End {
    pub axis: usize,
    pub positive: bool,
    pub at: [f32; 2],
    pub depth: f32,
}

/// Every end of the three axes under `orientation` (world to view, as
/// rows), farthest first, so that drawn in order the near ones lie on top.
pub fn ends(orientation: [[f32; 3]; 3]) -> Vec<End> {
    let mut ends: Vec<End> = (0..3)
        .flat_map(|axis| {
            // The axis in view space is the orientation's column.
            let v = [0, 1, 2].map(|r| orientation[r][axis]);
            [true, false].map(|positive| {
                let s = if positive { 1.0 } else { -1.0 };
                End {
                    axis,
                    positive,
                    at: [s * v[0] * REACH, -s * v[1] * REACH],
                    depth: s * v[2],
                }
            })
        })
        .collect();
    ends.sort_by(|a, b| a.depth.total_cmp(&b.depth));
    ends
}

/// The orientation that looks straight down an axis at the structure, that
/// end towards the eye. Y is kept up where it can be; looking down y, z is
/// up when +y faces the eye and down when -y does, so that x reads across to
/// the right in both, as a top and a bottom view do.
pub fn facing(axis: usize, positive: bool) -> [[f32; 3]; 3] {
    let s = if positive { 1.0 } else { -1.0 };
    let unit = |i: usize, k: f32| {
        let mut v = [0.0; 3];
        v[i] = k;
        v
    };
    let toward = unit(axis, s);
    let up = if axis == 1 { unit(2, -s) } else { unit(1, 1.0) };
    let across = [
        up[1] * toward[2] - up[2] * toward[1],
        up[2] * toward[0] - up[0] * toward[2],
        up[0] * toward[1] - up[1] * toward[0],
    ];
    [across, up, toward]
}

const NAMES: [&str; 3] = ["X", "Y", "Z"];

/// How solid an axis is drawn: faint and dashed pointing away, solid
/// pointing out, passing smoothly from one to the other as it crosses the
/// screen.
fn solidity(depth: f32) -> f32 {
    let t = ((depth + 0.15) / 0.3).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// The axes for a view turned by `orientation`. `pick` is told which axis
/// was clicked, and whether to look from its positive end: from the
/// negative one when the command key, or control, is held.
pub fn view(
    orientation: [[f32; 3]; 3],
    pick: impl Fn(&(usize, bool), &mut Window, &mut App) + 'static,
) -> impl IntoElement {
    let pick = Rc::new(pick);
    let middle = SIZE * 0.5;
    // Only the positive ends are drawn; how far one leans in or out is
    // what tells front from back.
    let axes: Vec<End> = ends(orientation)
        .into_iter()
        .filter(|e| e.positive)
        .collect();
    let ink = |depth: f32| theme::faint().blend(theme::ink().opacity(solidity(depth)));

    let lines = {
        let axes = axes.clone();
        canvas(
            |_, _, _| {},
            move |bounds, _, window, _| {
                let at = |x: f32, y: f32| bounds.origin + gpui::point(px(x), px(y));
                for e in &axes {
                    let color = ink(e.depth);
                    let front = solidity(e.depth) >= 0.5;
                    let (x, y) = (middle + e.at[0], middle + e.at[1]);
                    let mut path = gpui::PathBuilder::stroke(px(1.25));
                    if !front {
                        path = path.dash_array(&[px(2.), px(2.)]);
                    }
                    path.move_to(at(middle, middle));
                    path.line_to(at(x, y));
                    if let Ok(path) = path.build() {
                        window.paint_path(path, color);
                    }
                    let r = 1.9 + 0.7 * e.depth;
                    let dot =
                        gpui::Bounds::new(at(x - r, y - r), gpui::size(px(2. * r), px(2. * r)));
                    window.paint_quad(if front {
                        gpui::fill(dot, color).corner_radii(px(r))
                    } else {
                        gpui::outline(dot, color, gpui::BorderStyle::Solid).corner_radii(px(r))
                    });
                }
            },
        )
        .absolute()
        .inset_0()
    };

    let mut frame = div()
        .debug_selector(|| "axes".into())
        .relative()
        .size(px(SIZE))
        .child(lines);
    for e in axes {
        // The letter just past the end of its axis, or above and to the
        // right of the middle when the axis points straight along the view.
        let len = (e.at[0] * e.at[0] + e.at[1] * e.at[1]).sqrt();
        let (dx, dy) = if len > 1.0 {
            (e.at[0] / len, e.at[1] / len)
        } else {
            (0.7, -0.7)
        };
        let (x, y) = (middle + e.at[0] + 7.5 * dx, middle + e.at[1] + 7.5 * dy);
        let color = ink(e.depth);
        let pick = pick.clone();
        let axis = e.axis;
        frame = frame.child(
            div()
                .debug_selector(move || format!("axis-{}", NAMES[e.axis].to_lowercase()))
                .absolute()
                .left(px(x - 7.))
                .top(px(y - 7.))
                .size(px(14.))
                .flex()
                .items_center()
                .justify_center()
                .font_family(theme::MONO)
                .text_size(px(9.))
                .text_color(color)
                .cursor_pointer()
                .hover(|s| s.text_color(theme::ink()))
                .on_mouse_down(MouseButton::Left, move |event, window, cx| {
                    // The viewport underneath would take it as a drag.
                    cx.stop_propagation();
                    let m = &event.modifiers;
                    pick(&(axis, !(m.platform || m.control)), window, cx);
                })
                .child(NAMES[e.axis]),
        );
    }
    frame
}

#[cfg(test)]
mod tests {
    use super::*;

    const IDENTITY: [[f32; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    fn end(ends: &[End], axis: usize, positive: bool) -> End {
        *ends
            .iter()
            .find(|e| e.axis == axis && e.positive == positive)
            .unwrap()
    }

    #[test]
    fn square_on_x_runs_right_y_up_and_z_out_at_the_eye() {
        let e = ends(IDENTITY);
        assert_eq!(end(&e, 0, true).at, [REACH, 0.0]);
        assert_eq!(end(&e, 1, true).at, [0.0, -REACH]);
        let z = end(&e, 2, true);
        assert_eq!((z.at, z.depth), ([0.0, 0.0], 1.0));
        assert_eq!(end(&e, 2, false).depth, -1.0);
        // Farthest first: -z, then the four in the screen, then +z.
        assert_eq!((e[0].axis, e[0].positive), (2, false));
        assert_eq!((e[5].axis, e[5].positive), (2, true));
    }

    #[test]
    fn the_ends_follow_the_view() {
        // Turned a quarter about the vertical: x goes into the screen and z
        // comes round to the right.
        let turned = [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]];
        let e = ends(turned);
        assert_eq!(end(&e, 0, true).depth, -1.0);
        assert_eq!(end(&e, 2, true).at, [REACH, 0.0]);
        assert!(e.windows(2).all(|w| w[0].depth <= w[1].depth));
    }

    #[test]
    fn facing_an_axis_puts_that_end_at_the_eye() {
        for axis in 0..3 {
            for positive in [true, false] {
                let m = facing(axis, positive);
                // A proper rotation: orthonormal rows, right-handed.
                for i in 0..3 {
                    for j in 0..3 {
                        let d: f32 = (0..3).map(|k| m[i][k] * m[j][k]).sum();
                        assert_eq!(d, if i == j { 1.0 } else { 0.0 }, "{axis} {positive}");
                    }
                }
                let [a, b, c] = m;
                let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                    + a[2] * (b[0] * c[1] - b[1] * c[0]);
                assert_eq!(det, 1.0);
                let e = end(&ends(m), axis, positive);
                assert_eq!((e.at, e.depth), ([0.0, 0.0], 1.0));
            }
        }
        // Looking down y, x still reads across to the right.
        for positive in [true, false] {
            assert_eq!(end(&ends(facing(1, positive)), 0, true).at, [REACH, 0.0]);
        }
    }

    #[test]
    fn solid_out_hollow_in() {
        assert_eq!(solidity(1.0), 1.0);
        assert_eq!(solidity(-1.0), 0.0);
        assert_eq!(solidity(0.0), 0.5);
    }
}
