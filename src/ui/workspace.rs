//! The window: a row of tabs across the top, and the tab in front beneath.
//!
//! Each tab is a viewer of its own, with its own structure, camera and job,
//! so a protein, a trajectory and a molecule can be open side by side and a
//! fold can run in one tab while the others are used. A tab for folding is
//! opened only by asking for one, from the + menu.

use std::time::Duration;

use gpui::{
    actions, deferred, div, prelude::*, px, relative, Animation, AnimationExt, App, Context,
    Entity, ExternalPaths, FocusHandle, Focusable, MouseButton, PathPromptOptions, SharedString,
    Task, Window,
};

use crate::theme;
use crate::ui::app::{menu, swallow, Mode, Viewer};

actions!(
    profold,
    [
        OpenFile,
        NewTab,
        NewFoldTab,
        CloseTab,
        NextTab,
        PreviousTab,
        Tab1,
        Tab2,
        Tab3,
        Tab4,
        Tab5,
        Tab6,
        Tab7,
        Tab8,
        Tab9,
    ]
);

/// Height of the tab row, which is also the window's title bar.
const TAB_BAR: f32 = 40.0;
/// Room for the traffic lights at the left of the bar.
const LIGHTS: f32 = 80.0;

pub fn bind_keys(cx: &mut App) {
    cx.bind_keys([
        gpui::KeyBinding::new("cmd-o", OpenFile, None),
        gpui::KeyBinding::new("cmd-t", NewTab, None),
        gpui::KeyBinding::new("cmd-w", CloseTab, None),
        gpui::KeyBinding::new("ctrl-tab", NextTab, None),
        gpui::KeyBinding::new("ctrl-shift-tab", PreviousTab, None),
        gpui::KeyBinding::new("cmd-}", NextTab, None),
        gpui::KeyBinding::new("cmd-{", PreviousTab, None),
        gpui::KeyBinding::new("cmd-1", Tab1, None),
        gpui::KeyBinding::new("cmd-2", Tab2, None),
        gpui::KeyBinding::new("cmd-3", Tab3, None),
        gpui::KeyBinding::new("cmd-4", Tab4, None),
        gpui::KeyBinding::new("cmd-5", Tab5, None),
        gpui::KeyBinding::new("cmd-6", Tab6, None),
        gpui::KeyBinding::new("cmd-7", Tab7, None),
        gpui::KeyBinding::new("cmd-8", Tab8, None),
        gpui::KeyBinding::new("cmd-9", Tab9, None),
    ]);
}

pub struct Workspace {
    tabs: Vec<Entity<Viewer>>,
    active: usize,
    /// Whether the + menu is open.
    adding: bool,
    focus: FocusHandle,
    /// The open dialog, while it is up.
    prompt: Option<Task<()>>,
}

impl Workspace {
    pub fn new(window: &mut Window, cx: &mut Context<Self>) -> Self {
        let mut this = Workspace {
            tabs: Vec::new(),
            active: 0,
            adding: false,
            focus: cx.focus_handle(),
            prompt: None,
        };
        this.add(Mode::View, window, cx);
        this
    }

    /// Open a new tab beside the one in front, and bring it forward.
    pub fn add(&mut self, mode: Mode, window: &mut Window, cx: &mut Context<Self>) -> Entity<Viewer> {
        let tab = cx.new(|cx| Viewer::new(mode, window, cx));
        // The tab's name and activity show in the bar, which is drawn here.
        cx.observe(&tab, |_, _, cx| cx.notify()).detach();
        let at = if self.tabs.is_empty() { 0 } else { self.active + 1 };
        self.tabs.insert(at, tab.clone());
        self.active = at;
        self.adding = false;
        tab.update(cx, |tab, cx| tab.focus_input(window, cx));
        cx.notify();
        tab
    }

    /// The tab in front.
    pub fn active(&self) -> Entity<Viewer> {
        self.tabs[self.active].clone()
    }

    fn activate(&mut self, i: usize, window: &mut Window, cx: &mut Context<Self>) {
        if i >= self.tabs.len() {
            return;
        }
        self.active = i;
        self.adding = false;
        let tab = self.tabs[i].clone();
        // An empty tab is there to be typed into; a full one, to be looked at.
        tab.update(cx, |tab, cx| {
            if tab.is_blank() {
                tab.focus_input(window, cx);
            } else {
                window.focus(&tab.focus_handle(cx));
            }
        });
        cx.notify();
    }

    /// Close a tab. The last one is not closed but emptied, so the window
    /// always has a tab to type into.
    fn close(&mut self, i: usize, window: &mut Window, cx: &mut Context<Self>) {
        if i >= self.tabs.len() {
            return;
        }
        self.tabs.remove(i);
        if self.tabs.is_empty() {
            self.add(Mode::View, window, cx);
            return;
        }
        if self.active > i || self.active >= self.tabs.len() {
            self.active = self.active.saturating_sub(1);
        }
        self.activate(self.active, window, cx);
    }

    /// The tab something opened from the window goes into: the one in front
    /// when it is a viewer, and a new viewer when it is a fold tab.
    fn viewer_for_opening(&mut self, window: &mut Window, cx: &mut Context<Self>) -> Entity<Viewer> {
        let tab = self.active();
        if tab.read(cx).mode() == Mode::View {
            tab
        } else {
            self.add(Mode::View, window, cx)
        }
    }

    fn prompt_open(&mut self, window: &mut Window, cx: &mut Context<Self>) {
        let chosen = cx.prompt_for_paths(PathPromptOptions {
            files: true,
            directories: false,
            multiple: true,
            prompt: Some("Open".into()),
        });
        let window_handle = window.window_handle();
        self.prompt = Some(cx.spawn(async move |this, cx| {
            let Ok(Ok(Some(paths))) = chosen.await else {
                return;
            };
            let _ = window_handle.update(cx, |_, window, cx| {
                let _ = this.update(cx, |this, cx| this.open_paths(paths, window, cx));
            });
        }));
    }

    /// Open files: the first in the tab in front, each further one in a tab
    /// of its own.
    pub fn open_paths(&mut self, paths: Vec<std::path::PathBuf>, window: &mut Window, cx: &mut Context<Self>) {
        for (n, path) in paths.into_iter().enumerate() {
            let tab = if n == 0 {
                self.viewer_for_opening(window, cx)
            } else {
                self.add(Mode::View, window, cx)
            };
            tab.update(cx, |tab, cx| tab.open_path(path, cx));
        }
    }

    fn step(&mut self, by: isize, window: &mut Window, cx: &mut Context<Self>) {
        let n = self.tabs.len() as isize;
        let next = (self.active as isize + by).rem_euclid(n) as usize;
        self.activate(next, window, cx);
    }

    /// ⌘1 to ⌘8 go to that tab; ⌘9 to the last, as in a browser.
    fn go_to(&mut self, n: usize, window: &mut Window, cx: &mut Context<Self>) {
        let i = if n == 9 { self.tabs.len() - 1 } else { n - 1 };
        self.activate(i, window, cx);
    }
}

impl Focusable for Workspace {
    fn focus_handle(&self, _: &App) -> FocusHandle {
        self.focus.clone()
    }
}

// ---------------------------------------------------------------------------
// View
// ---------------------------------------------------------------------------

impl Workspace {
    fn tab_bar(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let count = self.tabs.len();
        div()
            .id("tab-bar")
            .flex()
            .flex_none()
            .items_end()
            .gap(px(2.))
            .h(px(TAB_BAR))
            .pl(px(LIGHTS))
            .pr(px(12.))
            .bg(theme::SUNKEN)
            .border_b_1()
            .border_color(theme::BORDER)
            // The bar is the title bar: dragging it moves the window.
            .on_mouse_down(MouseButton::Left, |_, window, _| window.start_window_move())
            .children(
                self.tabs
                    .iter()
                    .enumerate()
                    .map(|(i, tab)| self.tab(i, tab, count, cx).into_any_element()),
            )
            .child(self.add_button(cx))
            .child(div().flex_1().h_full())
    }

    fn tab(&self, i: usize, tab: &Entity<Viewer>, count: usize, cx: &mut Context<Self>) -> impl IntoElement {
        let viewer = tab.read(cx);
        let title: SharedString = viewer.title().into();
        let activity = viewer.activity();
        let folding = viewer.mode() == Mode::Fold;
        let front = i == self.active;
        let group: SharedString = format!("tab-{i}").into();

        div()
            .id(("tab", i))
            .debug_selector(move || format!("tab-{i}"))
            .group(group.clone())
            .relative()
            .flex()
            .items_center()
            .gap(px(6.))
            .flex_1()
            .min_w(px(72.))
            .max_w(px(220.))
            .h(px(31.))
            // Hangs down over the bar's hairline, so the tab in front runs
            // into the toolbar beneath it.
            .mb(px(-1.))
            .pl(px(12.))
            .pr(px(6.))
            .rounded_t(px(8.))
            .cursor_pointer()
            .when(front, |d| {
                d.bg(theme::SURFACE)
                    .border_1()
                    .border_b_0()
                    .border_color(theme::BORDER)
            })
            .when(!front, |d| d.hover(|s| s.bg(theme::PAPER)))
            .on_mouse_down(
                MouseButton::Left,
                cx.listener(move |this, _, window, cx| {
                    cx.stop_propagation();
                    this.activate(i, window, cx);
                }),
            )
            .on_mouse_down(
                MouseButton::Middle,
                cx.listener(move |this, _, window, cx| {
                    cx.stop_propagation();
                    this.close(i, window, cx);
                }),
            )
            .when(folding, |d| {
                d.child(
                    div()
                        .flex_none()
                        .px(px(4.))
                        .rounded(px(3.))
                        .bg(if front { theme::SUNKEN } else { theme::BORDER })
                        .text_size(px(8.5))
                        .text_color(theme::MUTED)
                        .child("FOLD"),
                )
            })
            .child(
                div()
                    .flex_1()
                    .min_w_0()
                    .truncate()
                    .text_size(px(11.5))
                    .text_color(if front { theme::INK } else { theme::MUTED })
                    .child(match (folding, title.strip_prefix("Fold · ")) {
                        (true, Some(rest)) => SharedString::from(rest.to_string()),
                        (true, None) if title.as_ref() == "Fold" => "New fold".into(),
                        _ => title.clone(),
                    }),
            )
            .when_some(activity, |d, fraction| d.child(spinner(i, fraction)))
            .child(
                div()
                    .id(("close-tab", i))
                    .debug_selector(move || format!("close-tab-{i}"))
                    .flex_none()
                    .size(px(18.))
                    .flex()
                    .items_center()
                    .justify_center()
                    .rounded(px(4.))
                    .text_size(px(13.))
                    .when(front, |d| d.text_color(theme::FAINT))
                    .when(!front, |d| {
                        d.text_color(gpui::transparent_black())
                            .group_hover(group, |s| s.text_color(theme::FAINT))
                    })
                    .hover(|s| s.bg(theme::SUNKEN).text_color(theme::INK))
                    .when(count > 1 || !tab.read(cx).is_blank(), |d| d.child("×"))
                    .on_mouse_down(MouseButton::Left, |_, _, cx| cx.stop_propagation())
                    .on_click(cx.listener(move |this, _, window, cx| this.close(i, window, cx))),
            )
    }

    /// New tab. With folding built in, it asks which kind.
    fn add_button(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let open = self.adding;
        div()
            .relative()
            .flex_none()
            .mb(px(4.))
            .ml(px(4.))
            .child(
                div()
                    .id("new-tab")
                    .debug_selector(|| "new-tab".into())
                    .size(px(24.))
                    .flex()
                    .items_center()
                    .justify_center()
                    .rounded(px(6.))
                    .text_size(px(16.))
                    .text_color(theme::MUTED)
                    .cursor_pointer()
                    .when(open, |d| d.bg(theme::BORDER).text_color(theme::INK))
                    .hover(|s| s.bg(theme::BORDER).text_color(theme::INK))
                    .tooltip(|_, cx| cx.new(|_| crate::ui::app::Tip("New tab  ⌘T".into())).into())
                    .on_mouse_down(MouseButton::Left, swallow)
                    .on_click(cx.listener(|this, _, window, cx| {
                        if cfg!(feature = "simplefold") {
                            this.adding = !this.adding;
                            cx.notify();
                        } else {
                            this.add(Mode::View, window, cx);
                        }
                    }))
                    .child("+"),
            )
            .when(open, |d| d.child(deferred(self.add_menu(cx)).with_priority(3)))
    }

    fn add_menu(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let entry = |id: &'static str, name: &'static str, note: &'static str, keys: &'static str| {
            div()
                .id(id)
                .debug_selector(move || id.into())
                .flex()
                .items_center()
                .justify_between()
                .gap(px(12.))
                .px(px(9.))
                .py(px(6.))
                .rounded(px(5.))
                .cursor_pointer()
                .hover(|s| s.bg(theme::SUNKEN))
                .child(
                    div()
                        .flex()
                        .flex_col()
                        .gap(px(1.))
                        .child(div().text_size(px(12.)).text_color(theme::INK).child(name))
                        .child(div().text_size(px(10.)).text_color(theme::FAINT).child(note)),
                )
                .child(
                    div()
                        .text_size(px(10.5))
                        .text_color(theme::FAINT)
                        .child(keys),
                )
        };
        menu(250.)
            .top(px(28.))
            .left(px(0.))
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.adding = false;
                cx.notify();
            }))
            .child(
                entry("new-viewer", "Viewer", "search, fetch and open structures", "⌘T")
                    .on_click(cx.listener(|this, _, window, cx| {
                        this.add(Mode::View, window, cx);
                    })),
            )
            .child(
                entry("new-fold", "Fold a protein", "predict a structure with SimpleFold", "")
                    .on_click(cx.listener(|this, _, window, cx| {
                        this.add(Mode::Fold, window, cx);
                    })),
            )
    }
}

/// A job's progress on its tab: a bar filling where the job knows how far it
/// is, and a dot that breathes where it does not.
fn spinner(i: usize, fraction: Option<f32>) -> gpui::AnyElement {
    match fraction {
        Some(f) => div()
            .flex_none()
            .w(px(22.))
            .h(px(3.))
            .rounded_full()
            .bg(theme::BORDER)
            .overflow_hidden()
            .child(div().h_full().w(relative(f.clamp(0.05, 1.0))).bg(theme::INK))
            .into_any_element(),
        None => div()
            .id(("tab-busy", i))
            .flex_none()
            .size(px(6.))
            .rounded_full()
            .bg(theme::INK)
            .with_animation(
                ("breathe", i),
                Animation::new(Duration::from_millis(1100)).repeat(),
                |el, t| el.opacity(0.25 + 0.75 * (t * std::f32::consts::PI).sin()),
            )
            .into_any_element(),
    }
}

impl Render for Workspace {
    fn render(&mut self, _: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        div()
            .flex()
            .flex_col()
            .size_full()
            .bg(theme::PAPER)
            .text_color(theme::TEXT)
            .track_focus(&self.focus)
            .on_action(cx.listener(|this, _: &OpenFile, window, cx| this.prompt_open(window, cx)))
            .on_action(cx.listener(|this, _: &NewTab, window, cx| {
                this.add(Mode::View, window, cx);
            }))
            .when(cfg!(feature = "simplefold"), |d| {
                d.on_action(cx.listener(|this, _: &NewFoldTab, window, cx| {
                    this.add(Mode::Fold, window, cx);
                }))
            })
            .on_action(cx.listener(|this, _: &CloseTab, window, cx| {
                this.close(this.active, window, cx)
            }))
            .on_action(cx.listener(|this, _: &NextTab, window, cx| this.step(1, window, cx)))
            .on_action(cx.listener(|this, _: &PreviousTab, window, cx| this.step(-1, window, cx)))
            .on_action(cx.listener(|this, _: &Tab1, window, cx| this.go_to(1, window, cx)))
            .on_action(cx.listener(|this, _: &Tab2, window, cx| this.go_to(2, window, cx)))
            .on_action(cx.listener(|this, _: &Tab3, window, cx| this.go_to(3, window, cx)))
            .on_action(cx.listener(|this, _: &Tab4, window, cx| this.go_to(4, window, cx)))
            .on_action(cx.listener(|this, _: &Tab5, window, cx| this.go_to(5, window, cx)))
            .on_action(cx.listener(|this, _: &Tab6, window, cx| this.go_to(6, window, cx)))
            .on_action(cx.listener(|this, _: &Tab7, window, cx| this.go_to(7, window, cx)))
            .on_action(cx.listener(|this, _: &Tab8, window, cx| this.go_to(8, window, cx)))
            .on_action(cx.listener(|this, _: &Tab9, window, cx| this.go_to(9, window, cx)))
            // Files dropped anywhere on the window are opened.
            .on_drop(cx.listener(|this, paths: &ExternalPaths, window, cx| {
                this.open_paths(paths.paths().to_vec(), window, cx);
            }))
            .child(self.tab_bar(cx))
            .child(
                div()
                    .flex_1()
                    .min_h_0()
                    .overflow_hidden()
                    .child(self.tabs[self.active].clone()),
            )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gpui::{Modifiers, TestAppContext, VisualTestContext};

    fn click(cx: &mut VisualTestContext, selector: &'static str) {
        let bounds = cx
            .debug_bounds(selector)
            .unwrap_or_else(|| panic!("{selector} is not on screen"));
        cx.simulate_click(bounds.center(), Modifiers::none());
        cx.run_until_parked();
    }

    fn open(cx: &mut TestAppContext) -> (Entity<Workspace>, &mut VisualTestContext) {
        cx.update(|cx| {
            crate::ui::input::bind_keys(cx);
            bind_keys(cx);
        });
        let (view, cx) = cx.add_window_view(|window, cx| Workspace::new(window, cx));
        cx.run_until_parked();
        (view, cx)
    }

    fn tabs(view: &Entity<Workspace>, cx: &mut VisualTestContext) -> (usize, usize) {
        view.update(cx, |this, _| (this.tabs.len(), this.active))
    }

    #[gpui::test]
    fn tabs_open_switch_and_close(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        assert_eq!(tabs(&view, cx), (1, 0));

        cx.simulate_keystrokes("cmd-t");
        cx.simulate_keystrokes("cmd-t");
        assert_eq!(tabs(&view, cx), (3, 2));

        // Each tab holds its own structure.
        view.update(cx, |this, cx| {
            let first = this.tabs[0].clone();
            first.update(cx, |tab, cx| tab.open_path("samples/caffeine.sdf".into(), cx));
        });
        cx.run_until_parked();

        click(cx, "tab-0");
        assert_eq!(tabs(&view, cx), (3, 0));
        cx.simulate_keystrokes("cmd-}");
        assert_eq!(tabs(&view, cx), (3, 1));
        cx.simulate_keystrokes("cmd-9");
        assert_eq!(tabs(&view, cx), (3, 2));
        cx.simulate_keystrokes("cmd-1");
        assert_eq!(tabs(&view, cx), (3, 0));

        click(cx, "close-tab-0");
        assert_eq!(tabs(&view, cx), (2, 0));
        cx.simulate_keystrokes("cmd-w");
        cx.simulate_keystrokes("cmd-w");
        // The last tab is emptied rather than closed.
        assert_eq!(tabs(&view, cx), (1, 0));
    }

    #[cfg(feature = "simplefold")]
    #[gpui::test]
    fn folding_is_its_own_tab_from_the_plus_menu(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        assert!(cx.debug_bounds("model").is_none(), "a viewer shows fold settings");
        click(cx, "new-tab");
        click(cx, "new-fold");
        let mode = view.update(cx, |this, cx| this.active().read(cx).mode());
        assert_eq!(mode, Mode::Fold);
        assert!(cx.debug_bounds("model").is_some(), "the fold tab has no model setting");
        // Something opened from the window goes to a viewer, not the fold tab.
        cx.update(|window, cx| {
            view.update(cx, |this, cx| this.open_paths(vec!["samples/caffeine.sdf".into()], window, cx))
        });
        let (count, mode) = view.update(cx, |this, cx| (this.tabs.len(), this.active().read(cx).mode()));
        assert_eq!((count, mode), (3, Mode::View));
    }
}
