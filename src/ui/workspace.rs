//! The window: a row of tabs across the top, and the tab in front beneath.
//!
//! Each tab is a viewer of its own, with its own structure, camera and job,
//! so a protein, a trajectory and a molecule can be open side by side and a
//! fold can run in one tab while the others are used. A tab for folding is
//! opened only by asking for one, from the + menu.

use std::time::Duration;

use gpui::{
    actions, deferred, div, prelude::*, px, relative, Animation, AnimationExt, App, Context,
    Entity, ExternalPaths, FocusHandle, Focusable, Menu, MenuItem, MouseButton, OsAction,
    PathPromptOptions, SharedString, SystemMenuType, Task, Window,
};

use crate::render::{Repr, Style};

use crate::settings;
use crate::theme;
use crate::ui::app::{divide, item, menu, swallow, Mode, Tip, Viewer};
use crate::update::{self, Status};

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
        CheckForUpdates,
        About,
        Quit,
        HideApp,
        HideOthers,
        ShowAll,
        SaveImage,
        ToggleDark,
        ShowCartoon,
        ShowTube,
        ShowBallAndStick,
        ShowLicorice,
        ShowSpacefill,
        ShowSurface,
        ShowShortcuts,
        Undo,
        Redo,
        ToggleConsole,
        StyleIllustrative,
        StyleFlat,
        StyleSoft,
        StyleStudio,
        Recentre,
        ResetView,
        ToggleFullScreen,
        Minimize,
        Zoom,
        OpenHelp,
        ReportIssue,
    ]
);

const REPO: &str = "https://github.com/odexion/profold";

/// Height of the tab row, which is also the window's title bar.
const TAB_BAR: f32 = 40.0;
/// Room for the traffic lights at the left of the bar. Elsewhere the
/// window's buttons are at the right, drawn by the bar itself.
const LIGHTS: f32 = if cfg!(target_os = "macos") { 80.0 } else { 8.0 };

pub fn bind_keys(cx: &mut App) {
    cx.bind_keys([
        gpui::KeyBinding::new("secondary-o", OpenFile, None),
        gpui::KeyBinding::new("secondary-t", NewTab, None),
        gpui::KeyBinding::new("secondary-shift-t", NewFoldTab, None),
        gpui::KeyBinding::new("secondary-w", CloseTab, None),
        gpui::KeyBinding::new("ctrl-tab", NextTab, None),
        gpui::KeyBinding::new("ctrl-shift-tab", PreviousTab, None),
        gpui::KeyBinding::new("secondary-}", NextTab, None),
        gpui::KeyBinding::new("secondary-{", PreviousTab, None),
        gpui::KeyBinding::new("secondary-1", Tab1, None),
        gpui::KeyBinding::new("secondary-2", Tab2, None),
        gpui::KeyBinding::new("secondary-3", Tab3, None),
        gpui::KeyBinding::new("secondary-4", Tab4, None),
        gpui::KeyBinding::new("secondary-5", Tab5, None),
        gpui::KeyBinding::new("secondary-6", Tab6, None),
        gpui::KeyBinding::new("secondary-7", Tab7, None),
        gpui::KeyBinding::new("secondary-8", Tab8, None),
        gpui::KeyBinding::new("secondary-9", Tab9, None),
        gpui::KeyBinding::new("secondary-q", Quit, None),
        #[cfg(not(target_os = "macos"))]
        gpui::KeyBinding::new("alt-f4", Quit, None),
        gpui::KeyBinding::new("secondary-shift-s", SaveImage, None),
        gpui::KeyBinding::new("secondary-shift-d", ToggleDark, None),
        // ⌘1 to ⌘9 are the tabs, so the views take ⌥⌘ and the styles ⌃⌘.
        gpui::KeyBinding::new("alt-secondary-1", ShowCartoon, None),
        gpui::KeyBinding::new("alt-secondary-2", ShowTube, None),
        gpui::KeyBinding::new("alt-secondary-3", ShowBallAndStick, None),
        gpui::KeyBinding::new("alt-secondary-4", ShowLicorice, None),
        gpui::KeyBinding::new("alt-secondary-5", ShowSpacefill, None),
        gpui::KeyBinding::new("alt-secondary-6", ShowSurface, None),
        gpui::KeyBinding::new("secondary-/", ShowShortcuts, None),
        gpui::KeyBinding::new("secondary-z", Undo, None),
        gpui::KeyBinding::new("secondary-k", ToggleConsole, None),
        gpui::KeyBinding::new("secondary-shift-z", Redo, None),
        gpui::KeyBinding::new("ctrl-secondary-1", StyleIllustrative, None),
        gpui::KeyBinding::new("ctrl-secondary-2", StyleFlat, None),
        gpui::KeyBinding::new("ctrl-secondary-3", StyleSoft, None),
        gpui::KeyBinding::new("ctrl-secondary-4", StyleStudio, None),
        gpui::KeyBinding::new("secondary-0", ResetView, None),
    ]);
    #[cfg(target_os = "macos")]
    cx.bind_keys([
        gpui::KeyBinding::new("cmd-h", HideApp, None),
        gpui::KeyBinding::new("alt-cmd-h", HideOthers, None),
        gpui::KeyBinding::new("cmd-m", Minimize, None),
        gpui::KeyBinding::new("ctrl-cmd-f", ToggleFullScreen, None),
    ]);
    cx.on_action(|_: &Quit, cx: &mut App| cx.quit());
    cx.on_action(|_: &HideApp, cx: &mut App| cx.hide());
    cx.on_action(|_: &HideOthers, cx: &mut App| cx.hide_other_apps());
    cx.on_action(|_: &ShowAll, cx: &mut App| cx.unhide_other_apps());
    cx.on_action(|_: &OpenHelp, cx: &mut App| cx.open_url(&format!("{REPO}#readme")));
    cx.on_action(|_: &ReportIssue, cx: &mut App| cx.open_url(&format!("{REPO}/issues/new")));
}

/// The menu bar. Built again when the theme changes, so its item names the
/// way it would go.
pub fn menus() -> Vec<Menu> {
    use crate::ui::input;
    vec![
        Menu {
            name: "profold".into(),
            items: vec![
                MenuItem::action("About profold", About),
                MenuItem::action("Check for Updates…", CheckForUpdates),
                MenuItem::separator(),
                MenuItem::os_submenu("Services", SystemMenuType::Services),
                MenuItem::separator(),
                MenuItem::action("Hide profold", HideApp),
                MenuItem::action("Hide Others", HideOthers),
                MenuItem::action("Show All", ShowAll),
                MenuItem::separator(),
                MenuItem::action("Quit profold", Quit),
            ],
        },
        Menu {
            name: "File".into(),
            items: {
                let mut items = vec![MenuItem::action("New Tab", NewTab)];
                if cfg!(feature = "simplefold") {
                    items.push(MenuItem::action("New Fold Tab", NewFoldTab));
                }
                items.extend([
                    MenuItem::action("Open…", OpenFile),
                    MenuItem::separator(),
                    MenuItem::action("Save Image…", SaveImage),
                    MenuItem::separator(),
                    MenuItem::action("Close Tab", CloseTab),
                ]);
                items
            },
        },
        Menu {
            name: "Edit".into(),
            items: vec![
                MenuItem::action("Undo", Undo),
                MenuItem::action("Redo", Redo),
                MenuItem::separator(),
                MenuItem::os_action("Cut", input::Cut, OsAction::Cut),
                MenuItem::os_action("Copy", input::Copy, OsAction::Copy),
                MenuItem::os_action("Paste", input::Paste, OsAction::Paste),
                MenuItem::os_action("Select All", input::SelectAll, OsAction::SelectAll),
            ],
        },
        Menu {
            name: "View".into(),
            items: vec![
                MenuItem::action(if theme::dark() { "Light Mode" } else { "Dark Mode" }, ToggleDark),
                MenuItem::separator(),
                MenuItem::action(Repr::Cartoon.label(), ShowCartoon),
                MenuItem::action(Repr::Tube.label(), ShowTube),
                MenuItem::action(Repr::BallAndStick.label(), ShowBallAndStick),
                MenuItem::action(Repr::Licorice.label(), ShowLicorice),
                MenuItem::action(Repr::Spacefill.label(), ShowSpacefill),
                MenuItem::action(Repr::Surface.label(), ShowSurface),
                MenuItem::separator(),
                MenuItem::action(Style::Illustrative.label(), StyleIllustrative),
                MenuItem::action(Style::Flat.label(), StyleFlat),
                MenuItem::action(Style::Soft.label(), StyleSoft),
                MenuItem::action(Style::Studio.label(), StyleStudio),
                MenuItem::separator(),
                MenuItem::action("Centre", Recentre),
                MenuItem::action("Reset View", ResetView),
                MenuItem::separator(),
                MenuItem::action("Keyboard Shortcuts", ShowShortcuts),
                MenuItem::action("Console", ToggleConsole),
                MenuItem::action("Enter Full Screen", ToggleFullScreen),
            ],
        },
        Menu {
            name: "Window".into(),
            items: vec![
                MenuItem::action("Minimize", Minimize),
                MenuItem::action("Zoom", Zoom),
                MenuItem::separator(),
                MenuItem::action("Next Tab", NextTab),
                MenuItem::action("Previous Tab", PreviousTab),
            ],
        },
        Menu {
            name: "Help".into(),
            items: vec![
                MenuItem::action("profold Help", OpenHelp),
                MenuItem::action("Report an Issue…", ReportIssue),
            ],
        },
    ]
}

pub struct Workspace {
    tabs: Vec<Entity<Viewer>>,
    active: usize,
    /// Whether the + menu is open.
    adding: bool,
    focus: FocusHandle,
    /// The open dialog, while it is up.
    prompt: Option<Task<()>>,
    /// The tabs as last written to `session.json`, so it is written only
    /// when they change.
    saved: String,
    /// Where the updater has got to, written from its own thread.
    update: update::Shared,
    /// Whether the version menu is open.
    about: bool,

    /// Redraws the bar while the updater works.
    update_poll: Option<Task<()>>,
}

/// A tab being dragged along the bar, by where it was.
#[derive(Clone, Copy)]
struct DraggedTab(usize);

/// What follows the pointer while a tab is dragged: its name on a card.
struct TabGhost(SharedString);

impl Render for TabGhost {
    fn render(&mut self, _: &mut Window, _: &mut Context<Self>) -> impl IntoElement {
        div()
            .px(px(12.))
            .py(px(6.))
            .rounded(px(7.))
            .bg(theme::surface())
            .border_1()
            .border_color(theme::border())
            .shadow_md()
            .text_size(px(11.5))
            .text_color(theme::ink())
            .child(self.0.clone())
    }
}

impl Workspace {
    pub fn new(window: &mut Window, cx: &mut Context<Self>) -> Self {
        let mut this = Workspace {
            tabs: Vec::new(),
            active: 0,
            adding: false,
            focus: cx.focus_handle(),
            prompt: None,
            saved: String::new(),
            update: std::sync::Arc::new(std::sync::Mutex::new(Status::Idle)),
            about: false,
            update_poll: None,
        };
        theme::set_dark(settings::Settings::load().dark);
        cx.set_menus(menus());
        this.add(Mode::View, window, cx);
        if update::due() {
            this.check_for_updates(cx);
        }
        this
    }

    /// Turn the window light or dark.
    fn toggle_theme(&mut self, cx: &mut Context<Self>) {
        theme::set_dark(!theme::dark());
        cx.set_menus(menus());
        self.retheme(cx);
    }

    fn show(&mut self, repr: Repr, cx: &mut Context<Self>) {
        self.active().update(cx, |tab, cx| tab.set_repr(repr, cx));
    }

    fn draw_in(&mut self, style: Style, cx: &mut Context<Self>) {
        self.active().update(cx, |tab, cx| tab.set_style(style, cx));
    }

    /// Redraw every tab on the new page.
    fn retheme(&mut self, cx: &mut Context<Self>) {
        for tab in &self.tabs {
            tab.update(cx, |tab, cx| tab.retheme(cx));
        }
        cx.notify();
    }

    fn theme_button(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let dark = theme::dark();
        div()
            .id("theme")
            .debug_selector(|| "theme".into())
            .flex_none()
            .mb(px(6.))
            .px(px(7.))
            .py(px(3.))
            .rounded(px(5.))
            .text_size(px(11.))
            .text_color(theme::faint())
            .cursor_pointer()
            .hover(|s| s.bg(theme::border()).text_color(theme::ink()))
            .on_mouse_down(MouseButton::Left, swallow)
            .tooltip(move |_, cx| {
                cx.new(|_| Tip(if dark { "Light mode" } else { "Dark mode" }.into())).into()
            })
            .on_click(cx.listener(|this, _, _, cx| this.toggle_theme(cx)))
            .child(if dark { "☀" } else { "☾" })
    }

    /// Look for a new release, and install it where this build can. It runs
    /// on a thread of its own; the bar is redrawn until it is done.
    pub fn check_for_updates(&mut self, cx: &mut Context<Self>) {
        {
            let mut status = self.update.lock().unwrap();
            // Installed already: only a restart is wanted.
            if status.busy() || matches!(*status, Status::Ready { .. }) {
                return;
            }
            *status = Status::Checking;
        }
        let status = self.update.clone();
        std::thread::spawn(move || update::run(status));
        let status = self.update.clone();
        self.update_poll = Some(cx.spawn(async move |this, cx| loop {
            cx.background_executor().timer(Duration::from_millis(120)).await;
            let busy = status.lock().unwrap().busy();
            if this.update(cx, |_, cx| cx.notify()).is_err() || !busy {
                return;
            }
        }));
        cx.notify();
    }

    /// Open a new tab beside the one in front, and bring it forward.
    pub fn add(&mut self, mode: Mode, window: &mut Window, cx: &mut Context<Self>) -> Entity<Viewer> {
        let tab = cx.new(|cx| Viewer::new(mode, window, cx));
        // The tab's name and activity show in the bar, which is drawn here,
        // and what it has open goes into the session.
        cx.observe(&tab, |this, _, cx| {
            this.save(cx);
            cx.notify();
        })
        .detach();
        if let Some(behind) = self.tabs.get(self.active).cloned() {
            behind.update(cx, |tab, cx| tab.pause(cx));
        }
        let at = if self.tabs.is_empty() { 0 } else { self.active + 1 };
        self.tabs.insert(at, tab.clone());
        self.active = at;
        self.adding = false;
        tab.update(cx, |tab, cx| tab.focus_input(window, cx));
        cx.notify();
        tab
    }

    /// Move a tab to another place in the bar, keeping the one in front.
    fn move_tab(&mut self, from: usize, to: usize, cx: &mut Context<Self>) {
        if from == to || from >= self.tabs.len() || to >= self.tabs.len() {
            return;
        }
        let front = self.tabs[self.active].clone();
        let tab = self.tabs.remove(from);
        self.tabs.insert(to, tab);
        self.active = self.tabs.iter().position(|t| *t == front).unwrap_or(0);
        self.save(cx);
        cx.notify();
    }

    /// What each tab has open, as JSON for `session.json`.
    fn session(&self, cx: &App) -> String {
        let tabs: Vec<serde_json::Value> = self
            .tabs
            .iter()
            .map(|tab| {
                let tab = tab.read(cx);
                serde_json::json!({
                    "mode": match tab.mode() { Mode::View => "view", Mode::Fold => "fold" },
                    "open": tab.origin(),
                })
            })
            .collect();
        serde_json::to_string_pretty(&serde_json::json!({ "tabs": tabs, "active": self.active }))
            .unwrap_or_default()
    }

    /// Write the tabs out if they have changed, to be reopened next launch.
    fn save(&mut self, cx: &App) {
        let json = self.session(cx);
        if json == self.saved {
            return;
        }
        if let Some(path) = settings::file("session.json") {
            if let Some(dir) = path.parent() {
                let _ = std::fs::create_dir_all(dir);
            }
            let _ = std::fs::write(path, &json);
        }
        self.saved = json;
    }

    /// Reopen the tabs of the last session. A tab whose file has gone, or
    /// whose download fails, opens empty. Returns whether there was one.
    pub fn restore(&mut self, window: &mut Window, cx: &mut Context<Self>) -> bool {
        let Some(json) = settings::file("session.json")
            .and_then(|p| std::fs::read_to_string(p).ok())
            .and_then(|t| serde_json::from_str::<serde_json::Value>(&t).ok())
        else {
            return false;
        };
        let entries = json.get("tabs").and_then(|t| t.as_array()).cloned().unwrap_or_default();
        // A lone empty tab is what a new window has anyway.
        if entries.is_empty() || (entries.len() == 1 && entries[0].get("open").is_none_or(|o| o.is_null())) {
            return false;
        }
        self.tabs.clear();
        self.active = 0;
        for entry in &entries {
            let mode = match entry.get("mode").and_then(|m| m.as_str()) {
                Some("fold") if cfg!(feature = "simplefold") => Mode::Fold,
                _ => Mode::View,
            };
            let tab = self.add(mode, window, cx);
            if let Some(origin) = entry.get("open").and_then(|o| o.as_str()) {
                let origin = origin.to_string();
                tab.update(cx, |tab, cx| tab.reopen(&origin, cx));
            }
        }
        let active = json.get("active").and_then(|a| a.as_u64()).unwrap_or(0) as usize;
        self.activate(active.min(self.tabs.len() - 1), window, cx);
        true
    }

    /// The tab in front.
    pub fn active(&self) -> Entity<Viewer> {
        self.tabs[self.active].clone()
    }

    fn activate(&mut self, i: usize, window: &mut Window, cx: &mut Context<Self>) {
        if i >= self.tabs.len() {
            return;
        }
        if i != self.active {
            if let Some(behind) = self.tabs.get(self.active).cloned() {
                behind.update(cx, |tab, cx| tab.pause(cx));
            }
        }
        self.active = i;
        self.adding = false;
        self.save(cx);
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
        self.save(cx);
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
    fn tab_bar(&self, window: &Window, cx: &mut Context<Self>) -> impl IntoElement {
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
            .relative()
            .bg(theme::sunken())
            // The bar is the title bar: dragging it moves the window.
            .on_mouse_down(MouseButton::Left, |_, window, _| window.start_window_move())
            // The hairline is drawn first, under the tabs, so the tab in front
            // covers it and runs into the toolbar beneath. A border would be
            // painted over the tabs.
            .child(
                div()
                    .absolute()
                    .left_0()
                    .right_0()
                    .bottom_0()
                    .h(px(1.))
                    .bg(theme::border()),
            )
            .children(
                self.tabs
                    .iter()
                    .enumerate()
                    .map(|(i, tab)| self.tab(i, tab, count, cx).into_any_element()),
            )
            .child(self.add_button(cx))
            .child(
                div()
                    .flex_1()
                    .h_full()
                    .window_control_area(gpui::WindowControlArea::Drag),
            )
            .child(self.theme_button(cx))
            .child(self.version_button(cx))
            .when(cfg!(not(target_os = "macos")), |d| d.children(window_buttons(window)))
    }

    /// The version, at the right of the bar, and the update when there is
    /// one waiting for a restart.
    fn version_button(&self, cx: &mut Context<Self>) -> impl IntoElement {
        let status = self.update.lock().unwrap().clone();
        let open = self.about;
        let label = match &status {
            Status::Checking => "Checking…".to_string(),
            Status::Downloading { fraction: Some(f), .. } => format!("Updating {:.0}%", f * 100.0),
            Status::Downloading { .. } | Status::Installing { .. } => "Updating…".to_string(),
            _ => format!("v{}", update::VERSION),
        };
        div()
            .relative()
            .flex()
            .flex_none()
            .items_center()
            .gap(px(6.))
            .mb(px(6.))
            .when(matches!(status, Status::Ready { .. }), |d| {
                d.child(
                    div()
                        .id("restart-to-update")
                        .debug_selector(|| "restart-to-update".into())
                        .px(px(10.))
                        .py(px(3.))
                        .rounded(px(6.))
                        .bg(theme::ink())
                        .text_color(theme::on_ink())
                        .text_size(px(11.))
                        .cursor_pointer()
                        .hover(|s| s.opacity(0.85))
                        .on_mouse_down(MouseButton::Left, swallow)
                        .on_click(|_, _, cx| update::restart(cx))
                        .child("Restart to update"),
                )
            })
            .when(matches!(status, Status::Available { .. }), |d| {
                d.child(div().size(px(6.)).rounded_full().bg(theme::ink()))
            })
            .child(
                div()
                    .id("version")
                    .debug_selector(|| "version".into())
                    .px(px(7.))
                    .py(px(3.))
                    .rounded(px(5.))
                    .text_size(px(10.5))
                    .text_color(theme::faint())
                    .cursor_pointer()
                    .when(open, |d| d.bg(theme::border()).text_color(theme::ink()))
                    .hover(|s| s.bg(theme::border()).text_color(theme::ink()))
                    .on_mouse_down(MouseButton::Left, swallow)
                    .on_click(cx.listener(|this, _, _, cx| {
                        this.about = !this.about;
                        cx.notify();
                    }))
                    .child(label),
            )
            .when(open, |d| d.child(deferred(self.about_menu(&status, cx)).with_priority(3)))
    }

    fn about_menu(&self, status: &Status, cx: &mut Context<Self>) -> impl IntoElement {
        let line = |text: String, color| {
            div().px(px(8.)).text_size(px(10.5)).text_color(color).child(text)
        };
        let (note, color) = match status {
            Status::Idle => (None, theme::muted()),
            Status::Checking => (Some("Checking for updates…".to_string()), theme::muted()),
            Status::UpToDate => (Some("You have the latest version.".to_string()), theme::muted()),
            Status::Available { version, .. } => (
                Some(if update::TARGET.is_none() {
                    format!("{version} is out. This build is from source, so it is not replaced.")
                } else {
                    format!("{version} is out, but this copy cannot be replaced; use the installer.")
                }),
                theme::ink(),
            ),
            Status::Downloading { version, fraction } => (
                Some(match fraction {
                    Some(f) => format!("Downloading {version}… {:.0}%", f * 100.0),
                    None => format!("Downloading {version}…"),
                }),
                theme::muted(),
            ),
            Status::Installing { version } => (Some(format!("Installing {version}…")), theme::muted()),
            Status::Ready { version } => (
                Some(format!("{version} is installed and starts next time.")),
                theme::ink(),
            ),
            Status::Failed(e) => (Some(format!("Could not update: {e}")), theme::danger()),
        };
        let busy = status.busy();
        let automatic = update::automatic();
        menu(270.)
            .top(px(26.))
            .left_auto()
            .right(px(0.))
            .on_mouse_down_out(cx.listener(|this, _, _, cx| {
                this.about = false;
                cx.notify();
            }))
            .child(
                div()
                    .px(px(8.))
                    .pt(px(5.))
                    .flex()
                    .items_baseline()
                    .gap(px(6.))
                    .child(div().text_size(px(13.)).text_color(theme::ink()).child("profold"))
                    .child(
                        div()
                            .text_size(px(11.))
                            .text_color(theme::muted())
                            .child(format!("version {}", update::VERSION)),
                    ),
            )
            .child(line(
                update::TARGET.unwrap_or("built from source").to_string(),
                theme::faint(),
            ))
            .when_some(note, |d, note| d.child(div().pt(px(4.)).child(line(note, color))))
            .child(divide())
            .when(matches!(status, Status::Ready { .. }), |d| {
                d.child(
                    item("about-restart", "Restart to update".into(), false)
                        .on_click(|_, _, cx| update::restart(cx)),
                )
            })
            .when_some(
                match status {
                    Status::Available { page, .. } => Some(page.clone()),
                    _ => None,
                },
                |d, page| {
                    d.child(
                        item("about-release", "Open the release page".into(), false)
                            .on_click(move |_, _, cx| cx.open_url(&page)),
                    )
                },
            )
            .child(
                item("check-for-updates", "Check for updates".into(), false)
                    .debug_selector(|| "check-for-updates".into())
                    .when(busy, |d| d.opacity(0.4))
                    .on_click(cx.listener(|this, _, _, cx| this.check_for_updates(cx))),
            )
            .child(
                item("check-automatically", "Check automatically".into(), automatic).on_click(
                    cx.listener(move |_, _, _, cx| {
                        update::set_automatic(!automatic);
                        cx.notify();
                    }),
                ),
            )
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
            .pl(px(12.))
            .pr(px(6.))
            .rounded_t(px(8.))
            .cursor_pointer()
            .when(front, |d| {
                d.bg(theme::surface())
                    .border_1()
                    .border_b_0()
                    .border_color(theme::border())
            })
            .when(!front, |d| d.hover(|s| s.bg(theme::paper())))
            // Dragged along the bar to reorder; dropped on another tab, it
            // takes that tab's place.
            .on_drag(DraggedTab(i), {
                let title = title.clone();
                move |_, _, _, cx| cx.new(|_| TabGhost(title.clone()))
            })
            .drag_over::<DraggedTab>(|style, _, _, _| style.bg(theme::border()))
            .on_drop(cx.listener(move |this, dragged: &DraggedTab, _, cx| {
                this.move_tab(dragged.0, i, cx);
            }))
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
                        .bg(if front { theme::sunken() } else { theme::border() })
                        .text_size(px(8.5))
                        .text_color(theme::muted())
                        .child("FOLD"),
                )
            })
            .child(
                div()
                    .flex_1()
                    .min_w_0()
                    .truncate()
                    .text_size(px(11.5))
                    .text_color(if front { theme::ink() } else { theme::muted() })
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
                    .when(front, |d| d.text_color(theme::faint()))
                    .when(!front, |d| {
                        d.text_color(gpui::transparent_black())
                            .group_hover(group, |s| s.text_color(theme::faint()))
                    })
                    .hover(|s| s.bg(theme::sunken()).text_color(theme::ink()))
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
                    .text_color(theme::muted())
                    .cursor_pointer()
                    .when(open, |d| d.bg(theme::border()).text_color(theme::ink()))
                    .hover(|s| s.bg(theme::border()).text_color(theme::ink()))
                    .tooltip(|_, cx| cx.new(|_| crate::ui::app::Tip(format!("New tab  {}", theme::keys("⌘T")).into())).into())
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
        let entry = |id: &'static str, name: &'static str, note: &'static str, keys: String| {
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
                .hover(|s| s.bg(theme::sunken()))
                .child(
                    div()
                        .flex()
                        .flex_col()
                        .gap(px(1.))
                        .child(div().text_size(px(12.)).text_color(theme::ink()).child(name))
                        .child(div().text_size(px(10.)).text_color(theme::faint()).child(note)),
                )
                .child(
                    div()
                        .text_size(px(10.5))
                        .text_color(theme::faint())
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
                entry("new-viewer", "Viewer", "search, fetch and open structures", theme::keys("⌘T"))
                    .on_click(cx.listener(|this, _, window, cx| {
                        this.add(Mode::View, window, cx);
                    })),
            )
            .child(
                entry("new-fold", "Fold a protein", "predict a structure with SimpleFold", theme::keys("⇧⌘T"))
                    .on_click(cx.listener(|this, _, window, cx| {
                        this.add(Mode::Fold, window, cx);
                    })),
            )
    }
}

/// Minimise, maximise and close, where the system does not draw them: always
/// on Windows, whose title bar the tab bar replaces, and on Linux when the
/// compositor leaves the frame to the app.
fn window_buttons(window: &Window) -> Option<impl IntoElement> {
    let drawn_here = cfg!(windows)
        || (cfg!(target_os = "linux")
            && matches!(window.window_decorations(), gpui::Decorations::Client { .. }));
    if !drawn_here {
        return None;
    }
    use gpui::WindowControlArea as Area;
    let button = |id: &'static str, glyph: &'static str, area: Area| {
        div()
            .id(id)
            .w(px(40.))
            .h(px(TAB_BAR))
            .flex()
            .items_center()
            .justify_center()
            .text_size(px(12.))
            .text_color(theme::muted())
            .window_control_area(area)
            .when(area == Area::Close, |d| {
                d.hover(|s| s.bg(theme::danger()).text_color(theme::on_ink()))
            })
            .when(area != Area::Close, |d| d.hover(|s| s.bg(theme::border()).text_color(theme::ink())))
            .on_mouse_down(MouseButton::Left, swallow)
            // Windows acts on these areas itself; Linux needs telling.
            .when(cfg!(target_os = "linux"), |d| {
                d.on_click(move |_, window, _| match area {
                    Area::Min => window.minimize_window(),
                    Area::Max => window.zoom_window(),
                    _ => window.remove_window(),
                })
            })
            .child(glyph)
    };
    let maximized = window.is_maximized();
    Some(
        div()
            .flex()
            .flex_none()
            .ml(px(6.))
            .mr(px(-12.))
            .child(button("window-min", "─", Area::Min))
            .child(button("window-max", if maximized { "❐" } else { "□" }, Area::Max))
            .child(button("window-close", "✕", Area::Close)),
    )
}

/// Which edge of the window a press at `p` is on, for a window that draws
/// its own frame.
fn resize_edge(p: gpui::Point<gpui::Pixels>, size: gpui::Size<gpui::Pixels>) -> Option<gpui::ResizeEdge> {
    use gpui::ResizeEdge as E;
    const GRIP: f32 = 6.0;
    let (x, y) = (f32::from(p.x), f32::from(p.y));
    let (w, h) = (f32::from(size.width), f32::from(size.height));
    let (left, right, top, bottom) = (x < GRIP, x > w - GRIP, y < GRIP, y > h - GRIP);
    Some(match (left, right, top, bottom) {
        (true, _, true, _) => E::TopLeft,
        (_, true, true, _) => E::TopRight,
        (true, _, _, true) => E::BottomLeft,
        (_, true, _, true) => E::BottomRight,
        (true, ..) => E::Left,
        (_, true, ..) => E::Right,
        (_, _, true, _) => E::Top,
        (.., true) => E::Bottom,
        _ => return None,
    })
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
            .bg(theme::border())
            .overflow_hidden()
            .child(div().h_full().w(relative(f.clamp(0.05, 1.0))).bg(theme::ink()))
            .into_any_element(),
        None => div()
            .id(("tab-busy", i))
            .flex_none()
            .size(px(6.))
            .rounded_full()
            .bg(theme::ink())
            .with_animation(
                ("breathe", i),
                Animation::new(Duration::from_millis(1100)).repeat(),
                |el, t| el.opacity(0.25 + 0.75 * (t * std::f32::consts::PI).sin()),
            )
            .into_any_element(),
    }
}

impl Render for Workspace {
    fn render(&mut self, window: &mut Window, cx: &mut Context<Self>) -> impl IntoElement {
        let client_side = matches!(window.window_decorations(), gpui::Decorations::Client { .. });
        div()
            .flex()
            .flex_col()
            .size_full()
            .bg(theme::paper())
            .text_color(theme::text())
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
            .on_action(cx.listener(|this, _: &CheckForUpdates, _, cx| {
                this.about = true;
                this.check_for_updates(cx);
            }))
            .on_action(cx.listener(|this, _: &About, _, cx| {
                this.about = true;
                cx.notify();
            }))
            .on_action(cx.listener(|this, _: &ToggleDark, _, cx| this.toggle_theme(cx)))
            .on_action(cx.listener(|this, _: &SaveImage, _, cx| {
                this.active().update(cx, |tab, cx| tab.save_image(2, cx))
            }))
            .on_action(cx.listener(|this, _: &ShowCartoon, _, cx| this.show(Repr::Cartoon, cx)))
            .on_action(cx.listener(|this, _: &ShowTube, _, cx| this.show(Repr::Tube, cx)))
            .on_action(cx.listener(|this, _: &ShowBallAndStick, _, cx| this.show(Repr::BallAndStick, cx)))
            .on_action(cx.listener(|this, _: &ShowLicorice, _, cx| this.show(Repr::Licorice, cx)))
            .on_action(cx.listener(|this, _: &ShowSpacefill, _, cx| this.show(Repr::Spacefill, cx)))
            .on_action(cx.listener(|this, _: &ShowSurface, _, cx| this.show(Repr::Surface, cx)))
            .on_action(cx.listener(|this, _: &Undo, _, cx| this.active().update(cx, |tab, cx| tab.undo(false, cx))))
            .on_action(cx.listener(|this, _: &ToggleConsole, window, cx| {
                this.active().update(cx, |tab, cx| tab.toggle_console(window, cx))
            }))
            .on_action(cx.listener(|this, _: &Redo, _, cx| this.active().update(cx, |tab, cx| tab.undo(true, cx))))
            .on_action(cx.listener(|this, _: &ShowShortcuts, _, cx| {
                this.active().update(cx, |tab, cx| tab.toggle_shortcuts(cx))
            }))
            .on_action(cx.listener(|this, _: &StyleIllustrative, _, cx| this.draw_in(Style::Illustrative, cx)))
            .on_action(cx.listener(|this, _: &StyleFlat, _, cx| this.draw_in(Style::Flat, cx)))
            .on_action(cx.listener(|this, _: &StyleSoft, _, cx| this.draw_in(Style::Soft, cx)))
            .on_action(cx.listener(|this, _: &StyleStudio, _, cx| this.draw_in(Style::Studio, cx)))
            .on_action(cx.listener(|this, _: &Recentre, _, cx| {
                this.active().update(cx, |tab, cx| tab.recentre(cx))
            }))
            .on_action(cx.listener(|this, _: &ResetView, _, cx| {
                this.active().update(cx, |tab, cx| tab.reset_view(cx))
            }))
            .on_action(|_: &ToggleFullScreen, window, _| window.toggle_fullscreen())
            .on_action(|_: &Minimize, window, _| window.minimize_window())
            .on_action(|_: &Zoom, window, _| window.zoom_window())
            // A window that draws its own frame is resized from its edges.
            .when(client_side, |d| {
                d.capture_any_mouse_down(|e, window, cx| {
                    if e.button != MouseButton::Left {
                        return;
                    }
                    if let Some(edge) = resize_edge(e.position, window.viewport_size()) {
                        cx.stop_propagation();
                        window.start_window_resize(edge);
                    }
                })
            })
            // Files dropped anywhere on the window are opened.
            .on_drop(cx.listener(|this, paths: &ExternalPaths, window, cx| {
                this.open_paths(paths.paths().to_vec(), window, cx);
            }))
            .child(self.tab_bar(window, cx))
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

    #[gpui::test]
    fn the_version_opens_a_menu_to_check_for_updates(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        // Tests never check on their own.
        assert_eq!(*view.read_with(cx, |w, _| w.update.clone()).lock().unwrap(), Status::Idle);
        assert!(cx.debug_bounds("restart-to-update").is_none());

        click(cx, "version");
        assert!(view.read_with(cx, |w, _| w.about));
        assert!(cx.debug_bounds("check-for-updates").is_some());

        // Once an update is in place, the bar offers the restart.
        view.update(cx, |w, cx| {
            *w.update.lock().unwrap() = Status::Ready { version: "9.9.9".into() };
            cx.notify();
        });
        cx.run_until_parked();
        assert!(cx.debug_bounds("restart-to-update").is_some());
    }

    #[gpui::test]
    fn the_moon_turns_the_window_dark_and_back(cx: &mut TestAppContext) {
        let (_view, cx) = open(cx);
        assert!(!theme::dark());
        click(cx, "theme");
        assert!(theme::dark());
        click(cx, "theme");
        assert!(!theme::dark());
        // And from the keyboard, as the View menu does it.
        cx.simulate_keystrokes("secondary-shift-d");
        assert!(theme::dark());
        cx.simulate_keystrokes("secondary-shift-d");
        assert!(!theme::dark());
    }

    fn tabs(view: &Entity<Workspace>, cx: &mut VisualTestContext) -> (usize, usize) {
        view.update(cx, |this, _| (this.tabs.len(), this.active))
    }

    #[gpui::test]
    fn tabs_open_switch_and_close(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        assert_eq!(tabs(&view, cx), (1, 0));

        cx.simulate_keystrokes("secondary-t");
        cx.simulate_keystrokes("secondary-t");
        assert_eq!(tabs(&view, cx), (3, 2));

        // Each tab holds its own structure.
        view.update(cx, |this, cx| {
            let first = this.tabs[0].clone();
            first.update(cx, |tab, cx| tab.open_path("samples/caffeine.sdf".into(), cx));
        });
        cx.run_until_parked();

        click(cx, "tab-0");
        assert_eq!(tabs(&view, cx), (3, 0));
        cx.simulate_keystrokes("secondary-}");
        assert_eq!(tabs(&view, cx), (3, 1));
        cx.simulate_keystrokes("secondary-9");
        assert_eq!(tabs(&view, cx), (3, 2));
        cx.simulate_keystrokes("secondary-1");
        assert_eq!(tabs(&view, cx), (3, 0));

        click(cx, "close-tab-0");
        assert_eq!(tabs(&view, cx), (2, 0));
        cx.simulate_keystrokes("secondary-w");
        cx.simulate_keystrokes("secondary-w");
        // The last tab is emptied rather than closed.
        assert_eq!(tabs(&view, cx), (1, 0));
    }

    #[gpui::test]
    fn tabs_reorder_and_describe_the_session(cx: &mut TestAppContext) {
        let (view, cx) = open(cx);
        cx.simulate_keystrokes("secondary-t secondary-t");
        view.update(cx, |this, cx| {
            this.tabs[0].update(cx, |tab, cx| tab.reopen("sample:caffeine.sdf", cx));
        });
        cx.run_until_parked();
        let first = view.update(cx, |this, _| this.tabs[0].clone());

        // Drag the first tab onto the last.
        let from = cx.debug_bounds("tab-0").unwrap().center();
        let to = cx.debug_bounds("tab-2").unwrap().center();
        cx.simulate_mouse_down(from, gpui::MouseButton::Left, Modifiers::none());
        cx.simulate_mouse_move(from + gpui::point(px(12.), px(0.)), gpui::MouseButton::Left, Modifiers::none());
        cx.simulate_mouse_move(to, gpui::MouseButton::Left, Modifiers::none());
        cx.simulate_mouse_up(to, gpui::MouseButton::Left, Modifiers::none());
        cx.run_until_parked();
        let moved = view.update(cx, |this, _| this.tabs[2].clone());
        assert_eq!(moved, first, "the dragged tab did not move");

        let session = view.update(cx, |this, cx| this.session(cx));
        assert!(session.contains("sample:caffeine.sdf"), "{session}");
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
