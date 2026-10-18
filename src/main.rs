//! profold — a molecular viewer for proteins, molecules, crystals and
//! simulation trajectories, with structure prediction as a plugin.

// A release build on Windows is a window, not a console program.
#![cfg_attr(all(windows, not(debug_assertions)), windows_subsystem = "windows")]

mod align;
mod bonds;
mod dssp;
mod elements;
mod fetch;
mod formats;
mod gpu;
mod job;
mod measure;
mod plugins;
mod protein;
mod render;
mod samples;
mod search;
mod settings;
mod surface;
mod theme;
mod ui;
mod update;

use gpui::{
    AppContext, point, px, size, App, Application, Bounds, TitlebarOptions, WindowBounds,
    WindowOptions,
};

fn main() {
    if std::env::args().nth(1).is_some_and(|a| a == "--version" || a == "-V") {
        #[cfg(windows)]
        unsafe {
            windows_sys::Win32::System::Console::AttachConsole(
                windows_sys::Win32::System::Console::ATTACH_PARENT_PROCESS,
            );
        }
        println!("profold {}", update::VERSION);
        return;
    }
    update::tidy();
    // Setting up the graphics card takes a moment; do it while the window
    // opens rather than on the first frame.
    std::thread::spawn(|| gpu::get());
    Application::new().run(|cx: &mut App| {
        #[cfg(target_os = "macos")]
        set_dock_icon();
        ui::input::bind_keys(cx);
        ui::workspace::bind_keys(cx);
        // Closing the window quits, except on a Mac, where the app stays in
        // the Dock as Mac apps do.
        #[cfg(not(target_os = "macos"))]
        cx.on_window_closed(|cx| {
            if cx.windows().is_empty() {
                cx.quit();
            }
        })
        .detach();
        cx.set_menus(ui::workspace::menus());

        let bounds = Bounds::centered(None, size(px(1180.), px(760.)), cx);
        cx.open_window(
            WindowOptions {
                window_bounds: Some(WindowBounds::Windowed(bounds)),
                titlebar: Some(TitlebarOptions {
                    title: Some("profold".into()),
                    appears_transparent: true,
                    traffic_light_position: Some(point(px(16.), px(15.))),
                }),
                // Linux asks the compositor for its own title bar; where it
                // has none to give (GNOME on Wayland), the tab bar draws the
                // window's buttons itself.
                #[cfg(target_os = "linux")]
                window_decorations: Some(gpui::WindowDecorations::Server),
                app_id: Some("dev.profold".into()),
                focus: true,
                ..Default::default()
            },
            |window, cx| {
                let view = cx.new(|cx| ui::workspace::Workspace::new(window, cx));
                // `profold 1CRN` or `profold traj.xyz` opens with that
                // loaded, and each further argument in a tab of its own. With
                // the SimpleFold plugin, `--fold` predicts an entry in a fold
                // tab instead of fetching it, and `--compare` does both and
                // lays one over the other.
                let args: Vec<String> = std::env::args().skip(1).collect();
                let verb = args.first().filter(|a| a.starts_with("--")).cloned();
                let targets = &args[usize::from(verb.is_some())..];
                view.update(cx, |this, cx| {
                    // `--fold` alone opens an empty fold tab.
                    #[cfg(feature = "simplefold")]
                    if verb.as_deref() == Some("--fold") && targets.is_empty() {
                        this.add(ui::app::Mode::Fold, window, cx);
                    }
                    // With nothing named, the tabs of last time come back.
                    if targets.is_empty() && verb.is_none() {
                        this.restore(window, cx);
                    }
                    for (n, target) in targets.iter().enumerate() {
                        let path = std::path::PathBuf::from(target);
                        match verb.as_deref() {
                            #[cfg(feature = "simplefold")]
                            Some("--compare") | Some("--fold") => {
                                let tab = this.add(ui::app::Mode::Fold, window, cx);
                                let compare = verb.as_deref() == Some("--compare");
                                tab.update(cx, |tab, cx| {
                                    if compare {
                                        tab.open_comparison(target, cx)
                                    } else {
                                        tab.open_fold(target, cx)
                                    }
                                });
                            }
                            _ => {
                                let tab = if n == 0 {
                                    this.active()
                                } else {
                                    this.add(ui::app::Mode::View, window, cx)
                                };
                                tab.update(cx, |tab, cx| {
                                    if path.is_file() {
                                        tab.open_path(path, cx)
                                    } else {
                                        tab.open(target, cx)
                                    }
                                });
                            }
                        }
                    }
                    // `PROFOLD_TYPE=caffeine` types into the field without
                    // pressing Return, for looking at the search results.
                    if let Ok(text) = std::env::var("PROFOLD_TYPE") {
                        this.active().update(cx, |tab, cx| tab.type_text(&text, cx));
                    }
                });
                view
            },
        )
        .unwrap();
        cx.activate(true);
    });
}

/// Puts the logo in the Dock. A bundled build gets it from Info.plist; this
/// covers `cargo run` and a bare binary.
#[cfg(target_os = "macos")]
fn set_dock_icon() {
    use objc2::{AllocAnyThread, MainThreadMarker};
    use objc2_app_kit::{NSApplication, NSImage};
    use objc2_foundation::NSData;

    let Some(mtm) = MainThreadMarker::new() else { return };
    let data = NSData::with_bytes(include_bytes!("../assets/icon.png"));
    if let Some(image) = NSImage::initWithData(NSImage::alloc(), &data) {
        unsafe { NSApplication::sharedApplication(mtm).setApplicationIconImage(Some(&image)) };
    }
}
