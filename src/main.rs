//! profold — a molecular viewer for proteins, molecules, crystals and
//! simulation trajectories, with structure prediction as a plugin.

mod align;
mod bonds;
mod dssp;
mod elements;
mod fetch;
mod formats;
mod job;
mod plugins;
mod protein;
mod render;
mod samples;
mod settings;
mod theme;
mod ui;

use gpui::{
    actions, AppContext, point, px, size, App, Application, Bounds, Menu, MenuItem, TitlebarOptions,
    WindowBounds, WindowOptions,
};

actions!(profold, [Quit]);

fn main() {
    Application::new().run(|cx: &mut App| {
        ui::input::bind_keys(cx);
        cx.bind_keys([
            gpui::KeyBinding::new("cmd-q", Quit, None),
            gpui::KeyBinding::new("cmd-o", ui::app::OpenFile, None),
        ]);
        cx.on_action(|_: &Quit, cx: &mut App| cx.quit());
        cx.set_menus(vec![
            Menu {
                name: "profold".into(),
                items: vec![MenuItem::action("Quit", Quit)],
            },
            Menu {
                name: "File".into(),
                items: vec![MenuItem::action("Open…", ui::app::OpenFile)],
            },
        ]);

        let bounds = Bounds::centered(None, size(px(1180.), px(760.)), cx);
        cx.open_window(
            WindowOptions {
                window_bounds: Some(WindowBounds::Windowed(bounds)),
                titlebar: Some(TitlebarOptions {
                    title: Some("profold".into()),
                    appears_transparent: true,
                    traffic_light_position: Some(point(px(16.), px(15.))),
                }),
                focus: true,
                ..Default::default()
            },
            |window, cx| {
                let view = cx.new(|cx| ui::app::Profold::new(window, cx));
                // `profold 1CRN` or `profold traj.xyz` opens with that
                // loaded. With the SimpleFold plugin, `--fold` predicts an
                // entry instead of fetching it, and `--compare` does both and
                // lays one over the other.
                let args: Vec<String> = std::env::args().skip(1).collect();
                let verb = args.first().filter(|a| a.starts_with("--")).cloned();
                if let Some(target) = args.get(usize::from(verb.is_some())) {
                    let target = target.clone();
                    view.update(cx, |this, cx| match verb.as_deref() {
                        #[cfg(feature = "simplefold")]
                        Some("--compare") => this.open_comparison(&target, cx),
                        #[cfg(feature = "simplefold")]
                        Some("--fold") => this.open_fold(&target, cx),
                        _ => {
                            let path = std::path::PathBuf::from(&target);
                            if path.is_file() {
                                this.open_path(path, cx);
                            } else {
                                this.open(&target, cx);
                            }
                        }
                    });
                }
                view
            },
        )
        .unwrap();
        cx.activate(true);
    });
}
