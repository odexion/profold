//! profold — a small desktop front end for Apple's SimpleFold.

mod align;
mod dssp;
mod fetch;
mod fold;
mod job;
mod protein;
mod render;
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
        cx.bind_keys([gpui::KeyBinding::new("cmd-q", Quit, None)]);
        cx.on_action(|_: &Quit, cx: &mut App| cx.quit());
        cx.set_menus(vec![Menu {
            name: "profold".into(),
            items: vec![MenuItem::action("Quit", Quit)],
        }]);

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
                // `profold 1CRN` opens with that loaded; `profold --compare
                // 1CRN` fetches it and folds its sequence alongside.
                let args: Vec<String> = std::env::args().skip(1).collect();
                let compare = args.first().is_some_and(|a| a == "--compare");
                if let Some(target) = args.get(usize::from(compare)) {
                    let target = target.clone();
                    view.update(cx, |this, cx| {
                        if compare {
                            this.open_comparison(&target, cx)
                        } else {
                            this.open(&target, cx)
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
