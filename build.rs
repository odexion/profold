//! Puts the icon into the Windows executable, for Explorer and the taskbar.
//! Elsewhere there is nothing to do: macOS takes it from the bundle and Linux
//! from the desktop entry the installer writes.

fn main() {
    println!("cargo:rerun-if-changed=assets/profold.ico");
    // Release builds name their target for the updater (see `update.rs`).
    println!("cargo:rerun-if-env-changed=PROFOLD_TARGET");

    #[cfg(windows)]
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("windows") {
        let mut res = winresource::WindowsResource::new();
        res.set_icon("assets/profold.ico");
        res.compile().expect("embedding the icon");
    }
}
