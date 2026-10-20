# AGENTS.md

Notes for coding agents working on profold.

## Rules

1. **Format before committing.** Run `cargo fmt` before every commit.
2. **Cross-platform by default.** Every feature must work on macOS, Linux and Windows (arm64 and x86-64). Don't use platform-only APIs without a fallback for the other platforms, and keep GPU work portable (wgpu, and GLES where possible).
3. **Test it yourself first.** Build the app, run it, and check the feature and UI with screenshots or computer use before showing it to the user. Fix what you find, then hand it over.
4. **Commit on `main`.** No feature branches unless asked. Commit only when asked, and don't push unless asked.
5. **Accuracy comes first.** profold is a scientific tool, so what it shows and computes must be correct: coordinates, units, bonds, periodic boundaries, per-atom properties, analysis results. Never trade correctness for speed or looks. Check results against known values or another tool (e.g. OVITO, VMD), and add tests for anything numerical.
6. **Make hot paths fast, without cheating.** Use parallelism (rayon), threading, the GPU (wgpu compute or render passes), better algorithms and data structures (BVHs, spatial grids, caching), and less copying and allocation. Don't buy speed by quietly dropping atoms or frames, lowering precision, skipping work the result depends on, approximating without telling the user, or special-casing benchmarks. If an approximation is truly worth it, make it opt-in and visible, and ask the user first. Measure before and after, because small per-frame jobs can get slower when split across threads.
7. **Keep changes surgical.** Touch only what the task needs, with no drive-by refactors or rewrites of unrelated code. Don't break existing UI or UX: layouts, menus, shortcuts, settings and workflows should keep working as before unless the user asks for a change. If you find a bug outside the task, don't fix it on the spot: tell the user what it is and where, and ask before fixing it.

## Workflow

- **After every change:** `cargo build`, then launch `./target/debug/profold` (with a sample from `samples/` if it fits) and leave it open for the user. Record the PID you start, and only kill your own instances, never the user's windows.
- **Clean settings:** `HOME=<scratch> ./target/debug/profold ...` starts with default settings and leaves `~/.profold` alone.
- **Driving the app:** each CLI argument opens a tab (`profold caffeine 4HHB traj.lammpstrj`), `--fold` opens a fold tab, and `PROFOLD_TYPE="text"` types into the search field. Headless GPUI tests in `src/ui/` cover clicks, keys and drags.
- **Screenshots:** capture only profold's window (`screencapture -l <window id>`). An inactive window stops repainting, so bring it to the front first.

## Performance

- Rendering is wgpu (`src/gpu.rs`, `gpu.wgsl`), with a CPU fallback in `render.rs` (`PROFOLD_CPU=1`).
- Measure in the app with `PROFOLD_SPIN=1 ./target/release/profold <file>`. Use `bench_gpu` for A/B tests and report before/after numbers.
- `gpu_matches_cpu` must stay green, and the CPU fallback must stay bit-identical.
- `data.denge` (LAMMPS, 42k atoms) is the benchmark file.

## Releases

- Bump the version in `Cargo.toml` and push a `vX.Y.Z` tag. CI builds 6 targets and attaches `install.sh` and `install.ps1`.
- Installs use a uv-style curl installer (`curl … | sh`, per-user location, no Apple Developer account, ad-hoc signing). The app updates itself, checking at most once a day.
- To test the installer and updater locally, serve archives with `python3 -m http.server` and set `PROFOLD_DOWNLOAD_BASE` / `PROFOLD_RELEASES_URL`.

## Roadmap

`docs/ROADMAP.md` ranks the features missing compared with OVITO and VMD. Work down it in order and update it as items land.
