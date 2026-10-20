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
- Ray tracing (`src/trace.rs`) runs on the card through a compute shader (`trace.wgsl`, driven from `gpu_trace.rs`), and on the CPU where the card can't. The shader follows the CPU loop step for step. `gpu_traces_what_the_cpu_traces` must stay green, and `PROFOLD_TRACE_CPU=1` keeps tracing on the CPU for A/B tests. Time it with `cargo test --release trace_bench -- --ignored --nocapture`, with `PROFOLD_TRACE_OUT=<dir>` to write the pictures out for comparing.
- `data.denge` (LAMMPS, 42k atoms) is the benchmark file. `dump.lammpstrj` (43k atoms × 801 frames, 1.2 GB, gitignored) is the large trajectory.
- Measure the other operations with `cargo test --release ops_bench -- --ignored --nocapture`, and one operation on the large dump with `PROFOLD_OP=<op> cargo test --release dump_bench -- --ignored --nocapture` (`src/bench.rs`). Run that under `/usr/bin/time -l` to get its peak memory.

### How to optimize

1. **Measure first, at real sizes.** Time each operation on the benchmark files, and on the largest input it will realistically meet (long trajectories, big boxes, large grids), not only on samples. Record time and peak memory before changing anything. Take the best of a few runs. Report a before/after table.
2. **Read the code for how cost scales, not just how fast it runs now.** Something fast today can still hold memory that grows with frames × atoms, or redo work on every repaint. Find what grows with trajectory length, atom count or frame rate, and fix that even when today's numbers look fine.
3. **Keep memory flat in the long dimension.** Stream over frames instead of collecting them. Reduce as you go: running sums, Welford's mean and variance merged between threads by Chan's method, min/max, histograms. Borrow data where it lives (`Cow`, slices) instead of cloning a frame per step.
4. **Read big files in pieces.** Read large inputs in chunks cut at record boundaries, parse each chunk's records in parallel, and drop the text before reading the next. Read the next chunk on another thread while the current one parses, holding at most one chunk ahead.
5. **Find neighbours with flat cell lists.** For distance searches (contacts, SASA, bonds, `within`), sort atoms once into a dense grid of cubes and scan contiguous runs of neighbouring cells. Avoid hash maps of `Vec`s and per-atom allocation. Cap the number of cells relative to the atom count, so a sparse box doesn't blow up memory.
6. **Do less of the same work.** Cache what the UI asks for on every repaint, keyed on everything it depends on (structure, frames, selection, settings). Order early-exit checks so the likely answer comes first (nearest neighbour first; try first whichever neighbour last gave the answer). Do costly work only when the user asks for it, show its progress, and let it be cancelled.
7. **Parallelize the big loop, once.** Use rayon over frames or atoms, whichever is large. Don't nest small parallel jobs. A per-frame job that's too small gets slower when split, so measure it.
8. **Keep results exact, and prove it.** An optimization must give the same answer as the plain method. That means bit-for-bit where the arithmetic is unchanged, and within float rounding where the order of sums changes. For each one, add a test against a slow reference written out in full: brute force over all pairs, the two-pass formula, or the whole input read in one piece. Use edge cases: tiny chunks that cut records, a single frame, empty selections. Never trade accuracy for speed silently. Fewer sample points, fewer frames or lower precision is a choice to offer the user, not to make for them.

## Releases

- Bump the version in `Cargo.toml` and push a `vX.Y.Z` tag. CI builds 6 targets and attaches `install.sh` and `install.ps1`.
- Installs use a uv-style curl installer (`curl … | sh`, per-user location, no Apple Developer account, ad-hoc signing). The app updates itself, checking at most once a day.
- To test the installer and updater locally, serve archives with `python3 -m http.server` and set `PROFOLD_DOWNLOAD_BASE` / `PROFOLD_RELEASES_URL`.

## Roadmap

`docs/ROADMAP.md` ranks the features missing compared with OVITO and VMD. Work down it in order and update it as items land.
