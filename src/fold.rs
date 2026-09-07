//! Drives Apple's SimpleFold as a Python sidecar.
//!
//! SimpleFold is a Python package, so the app keeps a private virtualenv under
//! its application-support directory, installs `apple/ml-simplefold` into it
//! once, and then shells out to the `simplefold` CLI per prediction. Every step
//! streams its output back so the UI can show what is happening.

use anyhow::{anyhow, Context, Result};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::mpsc::Sender;

use crate::job::{self, Event};

pub const MODELS: [&str; 6] = [
    "simplefold_100M",
    "simplefold_360M",
    "simplefold_700M",
    "simplefold_1.1B",
    "simplefold_1.6B",
    "simplefold_3B",
];

pub fn model_label(name: &str) -> &str {
    name.trim_start_matches("simplefold_")
}

pub struct Engine {
    pub root: PathBuf,
}

fn which(bin: &str) -> Option<PathBuf> {
    // A GUI process does not inherit a login shell's PATH, so check the usual
    // install locations too.
    if let Ok(out) = Command::new("/usr/bin/which").arg(bin).output() {
        if out.status.success() {
            let p = PathBuf::from(String::from_utf8_lossy(&out.stdout).trim());
            if p.exists() {
                return Some(p);
            }
        }
    }
    let home = dirs::home_dir().unwrap_or_default();
    for dir in [
        home.join(".local/bin"),
        home.join(".cargo/bin"),
        PathBuf::from("/opt/homebrew/bin"),
        PathBuf::from("/usr/local/bin"),
        PathBuf::from("/usr/bin"),
    ] {
        let p = dir.join(bin);
        if p.exists() {
            return Some(p);
        }
    }
    None
}

impl Engine {
    /// Everything lives under `~/.profold` rather than Application Support:
    /// SimpleFold fetches its weights with `os.system(f"curl -L {url} -o
    /// {path}")`, unquoted, so any directory with a space in it silently
    /// breaks the download.
    pub fn new() -> Self {
        let root = dirs::home_dir()
            .unwrap_or_else(|| PathBuf::from("."))
            .join(".profold");
        Self { root }
    }

    pub fn repo(&self) -> PathBuf {
        self.root.join("ml-simplefold")
    }
    pub fn venv(&self) -> PathBuf {
        self.root.join("venv")
    }
    pub fn python(&self) -> PathBuf {
        self.venv().join("bin/python")
    }
    pub fn cli(&self) -> PathBuf {
        self.venv().join("bin/simplefold")
    }
    pub fn checkpoints(&self) -> PathBuf {
        self.root.join("checkpoints")
    }
    /// One stable output directory for every run. SimpleFold keeps its working
    /// cache — a 3.4 GB structure module and a 336 MB component dictionary —
    /// inside `output_dir/cache`, so a per-run directory would re-download
    /// both on every single fold.
    pub fn predictions(&self) -> PathBuf {
        self.root.join("predictions")
    }

    pub fn is_installed(&self) -> bool {
        self.cli().exists()
    }

    /// Clone the repository and build the virtualenv. Blocking; run off the
    /// main thread.
    pub fn install(&self, tx: &Sender<Event>) -> Result<()> {
        std::fs::create_dir_all(&self.root)?;

        let uv = which("uv").ok_or_else(|| {
            anyhow!("`uv` was not found. Install it with `brew install uv`, then try again.")
        })?;
        let git = which("git").ok_or_else(|| anyhow!("`git` was not found."))?;

        if !self.repo().join("pyproject.toml").exists() {
            let _ = std::fs::remove_dir_all(self.repo());
            stage(tx, "Cloning apple/ml-simplefold…");
            run(
                &git,
                &[
                    "clone",
                    "--depth",
                    "1",
                    "https://github.com/apple/ml-simplefold.git",
                    &self.repo().to_string_lossy(),
                ],
                &self.root,
                &[],
                tx,
            )?;
        }

        if !self.python().exists() {
            stage(tx, "Creating a Python 3.10 environment…");
            run(
                &uv,
                &[
                    "venv",
                    "--python",
                    "3.10",
                    &self.venv().to_string_lossy(),
                ],
                &self.root,
                &[],
                tx,
            )?;
        }

        let py = self.python().to_string_lossy().to_string();

        stage(tx, "Installing simplefold (this downloads a lot; be patient)…");
        run(
            &uv,
            &[
                "pip",
                "install",
                "--python",
                &py,
                "-e",
                &self.repo().to_string_lossy(),
            ],
            &self.root,
            &[],
            tx,
        )?;

        if cfg!(all(target_os = "macos", target_arch = "aarch64")) {
            stage(tx, "Installing the MLX backend for Apple silicon…");
            run(
                &uv,
                &["pip", "install", "--python", &py, "mlx==0.28.0"],
                &self.root,
                &[],
                tx,
            )?;
            stage(tx, "Installing ESM…");
            run(
                &uv,
                &[
                    "pip",
                    "install",
                    "--python",
                    &py,
                    "git+https://github.com/facebookresearch/esm.git",
                ],
                &self.root,
                &[],
                tx,
            )?;
        }

        if !self.is_installed() {
            return Err(anyhow!("installation finished but `simplefold` is missing"));
        }
        stage(tx, "Folding engine ready.");
        Ok(())
    }

    pub fn backend(&self) -> &'static str {
        if cfg!(all(target_os = "macos", target_arch = "aarch64")) {
            "mlx"
        } else {
            "torch"
        }
    }

    /// Has the checkpoint that pLDDT scoring needs already been fetched?
    pub fn has_plddt_weights(&self) -> bool {
        self.checkpoints().join("simplefold_1.6B.ckpt").exists()
    }

    /// Fold one sequence and return the path of the predicted structure.
    pub fn fold(
        &self,
        name: &str,
        sequence: &str,
        model: &str,
        steps: u32,
        plddt: bool,
        tx: &Sender<Event>,
    ) -> Result<PathBuf> {
        let run_dir = self.predictions();
        std::fs::create_dir_all(&run_dir)?;
        std::fs::create_dir_all(self.checkpoints())?;
        clear_previous(&run_dir)?;
        let started = std::time::SystemTime::now();

        let fasta = run_dir.join(format!("{name}.fasta"));
        let mut f = std::fs::File::create(&fasta)?;
        writeln!(f, ">{name}")?;
        for chunk in sequence.as_bytes().chunks(60) {
            writeln!(f, "{}", String::from_utf8_lossy(chunk))?;
        }
        drop(f);

        stage(
            tx,
            &format!(
                "Folding {} residues with {} on the {} backend…",
                sequence.len(),
                model_label(model),
                self.backend()
            ),
        );

        let steps = steps.to_string();
        let fasta_arg = fasta.to_string_lossy().to_string();
        let out_arg = run_dir.to_string_lossy().to_string();
        let mut args = vec![
            "--simplefold_model",
            model,
            "--num_steps",
            &steps,
            "--tau",
            "0.01",
            "--nsample_per_protein",
            "1",
            "--output_format",
            "pdb",
            "--backend",
            self.backend(),
            // Relative, and resolved against the working directory set below:
            // it must not contain a space, whatever the home directory is
            // called, because SimpleFold fetches weights through an unquoted
            // `os.system("curl ...")`.
            "--ckpt_dir",
            "checkpoints",
            "--fasta_path",
            &fasta_arg,
            "--output_dir",
            &out_arg,
        ];
        if plddt {
            args.push("--plddt");
        }

        run(
            &self.cli(),
            &args,
            &self.root,
            &[
                ("PYTHONUNBUFFERED", "1"),
                ("HF_HOME", &self.root.join("hf").to_string_lossy()),
            ],
            tx,
        )?;

        // Name the file we expect; fall back to whatever appeared, in case a
        // future version numbers its output differently.
        let expected = run_dir
            .join(format!("predictions_{model}"))
            .join(format!("{name}_sampled_0.pdb"));
        if expected.exists() {
            return Ok(expected);
        }
        newest_structure(&run_dir, started)
            .ok_or_else(|| anyhow!("SimpleFold produced no structure file"))
    }
}

/// SimpleFold predicts *every* record left in its output directory, not just
/// the FASTA it was pointed at, so a directory reused across runs quietly folds
/// stale inputs and can hand back the wrong protein. Clear the previous run out
/// while keeping the multi-gigabyte cache that lives in the same place.
fn clear_previous(dir: &Path) -> Result<()> {
    if dir.file_name().is_none_or(|n| n != "predictions") {
        return Err(anyhow!("refusing to clear {}", dir.display()));
    }
    for entry in std::fs::read_dir(dir)?.flatten() {
        if entry.file_name() == "cache" {
            continue;
        }
        let path = entry.path();
        let result = if path.is_dir() {
            std::fs::remove_dir_all(&path)
        } else {
            std::fs::remove_file(&path)
        };
        result.with_context(|| format!("clearing {}", path.display()))?;
    }
    Ok(())
}

/// The newest structure written since `since`, so that a failed run cannot
/// hand back the previous run's output.
fn newest_structure(dir: &Path, since: std::time::SystemTime) -> Option<PathBuf> {
    let mut best: Option<(std::time::SystemTime, PathBuf)> = None;
    let mut stack = vec![dir.to_path_buf()];
    while let Some(d) = stack.pop() {
        for entry in std::fs::read_dir(&d).ok()?.flatten() {
            let p = entry.path();
            if p.is_dir() {
                stack.push(p);
                continue;
            }
            let ext = p.extension().and_then(|e| e.to_str()).unwrap_or("");
            if ext != "pdb" && ext != "cif" {
                continue;
            }
            let t = entry
                .metadata()
                .and_then(|m| m.modified())
                .unwrap_or(std::time::UNIX_EPOCH);
            if t < since {
                continue;
            }
            if best.as_ref().is_none_or(|(bt, _)| t >= *bt) {
                best = Some((t, p));
            }
        }
    }
    best.map(|(_, p)| p)
}

/// Announce a phase in our own words, rather than the subprocess's.
fn stage(tx: &Sender<Event>, msg: &str) {
    let _ = tx.send(Event::Stage(msg.to_string()));
}

/// Run a command, streaming stdout and stderr as events.
fn run(
    program: &Path,
    args: &[&str],
    cwd: &Path,
    env: &[(&str, &str)],
    tx: &Sender<Event>,
) -> Result<()> {
    let mut cmd = Command::new(program);
    cmd.args(args)
        .current_dir(cwd)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .stdin(Stdio::null());
    for (k, v) in env {
        cmd.env(k, v);
    }

    let mut child = cmd
        .spawn()
        .with_context(|| format!("running {}", program.display()))?;

    let out = child.stdout.take().unwrap();
    let err = child.stderr.take().unwrap();
    let tx_out = tx.clone();
    let tx_err = tx.clone();
    let h1 = std::thread::spawn(move || pump(out, tx_out));
    let h2 = std::thread::spawn(move || pump(err, tx_err));

    let status = child.wait()?;
    let _ = h1.join();
    let _ = h2.join();

    if !status.success() {
        return Err(anyhow!(
            "{} exited with {}",
            program.file_name().unwrap_or_default().to_string_lossy(),
            status.code().unwrap_or(-1)
        ));
    }
    Ok(())
}

/// Split on both newlines and carriage returns so that progress bars, which
/// rewrite a single line, arrive as `Progress` rather than flooding the log.
fn pump<R: Read + Send + 'static>(r: R, tx: Sender<Event>) {
    let mut reader = BufReader::new(r);
    let mut chunk: Vec<u8> = Vec::new();
    loop {
        chunk.clear();
        match reader.read_until(b'\n', &mut chunk) {
            Ok(0) | Err(_) => break,
            Ok(_) => {}
        }
        let text = String::from_utf8_lossy(&chunk);
        // A progress bar rewrites one line with carriage returns; only the last
        // segment of such a chunk is current, and none of it belongs in the log.
        let rewritten = text.trim_end_matches('\n').contains('\r');
        for part in text.split('\r') {
            let line = part.trim_end();
            if line.trim().is_empty() {
                continue;
            }
            let event = if rewritten {
                Event::Progress(job::parse(line))
            } else {
                Event::Log(line.to_string())
            };
            if tx.send(event).is_err() {
                return;
            }
        }
    }
}
