//! Drives Apple's SimpleFold as a Python sidecar.
//!
//! SimpleFold is a Python package, so the app keeps a private virtualenv under
//! its application-support directory, installs `apple/ml-simplefold` into it
//! once, and then shells out to the `simplefold` CLI per prediction. Every step
//! streams its output back so the UI can show what is happening.

use anyhow::{anyhow, Context, Result};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::mpsc::Sender;

use super::download;
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

/// The Python the folding environment is built on. SimpleFold asks for 3.10
/// or later; 3.11 is the oldest whose SciPy wheels load on current macOS.
const PYTHON: &str = "3.11";

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

    /// The files a fold with `model` needs, and where SimpleFold looks for
    /// them: its checkpoints, ESM2 in torch's hub, and Boltz's component
    /// dictionary and structure module in the working cache.
    fn weights(&self, model: &str, plddt: bool) -> Vec<download::Wanted> {
        const APPLE: &str = "https://ml-site.cdn-apple.com/models/simplefold";
        const ESM: &str = "https://dl.fbaipublicfiles.com/fair-esm";
        const BOLTZ: &str = "https://huggingface.co/boltz-community/boltz-1/resolve/main";
        let want = |url: String, path: PathBuf, what: &str| download::Wanted {
            url,
            path,
            what: what.to_string(),
        };
        let ckpt = self.checkpoints();
        let hub = torch_hub().join("checkpoints");
        let cache = self.predictions().join("cache");
        // Torch's own downloads leave `.partial` files when interrupted,
        // which it never resumes; they are dead weight.
        if let Ok(entries) = std::fs::read_dir(&hub) {
            for entry in entries.flatten() {
                if entry.file_name().to_string_lossy().ends_with(".partial") {
                    let _ = std::fs::remove_file(entry.path());
                }
            }
        }

        let mut files = vec![want(
            format!("{APPLE}/{model}.ckpt"),
            ckpt.join(format!("{model}.ckpt")),
            &format!("SimpleFold {} weights", model_label(model)),
        )];
        if plddt {
            files.push(want(
                format!("{APPLE}/plddt_module_1.6B.ckpt"),
                ckpt.join("plddt.ckpt"),
                "the pLDDT module",
            ));
            files.push(want(
                format!("{APPLE}/simplefold_1.6B.ckpt"),
                ckpt.join("simplefold_1.6B.ckpt"),
                "SimpleFold 1.6B weights for pLDDT",
            ));
        }
        files.extend([
            want(
                format!("{ESM}/models/esm2_t36_3B_UR50D.pt"),
                hub.join("esm2_t36_3B_UR50D.pt"),
                "the ESM2-3B language model",
            ),
            want(
                format!("{ESM}/regression/esm2_t36_3B_UR50D-contact-regression.pt"),
                hub.join("esm2_t36_3B_UR50D-contact-regression.pt"),
                "ESM2-3B's contact regression",
            ),
            want(format!("{BOLTZ}/ccd.pkl"), cache.join("ccd.pkl"), "the chemical component dictionary"),
            want(
                format!("{BOLTZ}/boltz1_conf.ckpt"),
                cache.join("boltz1_conf.ckpt"),
                "the Boltz-1 structure module",
            ),
        ]);
        files
    }

    pub fn is_installed(&self) -> bool {
        self.cli().exists() && self.venv_is_current() && self.installed_marker().exists()
    }

    /// Written when an install has run to the end, so that one cut short,
    /// with the `simplefold` script in place but packages missing, is
    /// finished next time rather than trusted.
    fn installed_marker(&self) -> PathBuf {
        self.venv().join(".profold-installed")
    }

    /// Whether the environment runs the Python this version of profold sets
    /// up. Early versions used 3.10, whose last SciPy no longer loads on
    /// recent macOS (dyld rejects how its extensions lay out thread-local
    /// storage), so a 3.10 environment is rebuilt.
    fn venv_is_current(&self) -> bool {
        let Ok(cfg) = std::fs::read_to_string(self.venv().join("pyvenv.cfg")) else {
            return false;
        };
        cfg.lines()
            .filter_map(|l| l.split_once('='))
            .filter(|(k, _)| matches!(k.trim(), "version_info" | "version"))
            .map(|(_, v)| v.trim())
            .any(|v| v == PYTHON || v.starts_with(&format!("{PYTHON}.")))
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

        if self.python().exists() && !self.venv_is_current() {
            stage(tx, &format!("Rebuilding the Python environment on {PYTHON}…"));
            std::fs::remove_dir_all(self.venv())
                .with_context(|| format!("removing {}", self.venv().display()))?;
        }
        if !self.python().exists() {
            stage(tx, &format!("Creating a Python {PYTHON} environment…"));
            run(
                &uv,
                &[
                    "venv",
                    "--python",
                    PYTHON,
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

        if !self.cli().exists() {
            return Err(anyhow!("installation finished but `simplefold` is missing"));
        }
        std::fs::write(self.installed_marker(), PYTHON)?;
        self.ease_esm_load(tx);
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

    /// Is this model's checkpoint already on disk? SimpleFold fetches a
    /// missing one on first use, which is a wait worth knowing about before
    /// picking it.
    pub fn has_checkpoint(&self, model: &str) -> bool {
        self.checkpoints().join(format!("{model}.ckpt")).exists()
    }

    /// Has the checkpoint that pLDDT scoring needs already been fetched?
    pub fn has_plddt_weights(&self) -> bool {
        self.checkpoints().join("simplefold_1.6B.ckpt").exists()
    }

    /// SimpleFold's MLX path converts ESM2-3B in a single comprehension, so a
    /// torch copy and an MLX copy of the same eleven gigabytes are resident at
    /// once. That peak is what folding costs, whatever the sequence: the
    /// language model is loaded in full for a forty-residue chain. The
    /// conversion is rewritten here to free each tensor as it goes, which
    /// halves the peak and produces exactly the same weights.
    ///
    /// It is their file, so the edit is made only when the code is the code
    /// profold knows, and only once.
    fn ease_esm_load(&self, tx: &Sender<Event>) {
        const MARKER: &str = "# profold: one tensor at a time";
        const ORIGINAL: &str = concat!(
            "def initialize_esm_model(args, device):\n",
            "    # load ESM2 model\n",
            "    esm_model, esm_dict = esm_registry[\"esm2_3B\"]()\n",
            "    af2_to_esm = _af2_to_esm(esm_dict)\n",
            "\n",
            "    if args.backend == 'torch':\n",
            "        esm_model = esm_model.to(device)\n",
            "        af2_to_esm = af2_to_esm.to(device)\n",
            "    elif args.backend == 'mlx':\n",
            "        esm_model_mlx = ESM2MLX(num_layers=36, embed_dim=2560, attention_heads=40)\n",
            "        esm_state_dict_torch = esm_model.cpu().state_dict()\n",
            "\n",
            "        esm_state_dict_torch = {k: mx.array(v) for k, v in starmap(map_torch_to_mlx, esm_state_dict_torch.items()) if k is not None}\n",
            "        esm_model_mlx.update(tree_unflatten(list(esm_state_dict_torch.items())))\n",
            "        esm_model = esm_model_mlx\n",
            "    print(f\"pLM ESM-3B loaded with {args.backend} backend.\")\n",
        );
        const PATCHED: &str = concat!(
            "def initialize_esm_model(args, device):\n",
            "    # profold: one tensor at a time. Upstream builds a float32 torch copy of\n",
            "    # ESM2-3B and converts that to mlx, so two eleven-gigabyte copies are\n",
            "    # resident at once whatever is being folded. On the mlx backend the\n",
            "    # weights are read straight from the checkpoint instead, in short runs so\n",
            "    # the mapped pages do not pile up beside the model being built. Same\n",
            "    # weights, same precision, about half the peak.\n",
            "    if args.backend == 'mlx':\n",
            "        import re as _re\n",
            "        import esm as _esm\n",
            "\n",
            "        name = \"esm2_t36_3B_UR50D\"\n",
            "        url = f\"https://dl.fbaipublicfiles.com/fair-esm/models/{name}.pt\"\n",
            "        ckpt = os.path.join(torch.hub.get_dir(), \"checkpoints\", f\"{name}.pt\")\n",
            "        regression = ckpt.replace(\".pt\", \"-contact-regression.pt\")\n",
            "        for path, source in ((ckpt, url), (regression, url.replace(\"/models/\", \"/regression/\").replace(\".pt\", \"-contact-regression.pt\"))):\n",
            "            if not os.path.exists(path):\n",
            "                os.makedirs(os.path.dirname(path), exist_ok=True)\n",
            "                torch.hub.download_url_to_file(source, path)\n",
            "\n",
            "        esm_dict = _esm.data.Alphabet.from_architecture(\"ESM-1b\")\n",
            "        af2_to_esm = _af2_to_esm(esm_dict)\n",
            "        prefix = _re.compile(\"^encoder.sentence_encoder.|^encoder.\")\n",
            "\n",
            "        def read(names):\n",
            "            # The mapping is dropped with this frame, freeing what it touched.\n",
            "            weights = torch.load(ckpt, map_location=\"cpu\", mmap=True, weights_only=False)[\"model\"]\n",
            "            done = {}\n",
            "            for stored in names:\n",
            "                key, value = map_torch_to_mlx(prefix.sub(\"\", stored), weights[stored].float())\n",
            "                if key is not None:\n",
            "                    done[key] = mx.array(value)\n",
            "            return done\n",
            "\n",
            "        stored_names = list(\n",
            "            torch.load(ckpt, map_location=\"cpu\", mmap=True, weights_only=False)[\"model\"].keys()\n",
            "        )\n",
            "        converted = {}\n",
            "        for start in range(0, len(stored_names), 48):\n",
            "            converted.update(read(stored_names[start:start + 48]))\n",
            "        for stored, value in torch.load(regression, map_location=\"cpu\", weights_only=False)[\"model\"].items():\n",
            "            key, value = map_torch_to_mlx(prefix.sub(\"\", stored), value.float())\n",
            "            if key is not None:\n",
            "                converted[key] = mx.array(value)\n",
            "\n",
            "        esm_model = ESM2MLX(num_layers=36, embed_dim=2560, attention_heads=40)\n",
            "        esm_model.update(tree_unflatten(list(converted.items())))\n",
            "        del converted\n",
            "\n",
            "        print(\"pLM ESM-3B loaded with mlx backend.\")\n",
            "        esm_model.eval()\n",
            "        return esm_model, esm_dict, af2_to_esm\n",
            "\n",
            "    # load ESM2 model\n",
            "    esm_model, esm_dict = esm_registry[\"esm2_3B\"]()\n",
            "    af2_to_esm = _af2_to_esm(esm_dict)\n",
            "\n",
            "    if args.backend == 'torch':\n",
            "        esm_model = esm_model.to(device)\n",
            "        af2_to_esm = af2_to_esm.to(device)\n",
            "    print(f\"pLM ESM-3B loaded with {args.backend} backend.\")\n",
        );

        // The folding checkpoint is converted the same way, which is eleven
        // gigabytes of its own at the 3B model.
        const ORIGINAL_MODEL: &str = concat!(
            "        mlx_state_dict = {k: mx.array(v) for k, v in ",
            "starmap(map_torch_to_mlx, checkpoint.items()) if k is not None}\n",
            "        model.update(tree_unflatten(list(mlx_state_dict.items())))\n",
        );
        const PATCHED_MODEL: &str = concat!(
            "        mlx_state_dict = {}\n",
            "        for name in list(checkpoint.keys()):\n",
            "            key, value = map_torch_to_mlx(name, checkpoint.pop(name))\n",
            "            if key is not None:\n",
            "                mlx_state_dict[key] = mx.array(value)\n",
            "            del value\n",
            "        model.update(tree_unflatten(list(mlx_state_dict.items())))\n",
            "        del mlx_state_dict\n",
        );

        let path = self.repo().join("src/simplefold/inference.py");
        let Ok(source) = std::fs::read_to_string(&path) else {
            return;
        };
        if source.contains(MARKER) {
            return;
        }
        if !source.contains(ORIGINAL) {
            let _ = tx.send(Event::Log(
                "SimpleFold loads ESM differently than expected; leaving its code alone.".into(),
            ));
            return;
        }
        let eased = source
            .replace(ORIGINAL, PATCHED)
            .replace(ORIGINAL_MODEL, PATCHED_MODEL);
        match std::fs::write(&path, eased) {
            Ok(()) => {
                let _ = tx.send(Event::Log(
                    "Eased SimpleFold's model loading to one tensor at a time.".into(),
                ));
            }
            Err(e) => {
                let _ = tx.send(Event::Log(format!("Could not ease the model loading: {e}")));
            }
        }
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
        self.ease_esm_load(tx);
        // Everything SimpleFold would fetch for itself, fetched first and
        // whole; see download.rs.
        download::ensure(&self.weights(model, plddt), &self.root.join("downloads.json"), tx)?;
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
/// Where torch keeps what it downloads: `$TORCH_HOME/hub`, or the hub under
/// the cache directory, as `torch.hub.get_dir` works it out.
fn torch_hub() -> PathBuf {
    if let Some(home) = std::env::var_os("TORCH_HOME") {
        return PathBuf::from(home).join("hub");
    }
    let cache = std::env::var_os("XDG_CACHE_HOME")
        .map(PathBuf::from)
        .or_else(|| dirs::home_dir().map(|h| h.join(".cache")))
        .unwrap_or_else(|| PathBuf::from(".cache"));
    cache.join("torch").join("hub")
}

fn stage(tx: &Sender<Event>, msg: &str) {
    STAGE.with_borrow_mut(|s| *s = msg.to_string());
    let _ = tx.send(Event::Stage(msg.to_string()));
}

thread_local! {
    /// The stage the job on this thread last announced, to go back to once a
    /// download inside a step is over.
    static STAGE: std::cell::RefCell<String> = Default::default();
}

/// The total size, when a progress line is a download rather than work:
/// curl's transfer table, or a tqdm bar counting bytes ("31M/5.29G").
fn download_size(line: &str) -> Option<String> {
    let f: Vec<&str> = line.split_whitespace().collect();
    if f.len() >= 10
        && f[f.len() - 2].matches(':').count() == 2
        && f[f.len() - 3].matches(':').count() == 2
        && f[0].parse::<f32>().is_ok()
    {
        return Some(f[1].to_string());
    }
    let bytes = |s: &str| {
        let s = s.trim_end_matches('B');
        s.ends_with(['k', 'K', 'M', 'G', 'T']) && s[..s.len() - 1].parse::<f32>().is_ok()
    };
    let counts = line.rsplit('|').next()?.split_whitespace().next()?;
    let (done, total) = counts.split_once('/')?;
    (bytes(done) && bytes(total)).then(|| total.trim_end_matches('B').to_string() + "B")
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
    let label = STAGE.with_borrow(|s| s.clone());
    let label_err = label.clone();
    let h1 = std::thread::spawn(move || pump(out, tx_out, label));
    let h2 = std::thread::spawn(move || pump(err, tx_err, label_err));

    // Wait, asking now and then whether anyone is still listening: when the
    // tab that started the job is closed, its end of the channel goes, and
    // a fold nobody will see is stopped rather than left to run for minutes.
    let status = loop {
        if let Some(status) = child.try_wait()? {
            break status;
        }
        if tx.send(Event::Heartbeat).is_err() {
            let _ = child.kill();
            let _ = child.wait();
            return Err(anyhow!("cancelled"));
        }
        std::thread::sleep(std::time::Duration::from_millis(250));
    };
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
///
/// A step can stop to download what it needs (SimpleFold fetches ESM2's
/// weights, 5.3 GB, the first time it folds), and saying "Folding" through
/// hours of that would be untrue. So a bar that counts bytes changes the
/// stage to a download, and the first bar after it that does not puts back
/// `label`, what the step was announced as.
fn pump<R: Read + Send + 'static>(mut r: R, tx: Sender<Event>, label: String) {
    let mut downloading = false;
    // A bar rewrites its line with carriage returns and may never end one
    // with a newline until it finishes, so the output is cut at either as it
    // arrives rather than read a line at a time.
    let mut line: Vec<u8> = Vec::new();
    let mut buf = [0u8; 8192];
    // Whether the line under way has been rewritten, which makes it a bar.
    let mut rewritten = false;
    loop {
        let n = match r.read(&mut buf) {
            Ok(0) | Err(_) => break,
            Ok(n) => n,
        };
        for &b in &buf[..n] {
            if b != b'\r' && b != b'\n' {
                line.push(b);
                continue;
            }
            let bar = rewritten || b == b'\r';
            let text = String::from_utf8_lossy(&line).trim_end().to_string();
            line.clear();
            rewritten = b == b'\r';
            if text.trim().is_empty() {
                continue;
            }
            if bar {
                match download_size(&text) {
                    Some(size) if !downloading => {
                        downloading = true;
                        let _ = tx.send(Event::Stage(format!(
                            "Downloading model weights ({size}), needed once"
                        )));
                    }
                    None if downloading => {
                        downloading = false;
                        if !label.is_empty() {
                            let _ = tx.send(Event::Stage(label.clone()));
                        }
                    }
                    _ => {}
                }
            }
            let event = if bar {
                Event::Progress(job::parse(&text))
            } else {
                Event::Log(text)
            };
            if tx.send(event).is_err() {
                return;
            }
        }
    }
    let rest = String::from_utf8_lossy(&line).trim().to_string();
    if !rest.is_empty() {
        let _ = tx.send(Event::Log(rest));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// tqdm rewrites one line with carriage returns and ends it only when
    /// done; each rewrite must arrive as it happens, and a download must say
    /// it is one.
    #[test]
    fn bars_arrive_without_newlines() {
        let out = "Loading\n  1%|  | 30.8M/5.29G [01:03<4:31:57, 346kB/s]\r  2%|  | 61.0M/5.29G [01:30<4:00:00, 350kB/s]\r";
        let (tx, rx) = std::sync::mpsc::channel::<Event>();
        pump(out.as_bytes(), tx, "Folding".into());
        let events: Vec<String> = rx
            .try_iter()
            .map(|e| match e {
                Event::Stage(s) => format!("stage {s}"),
                Event::Log(l) => format!("log {l}"),
                Event::Progress(p) => format!("progress {:?}", p.fraction),
                _ => String::new(),
            })
            .collect();
        assert_eq!(events[0], "log Loading");
        assert_eq!(events[1], "stage Downloading model weights (5.29GB), needed once");
        assert_eq!(events[2], "progress Some(0.01)");
        assert_eq!(events[3], "progress Some(0.02)");
        assert_eq!(events.len(), 4, "{events:?}");
    }

    #[test]
    fn downloads_are_told_from_work() {
        let esm = "  1%|          | 30.8M/5.29G [01:03<4:31:57, 346kB/s]";
        assert_eq!(download_size(esm).as_deref(), Some("5.29GB"));
        let curl = " 42  400M   42  168M    0     0  12.1M      0  0:00:33  0:00:13  0:00:20 12.5M";
        assert_eq!(download_size(curl).as_deref(), Some("400M"));
        let fold = " 30%|###       | 60/200 [00:07<00:16,  8.43it/s]";
        assert_eq!(download_size(fold), None);
    }

    /// A closed tab drops its end of the channel; the job's subprocess must
    /// be stopped then, not left running.
    #[test]
    fn a_job_nobody_listens_to_is_stopped() {
        let (tx, rx) = std::sync::mpsc::channel::<Event>();
        drop(rx);
        let started = std::time::Instant::now();
        let result = run(Path::new("/bin/sleep"), &["30"], Path::new("/"), &[], &tx);
        assert!(result.is_err_and(|e| e.to_string() == "cancelled"));
        assert!(started.elapsed() < std::time::Duration::from_secs(5));
    }
}
