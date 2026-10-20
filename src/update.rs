//! Updates from the project's GitHub releases.
//!
//! A release carries one archive per platform, `profold-<target>.tar.gz`, and
//! its SHA-256 beside it. The updater asks GitHub for the latest release,
//! downloads the archive for the target this build was made for, checks it,
//! and puts the new app where the running one is: the whole `profold.app` on
//! macOS, the binary on Linux and Windows. The running process keeps its old
//! files open, so it carries on until it is restarted into the new ones.
//!
//! Only release builds replace themselves. CI sets `PROFOLD_TARGET` when it
//! builds one; a build from source has no target, and is told about a new
//! release without being overwritten by it.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::{Duration, SystemTime, UNIX_EPOCH};

use anyhow::{anyhow, bail, Context, Result};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// The Rust target a release build was made for, set by CI.
pub const TARGET: Option<&str> = option_env!("PROFOLD_TARGET");

pub const REPO: &str = "odexion/profold";

/// Automatic checks are at most this often.
const EVERY: Duration = Duration::from_secs(24 * 60 * 60);

/// Archives run to a few tens of megabytes; this is a ceiling, not a guess.
const MAX_ARCHIVE: u64 = 512 * 1024 * 1024;

/// Where the updater has got to.
#[derive(Clone, Debug, PartialEq)]
pub enum Status {
    Idle,
    Checking,
    UpToDate,
    /// A newer release this build cannot install itself: a build from
    /// source, or one installed somewhere it cannot write.
    Available {
        version: String,
        page: String,
    },
    Downloading {
        version: String,
        fraction: Option<f32>,
    },
    Installing {
        version: String,
    },
    /// Installed; it takes over at the next start.
    Ready {
        version: String,
    },
    Failed(String),
}

impl Status {
    pub fn busy(&self) -> bool {
        matches!(
            self,
            Status::Checking | Status::Downloading { .. } | Status::Installing { .. }
        )
    }
}

pub type Shared = Arc<Mutex<Status>>;

fn set(status: &Shared, s: Status) {
    *status.lock().unwrap() = s;
}

// ---------------------------------------------------------------------------
// Preferences
// ---------------------------------------------------------------------------

/// Kept apart from `settings.json`, which each tab writes in full.
fn prefs_path() -> Option<PathBuf> {
    crate::settings::file("update.json")
}

fn prefs() -> Value {
    prefs_path()
        .and_then(|p| std::fs::read_to_string(p).ok())
        .and_then(|t| serde_json::from_str(&t).ok())
        .unwrap_or_else(|| json!({}))
}

fn save_prefs(v: &Value) {
    let Some(p) = prefs_path() else { return };
    if let Some(dir) = p.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let _ = std::fs::write(p, serde_json::to_string_pretty(v).unwrap_or_default());
}

/// Whether to look for a new release on launch. On unless turned off.
pub fn automatic() -> bool {
    prefs()
        .get("automatic")
        .and_then(Value::as_bool)
        .unwrap_or(true)
}

pub fn set_automatic(on: bool) {
    let mut v = prefs();
    v["automatic"] = json!(on);
    save_prefs(&v);
}

fn now() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

/// Is an automatic check due?
pub fn due() -> bool {
    if cfg!(test) || std::env::var_os("PROFOLD_NO_UPDATE").is_some() {
        return false;
    }
    let last = prefs().get("checked").and_then(Value::as_u64).unwrap_or(0);
    automatic() && now().saturating_sub(last) >= EVERY.as_secs()
}

fn mark_checked() {
    let mut v = prefs();
    v["checked"] = json!(now());
    save_prefs(&v);
}

// ---------------------------------------------------------------------------
// Versions and releases
// ---------------------------------------------------------------------------

/// `v1.2.3` or `1.2.3` as numbers. A pre-release suffix is ignored, since
/// GitHub's latest release is never a pre-release.
pub fn parse_version(s: &str) -> Option<(u64, u64, u64)> {
    let s = s.trim().trim_start_matches('v');
    let core = s.split(['-', '+']).next()?;
    let mut parts = core.split('.').map(|p| p.parse::<u64>().ok());
    let v = (
        parts.next()??,
        parts.next().unwrap_or(Some(0))?,
        parts.next().unwrap_or(Some(0))?,
    );
    Some(v)
}

pub fn is_newer(candidate: &str, current: &str) -> bool {
    match (parse_version(candidate), parse_version(current)) {
        (Some(a), Some(b)) => a > b,
        _ => false,
    }
}

pub fn archive_name(target: &str) -> String {
    format!("profold-{target}.tar.gz")
}

#[derive(Debug)]
pub struct Release {
    pub version: String,
    pub page: String,
    /// The archive for this build's target, and its checksum file.
    pub archive: Option<(String, String)>,
}

fn agent_get(url: &str) -> Result<ureq::http::Response<ureq::Body>> {
    ureq::get(url)
        .header("user-agent", &format!("profold/{VERSION}"))
        .header("accept", "application/vnd.github+json")
        .call()
        .with_context(|| format!("requesting {url}"))
}

pub fn latest() -> Result<Release> {
    // Overridable, to try an update against a release served locally.
    let url = std::env::var("PROFOLD_RELEASES_URL")
        .unwrap_or_else(|_| format!("https://api.github.com/repos/{REPO}/releases/latest"));
    let text = agent_get(&url)?
        .body_mut()
        .read_to_string()
        .context("reading the release")?;
    parse_release(&text, TARGET)
}

fn parse_release(text: &str, target: Option<&str>) -> Result<Release> {
    let v: Value = serde_json::from_str(text).context("the release is not JSON")?;
    let tag = v
        .get("tag_name")
        .and_then(Value::as_str)
        .ok_or_else(|| anyhow!("the release has no tag"))?;
    let page = v
        .get("html_url")
        .and_then(Value::as_str)
        .unwrap_or("https://github.com/odexion/profold/releases/latest")
        .to_string();
    let assets: Vec<(&str, &str)> = v
        .get("assets")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(|a| {
            Some((
                a.get("name")?.as_str()?,
                a.get("browser_download_url")?.as_str()?,
            ))
        })
        .collect();
    let find = |name: &str| {
        assets
            .iter()
            .find(|(n, _)| *n == name)
            .map(|(_, u)| u.to_string())
    };
    let archive = target.and_then(|t| {
        let name = archive_name(t);
        Some((find(&name)?, find(&format!("{name}.sha256"))?))
    });
    Ok(Release {
        version: tag.trim_start_matches('v').to_string(),
        page,
        archive,
    })
}

fn download(url: &str, mut progress: impl FnMut(Option<f32>)) -> Result<Vec<u8>> {
    let mut resp = agent_get(url)?;
    let total: Option<u64> = resp
        .headers()
        .get("content-length")
        .and_then(|v| v.to_str().ok())
        .and_then(|v| v.parse().ok());
    let mut reader = resp.body_mut().with_config().limit(MAX_ARCHIVE).reader();
    let mut out = Vec::with_capacity(total.unwrap_or(1 << 20) as usize);
    let mut buf = vec![0u8; 256 * 1024];
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        out.extend_from_slice(&buf[..n]);
        progress(
            total
                .filter(|t| *t > 0)
                .map(|t| out.len() as f32 / t as f32),
        );
    }
    Ok(out)
}

/// The first word of a `sha256sum` line is the digest.
fn verify(bytes: &[u8], sums: &str) -> Result<()> {
    let want = sums
        .split_whitespace()
        .next()
        .ok_or_else(|| anyhow!("the checksum file is empty"))?
        .to_ascii_lowercase();
    let got: String = Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();
    if got != want {
        bail!("the download is damaged (checksum mismatch)");
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Installing
// ---------------------------------------------------------------------------

/// What the updater replaces.
#[derive(Clone, Debug, PartialEq)]
pub enum Layout {
    /// A macOS app bundle, replaced whole.
    Bundle(PathBuf),
    /// A bare executable.
    Binary(PathBuf),
}

impl Layout {
    pub fn path(&self) -> &Path {
        match self {
            Layout::Bundle(p) | Layout::Binary(p) => p,
        }
    }

    /// The name of what replaces it inside a release archive.
    fn entry(&self) -> &'static str {
        match self {
            Layout::Bundle(_) => "profold.app",
            Layout::Binary(_) if cfg!(windows) => "profold.exe",
            Layout::Binary(_) => "profold",
        }
    }
}

/// Worked out once, at start: on Linux the running executable's path reads
/// `… (deleted)` once it has been replaced.
pub fn layout() -> Option<&'static Layout> {
    static LAYOUT: OnceLock<Option<Layout>> = OnceLock::new();
    LAYOUT
        .get_or_init(|| {
            let exe = std::env::current_exe().ok()?;
            // Through any symlink to the real file; not on Windows, where
            // that gives a `\\?\` path other programs choke on.
            let exe = if cfg!(windows) {
                exe
            } else {
                exe.canonicalize().ok()?
            };
            let bundle = exe
                .parent()
                .and_then(Path::parent)
                .and_then(Path::parent)
                .filter(|p| p.extension().is_some_and(|e| e == "app"));
            Some(match bundle {
                Some(app) if cfg!(target_os = "macos") => Layout::Bundle(app.to_path_buf()),
                _ => Layout::Binary(exe),
            })
        })
        .as_ref()
}

fn staging(layout: &Layout) -> PathBuf {
    // Beside what it replaces, so the swap is a rename on one file system.
    layout.path().with_file_name(".profold-update")
}

fn aside(layout: &Layout) -> PathBuf {
    let path = layout.path();
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    path.with_file_name(format!(".{name}.old"))
}

fn remove(path: &Path) -> std::io::Result<()> {
    match std::fs::symlink_metadata(path) {
        Ok(m) if m.is_dir() => std::fs::remove_dir_all(path),
        Ok(_) => std::fs::remove_file(path),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(e),
    }
}

/// Put the app in `archive` in place of the one at `layout`.
pub fn install(layout: &Layout, archive: &[u8]) -> Result<()> {
    let stage = staging(layout);
    remove(&stage).ok();
    std::fs::create_dir_all(&stage).with_context(|| {
        format!(
            "cannot write to {}; reinstall profold with the installer",
            stage.parent().unwrap_or(&stage).display()
        )
    })?;
    let result = (|| {
        tar::Archive::new(flate2::read::GzDecoder::new(archive))
            .unpack(&stage)
            .context("unpacking the update")?;
        let new = stage.join(layout.entry());
        if !new.exists() {
            bail!("the update has no {}", layout.entry());
        }
        swap(layout.path(), &new, &aside(layout))
    })();
    remove(&stage).ok();
    result
}

/// Move `current` aside to `old`, and `new` into its place, putting
/// `current` back if the second step fails.
fn swap(current: &Path, new: &Path, old: &Path) -> Result<()> {
    remove(old).ok();
    std::fs::rename(current, old).with_context(|| format!("moving {} aside", current.display()))?;
    if let Err(e) = std::fs::rename(new, current) {
        let _ = std::fs::rename(old, current);
        return Err(e).with_context(|| format!("putting the update at {}", current.display()));
    }
    // Windows will not delete an executable that is running; it goes at the
    // next start instead.
    remove(old).ok();
    Ok(())
}

/// Clear what an earlier update left behind.
pub fn tidy() {
    if let Some(layout) = layout() {
        remove(&staging(layout)).ok();
        remove(&aside(layout)).ok();
    }
}

/// Check for a release and, where this build can, install it. Progress goes
/// into `status`, which the window reads.
pub fn run(status: Shared) {
    set(&status, Status::Checking);
    let result = (|| -> Result<Status> {
        let release = latest()?;
        mark_checked();
        if !is_newer(&release.version, VERSION) {
            return Ok(Status::UpToDate);
        }
        let version = release.version.clone();
        let (Some((url, sums_url)), Some(layout)) = (release.archive.as_ref(), layout()) else {
            return Ok(Status::Available {
                version,
                page: release.page,
            });
        };
        set(
            &status,
            Status::Downloading {
                version: version.clone(),
                fraction: None,
            },
        );
        let sums = String::from_utf8(download(sums_url, |_| {})?)?;
        let mut last = std::time::Instant::now();
        let bytes = download(url, |fraction| {
            if last.elapsed() > Duration::from_millis(100) {
                last = std::time::Instant::now();
                set(
                    &status,
                    Status::Downloading {
                        version: version.clone(),
                        fraction,
                    },
                );
            }
        })?;
        verify(&bytes, &sums)?;
        set(
            &status,
            Status::Installing {
                version: version.clone(),
            },
        );
        install(layout, &bytes)?;
        Ok(Status::Ready { version })
    })();
    set(
        &status,
        result.unwrap_or_else(|e| Status::Failed(format!("{e:#}"))),
    );
}

/// Start the new version once this one has gone. GPUI's own restart does
/// this on macOS and Windows; on Linux it does not quote the path.
pub fn restart(cx: &mut gpui::App) {
    #[cfg(target_os = "linux")]
    if let Some(layout) = layout() {
        use std::os::unix::process::CommandExt;
        let spawned = std::process::Command::new("/bin/sh")
            .arg("-c")
            .arg(r#"while kill -0 "$0" 2>/dev/null; do sleep 0.1; done; exec "$1""#)
            .arg(std::process::id().to_string())
            .arg(layout.path())
            .process_group(0)
            .spawn();
        if spawned.is_ok() {
            cx.quit();
            return;
        }
    }
    if let Some(layout) = layout() {
        cx.set_restart_path(layout.path().to_path_buf());
    }
    cx.restart();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn versions_compare_as_numbers() {
        assert_eq!(parse_version("v0.10.2"), Some((0, 10, 2)));
        assert_eq!(parse_version("1.2"), Some((1, 2, 0)));
        assert_eq!(parse_version("1.2.3-rc.1"), Some((1, 2, 3)));
        assert_eq!(parse_version("nightly"), None);
        assert!(is_newer("v0.10.0", "0.9.9"));
        assert!(!is_newer("0.1.0", "0.1.0"));
        assert!(!is_newer("garbage", "0.1.0"));
    }

    #[test]
    fn a_release_names_the_archive_for_its_target() {
        let text = r#"{
            "tag_name": "v0.2.0",
            "html_url": "https://github.com/odexion/profold/releases/tag/v0.2.0",
            "assets": [
                {"name": "profold-x86_64-unknown-linux-gnu.tar.gz", "browser_download_url": "https://x/a"},
                {"name": "profold-x86_64-unknown-linux-gnu.tar.gz.sha256", "browser_download_url": "https://x/a.sha256"},
                {"name": "profold-aarch64-apple-darwin.tar.gz", "browser_download_url": "https://x/b"}
            ]
        }"#;
        let r = parse_release(text, Some("x86_64-unknown-linux-gnu")).unwrap();
        assert_eq!(r.version, "0.2.0");
        assert_eq!(
            r.archive,
            Some(("https://x/a".into(), "https://x/a.sha256".into()))
        );
        // No checksum, no archive: an unverifiable download is not installed.
        let r = parse_release(text, Some("aarch64-apple-darwin")).unwrap();
        assert_eq!(r.archive, None);
        assert_eq!(parse_release(text, None).unwrap().archive, None);
    }

    #[test]
    fn a_damaged_download_is_refused() {
        let sum: String = Sha256::digest(b"hello")
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        assert!(verify(b"hello", &format!("{sum}  profold.tar.gz\n")).is_ok());
        assert!(verify(b"hellp", &sum).is_err());
    }

    fn archive_of(name: &str, body: &[u8]) -> Vec<u8> {
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        {
            let mut tar = tar::Builder::new(&mut gz);
            let mut header = tar::Header::new_gnu();
            header.set_size(body.len() as u64);
            header.set_mode(0o755);
            header.set_cksum();
            tar.append_data(&mut header, name, body).unwrap();
            tar.finish().unwrap();
        }
        gz.finish().unwrap()
    }

    #[test]
    fn an_update_replaces_the_binary_and_leaves_nothing_behind() {
        let dir = std::env::temp_dir().join(format!("profold-update-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let exe = dir.join(if cfg!(windows) {
            "profold.exe"
        } else {
            "profold"
        });
        std::fs::write(&exe, b"old").unwrap();
        let layout = Layout::Binary(exe.clone());

        install(&layout, &archive_of(layout.entry(), b"new")).unwrap();
        assert_eq!(std::fs::read(&exe).unwrap(), b"new");
        let left: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .map(|e| e.unwrap().file_name())
            .collect();
        assert_eq!(left.len(), 1, "left behind: {left:?}");

        // An archive without the app in it changes nothing.
        assert!(install(&layout, &archive_of("other", b"x")).is_err());
        assert_eq!(std::fs::read(&exe).unwrap(), b"new");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn an_update_replaces_a_whole_bundle() {
        let dir = std::env::temp_dir().join(format!("profold-bundle-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let app = dir.join("profold.app");
        std::fs::create_dir_all(app.join("Contents/MacOS")).unwrap();
        std::fs::write(app.join("Contents/MacOS/profold"), b"old").unwrap();
        std::fs::write(app.join("Contents/stale"), b"gone after").unwrap();

        let archive = archive_of("profold.app/Contents/MacOS/profold", b"new");
        install(&Layout::Bundle(app.clone()), &archive).unwrap();
        assert_eq!(
            std::fs::read(app.join("Contents/MacOS/profold")).unwrap(),
            b"new"
        );
        assert!(!app.join("Contents/stale").exists());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
