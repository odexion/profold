//! Downloads that are all or nothing.
//!
//! SimpleFold fetches its weights with `curl -o` straight to the final file,
//! and its helpers with `urlretrieve` the same way, so an interrupted fold
//! leaves a truncated checkpoint where a whole one is expected. The next run
//! finds the file, believes it, and fails to load it, or worse. So everything
//! a fold needs is fetched here first, and SimpleFold finds it in place.
//!
//! A file is written beside its destination as `<name>.part`, checked against
//! the size the server declared, flushed to disk, and only then renamed into
//! place. A rename within a directory is atomic, so the destination is either
//! absent or whole. An interrupted download keeps its `.part` and picks up
//! where it stopped, with an HTTP range request, the next time. Several
//! gigabytes are not fetched twice because a tab was closed.

use std::collections::BTreeMap;
use std::fs::{File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc::Sender;
use std::time::{Duration, Instant};

use anyhow::{anyhow, bail, Context, Result};

use crate::job::{self, Event, Progress};

/// A file a fold needs: where it comes from, where it goes, and what to call
/// it while it downloads.
pub struct Wanted {
    pub url: String,
    pub path: PathBuf,
    pub what: String,
}

/// Make sure each file is present and whole, downloading what is missing.
pub fn ensure(files: &[Wanted], ledger: &Path, tx: &Sender<Event>) -> Result<()> {
    for file in files {
        fetch(file, ledger, tx)?;
    }
    Ok(())
}

/// The sizes of the files known to be whole, kept in `ledger`, so that a
/// file checked once is not checked against the server on every fold.
fn read_ledger(ledger: &Path) -> BTreeMap<String, u64> {
    std::fs::read_to_string(ledger)
        .ok()
        .and_then(|t| serde_json::from_str::<serde_json::Value>(&t).ok())
        .and_then(|v| v.as_object().cloned())
        .map(|m| {
            m.into_iter()
                .filter_map(|(k, v)| Some((k, v.as_u64()?)))
                .collect()
        })
        .unwrap_or_default()
}

fn record(ledger: &Path, path: &Path, size: u64) {
    let mut known = read_ledger(ledger);
    known.insert(path.to_string_lossy().into_owned(), size);
    let json = serde_json::to_string_pretty(&known).unwrap_or_default();
    // The ledger is replaced the same way as the files it describes.
    let part = ledger.with_extension("json.part");
    if std::fs::write(&part, json).is_ok() {
        let _ = std::fs::rename(&part, ledger);
    }
}

fn part_of(path: &Path) -> PathBuf {
    let mut name = path.file_name().unwrap_or_default().to_os_string();
    name.push(".part");
    path.with_file_name(name)
}

fn agent() -> ureq::Agent {
    ureq::Agent::config_builder()
        .timeout_connect(Some(Duration::from_secs(30)))
        // A stalled transfer fails rather than hanging the tab for ever; its
        // `.part` is kept and resumed next time.
        .timeout_recv_body(None)
        .timeout_recv_response(Some(Duration::from_secs(60)))
        .build()
        .new_agent()
}

/// The size the server says the file has.
fn remote_size(url: &str) -> Result<u64> {
    let resp = agent()
        .head(url)
        .header("user-agent", "profold/0.1")
        .call()
        .with_context(|| format!("asking {url} for its size"))?;
    header_u64(&resp, "content-length").ok_or_else(|| anyhow!("{url} did not report a file size"))
}

fn header_u64<B>(resp: &ureq::http::Response<B>, name: &str) -> Option<u64> {
    resp.headers().get(name)?.to_str().ok()?.trim().parse().ok()
}

fn fetch(file: &Wanted, ledger: &Path, tx: &Sender<Event>) -> Result<()> {
    let path = &file.path;
    let part = part_of(path);

    if let Ok(meta) = std::fs::metadata(path) {
        let have = meta.len();
        if read_ledger(ledger).get(&*path.to_string_lossy()) == Some(&have) {
            return Ok(());
        }
        // There from before profold checked downloads, perhaps from an
        // interrupted curl. Ask the server how big it should be.
        match remote_size(&file.url) {
            Ok(size) if size == have => {
                record(ledger, path, have);
                return Ok(());
            }
            // A truncated file is a good start: it becomes the part, and the
            // rest is fetched after it.
            Ok(size) if have < size => {
                let _ = tx.send(Event::Log(format!(
                    "{} is incomplete ({} of {}); finishing it.",
                    path.display(),
                    job::size(have),
                    job::size(size)
                )));
                std::fs::rename(path, &part)?;
            }
            Ok(_) => std::fs::remove_file(path)?,
            // Offline: a fold that has everything on disk should still run.
            Err(_) => return Ok(()),
        }
    }

    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let start = std::fs::metadata(&part).map(|m| m.len()).unwrap_or(0);
    let mut request = agent()
        .get(&file.url)
        .header("user-agent", "profold/0.1")
        // Compressed bytes would not line up with the file's own.
        .header("accept-encoding", "identity");
    if start > 0 {
        request = request.header("range", format!("bytes={start}-"));
    }
    let mut resp = request
        .call()
        .with_context(|| format!("downloading {}", file.what))?;

    // 206 continues the part; a server that ignores the range sends it all
    // again with 200, and the part starts over.
    let resumed = resp.status() == 206;
    let offset = if resumed { start } else { 0 };
    let total = if resumed {
        resp.headers()
            .get("content-range")
            .and_then(|v| v.to_str().ok())
            .and_then(|v| v.rsplit('/').next())
            .and_then(|v| v.trim().parse::<u64>().ok())
    } else {
        header_u64(&resp, "content-length")
    };

    let _ = tx.send(Event::Stage(match total {
        Some(t) => format!("Downloading {} ({}), needed once", file.what, job::size(t)),
        None => format!("Downloading {}, needed once", file.what),
    }));
    if resumed {
        let _ = tx.send(Event::Log(format!(
            "Resuming {} from {}.",
            file.what,
            job::size(start)
        )));
    }

    let mut out: File = if resumed {
        OpenOptions::new().append(true).open(&part)?
    } else {
        File::create(&part)?
    };
    let mut reader = resp.body_mut().with_config().limit(u64::MAX).reader();
    let mut buf = vec![0u8; 1 << 20];
    let mut written = offset;
    let mut last = Instant::now();
    loop {
        let n = reader
            .read(&mut buf)
            .with_context(|| format!("downloading {}", file.what))?;
        if n == 0 {
            break;
        }
        out.write_all(&buf[..n])?;
        written += n as u64;
        if last.elapsed() > Duration::from_millis(120) {
            last = Instant::now();
            // Nobody listening means the tab was closed: stop, and leave the
            // part to be resumed.
            if tx
                .send(Event::Progress(Progress::bytes(written, total)))
                .is_err()
            {
                out.flush()?;
                bail!("cancelled");
            }
        }
    }
    out.sync_all()?;
    drop(out);

    if let Some(total) = total {
        if written != total {
            bail!(
                "{} stopped at {} of {}; it will resume next time",
                file.what,
                job::size(written),
                job::size(total)
            );
        }
    }
    std::fs::rename(&part, path)?;
    record(ledger, path, written);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parts_sit_beside_their_file() {
        assert_eq!(
            part_of(Path::new("/x/simplefold_100M.ckpt")),
            PathBuf::from("/x/simplefold_100M.ckpt.part")
        );
    }

    #[test]
    fn the_ledger_round_trips() {
        let dir = std::env::temp_dir().join(format!("profold-ledger-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let ledger = dir.join("downloads.json");
        record(&ledger, Path::new("/a/b.ckpt"), 42);
        record(&ledger, Path::new("/a/c.ckpt"), 7);
        let known = read_ledger(&ledger);
        assert_eq!(known.get("/a/b.ckpt"), Some(&42));
        assert_eq!(known.get("/a/c.ckpt"), Some(&7));
        let _ = std::fs::remove_dir_all(dir);
    }

    /// Live: fetch a small file in two halves, as an interrupted download
    /// would, and check it lands whole. `cargo test resumes_live -- --ignored`.
    #[test]
    #[ignore]
    fn resumes_live() {
        let url = "https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t36_3B_UR50D-contact-regression.pt";
        let dir = std::env::temp_dir().join(format!("profold-dl-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("regression.pt");
        let ledger = dir.join("downloads.json");
        // The first half, as a part left by an interrupted run.
        let mut whole = Vec::new();
        agent()
            .get(url)
            .call()
            .unwrap()
            .body_mut()
            .as_reader()
            .read_to_end(&mut whole)
            .unwrap();
        std::fs::write(part_of(&path), &whole[..whole.len() / 2]).unwrap();

        let (tx, rx) = std::sync::mpsc::channel();
        let wanted = Wanted {
            url: url.into(),
            path: path.clone(),
            what: "test".into(),
        };
        ensure(&[wanted], &ledger, &tx).unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), whole);
        assert!(!part_of(&path).exists());
        assert!(rx
            .try_iter()
            .any(|e| matches!(e, Event::Log(l) if l.starts_with("Resuming"))));
        let _ = std::fs::remove_dir_all(dir);
    }
}
