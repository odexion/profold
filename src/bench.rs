//! How long the operations on a structure take, at the sizes they meet:
//! `cargo test --release ops_bench -- --ignored --nocapture`. For A/B tests
//! of the periodic tools, the analysis, properties, selections and the
//! isosurface; the renderer has its own in `render.rs`.

use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::Arc;
use std::time::{Duration, Instant};

use crate::analysis::{measure, Measure};
use crate::protein::{Property, Snapshot, Structure};

/// The best of a few runs, which is what is not noise.
fn time(label: &str, runs: usize, mut f: impl FnMut()) -> Duration {
    let mut best = Duration::MAX;
    for _ in 0..runs {
        let t = Instant::now();
        f();
        best = best.min(t.elapsed());
    }
    eprintln!("{label:<44} {best:>10.2?}");
    best
}

/// `s` with `n` frames, each atom jiggled about where it was by a fixed
/// pseudo-random walk, so runs are alike.
fn jiggled(s: &Structure, n: usize, step: f32) -> Structure {
    let mut seed = 0x2545_f491_4f6c_dd1du64;
    let mut next = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 40) as f32 / (1u64 << 24) as f32 - 0.5
    };
    let mut pos: Vec<[f32; 3]> = s.atoms.iter().map(|a| a.pos).collect();
    let frames: Vec<Snapshot> = (0..n)
        .map(|_| {
            for p in pos.iter_mut() {
                for x in p.iter_mut() {
                    *x += next() * step;
                }
            }
            Snapshot {
                step: None,
                pos: pos.clone(),
                cell: s.cell,
            }
        })
        .collect();
    let mut out = s.with_frames(frames).expect("as many atoms");
    out.cell = s.cell;
    out
}

#[test]
#[ignore]
fn ops_bench() {
    let polymer = jiggled(
        &crate::formats::load_bytes(
            "4hhb.pdb.gz",
            include_bytes!("../samples/4hhb.pdb.gz").to_vec(),
        )
        .unwrap(),
        200,
        0.15,
    );
    let denge =
        std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/data.denge")).expect("data.denge");
    let mut sim = jiggled(
        &crate::formats::load_bytes("data.denge", denge).unwrap(),
        50,
        0.3,
    );
    let n = sim.atoms.len();
    let per_frame: Vec<Vec<f32>> = sim
        .frames
        .iter()
        .map(|f| f.pos.iter().map(|p| p[0] - p[2]).collect())
        .collect();
    sim.props = vec![Property {
        name: "c_pe".into(),
        values: per_frame[0].clone(),
        frames: Arc::new(per_frame),
    }];
    eprintln!(
        "\nhaemoglobin: {} atoms × {} frames; data.denge: {} atoms × {} frames\n",
        polymer.atoms.len(),
        polymer.frame_count(),
        n,
        sim.frame_count()
    );

    let go = |s: &Structure, m: Measure, atoms: &[usize]| {
        measure(s, m, atoms, &AtomicBool::new(false), &AtomicUsize::new(0));
    };
    let all = |s: &Structure| (0..s.atoms.len()).collect::<Vec<_>>();
    let ca: Vec<usize> = polymer
        .residues
        .iter()
        .filter_map(|r| r.ca.map(|a| a as usize))
        .collect();
    for (m, runs) in [
        (Measure::Rmsd, 3),
        (Measure::Rg, 3),
        (Measure::Rmsf, 3),
        (Measure::Sasa, 1),
        (Measure::Contacts, 3),
        (Measure::Secondary, 1),
    ] {
        time(&format!("haemoglobin {} of Cα", m.label()), runs, || {
            go(&polymer, m, &ca)
        });
        time(
            &format!("haemoglobin {} of all atoms", m.label()),
            runs,
            || go(&polymer, m, &all(&polymer)),
        );
    }
    for (m, runs) in [
        (Measure::Rmsd, 3),
        (Measure::Rg, 3),
        (Measure::Rmsf, 3),
        (Measure::Sasa, 1),
    ] {
        time(
            &format!("data.denge {} of all atoms", m.label()),
            runs,
            || go(&sim, m, &all(&sim)),
        );
    }
    time("haemoglobin fit frames", 3, || {
        let mut s = polymer.clone();
        crate::analysis::fit_frames(&mut s).unwrap();
    });
    time("data.denge fit frames", 3, || {
        let mut s = sim.clone();
        crate::analysis::fit_frames(&mut s).unwrap();
    });

    time("data.denge wrap", 3, || {
        crate::pbc::wrap(&mut sim.clone()).unwrap()
    });
    time("data.denge make whole", 3, || {
        crate::pbc::make_whole(&mut sim.clone()).unwrap()
    });
    time("data.denge unwrap in time", 3, || {
        crate::pbc::unwrap_in_time(&mut sim.clone()).unwrap()
    });
    time("data.denge replicate 2×2×2", 3, || {
        crate::pbc::replicate(&sim, [2, 2, 2]).unwrap();
    });
    time("data.denge clone (the floor for the above)", 3, || {
        let _ = sim.clone();
    });

    time("data.denge property span, every frame", 3, || {
        let mut s = sim.clone();
        s.forget_derived();
        s.prop_span(0);
    });
    time("data.denge select c_pe < 0 and z > 0", 3, || {
        crate::select::Selection::parse("c_pe < 0 and z > 0")
            .unwrap()
            .atoms(&sim);
    });
    time("data.denge colour by property (build scene)", 3, || {
        let mut s = sim.clone();
        s.prop_color = Some(crate::protein::PropertyColor {
            index: 0,
            range: None,
        });
        crate::render::build_scene(
            &s,
            crate::render::Repr::Spacefill,
            crate::render::Coloring::Property,
            1.0,
            crate::render::Style::Soft,
            crate::render::Palette::Viridis,
            true,
            false,
            &crate::render::Appearance::default(),
        );
    });

    let k = 200;
    let volume = crate::volume::Volume {
        origin: [0.0; 3],
        axes: [[0.25, 0.0, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 0.25]],
        dims: [k; 3],
        data: (0..k * k * k)
            .map(|i| {
                let (x, y, z) = ((i % k) as f32, ((i / k) % k) as f32, (i / k / k) as f32);
                ((x * 0.11).sin() + (y * 0.07).cos() + (z * 0.05).sin()) * (z - 100.0) / 100.0
            })
            .collect(),
    };
    time("200³ volume default level", 3, || {
        volume.default_level();
    });
    time("200³ volume signed?", 3, || {
        volume.signed();
    });
    let level = volume.default_level();
    let mut facets = 0;
    time("200³ volume isosurface, one lobe", 3, || {
        facets = volume.isosurface(level, false).len()
    });
    eprintln!("  ({facets} triangles)");
}

/// One operation on the large dump in the repository's root, for its time
/// and, under `/usr/bin/time -l`, its peak memory:
/// `PROFOLD_OP=rmsf cargo test --release dump_bench -- --ignored --nocapture`.
/// `PROFOLD_OP` is load, same, rmsd, rmsf, rg, wrap, whole, unwrap, fit or span.
#[test]
#[ignore]
fn dump_bench() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/dump.lammpstrj");
    let op = std::env::var("PROFOLD_OP").unwrap_or_else(|_| "load".into());
    let t = Instant::now();
    let mut s =
        crate::formats::load_path(std::path::Path::new(path)).expect("dump.lammpstrj in the root");
    eprintln!(
        "load: {} atoms × {} frames in {:.2?}",
        s.atoms.len(),
        s.frame_count(),
        t.elapsed()
    );
    let all: Vec<usize> = (0..s.atoms.len()).collect();
    let go = |s: &Structure, m: Measure| {
        measure(s, m, &all, &AtomicBool::new(false), &AtomicUsize::new(0)).is_some()
    };
    let t = Instant::now();
    match op.as_str() {
        "load" => {}
        // Streamed and read whole, the same to the last bit.
        "same" => {
            let whole =
                crate::formats::load_bytes("dump.lammpstrj", std::fs::read(path).unwrap()).unwrap();
            assert_eq!(whole.frame_count(), s.frame_count());
            assert!(whole
                .frames
                .iter()
                .zip(s.frames.iter())
                .all(|(a, b)| a.pos == b.pos && a.step == b.step));
            eprintln!("streamed and whole agree on {} frames", s.frame_count());
        }
        "rmsd" => assert!(go(&s, Measure::Rmsd)),
        "rmsf" => assert!(go(&s, Measure::Rmsf)),
        "rg" => assert!(go(&s, Measure::Rg)),
        "wrap" => crate::pbc::wrap(&mut s).unwrap(),
        "whole" => crate::pbc::make_whole(&mut s).unwrap(),
        "unwrap" => crate::pbc::unwrap_in_time(&mut s).unwrap(),
        "fit" => crate::analysis::fit_frames(&mut s).unwrap(),
        "span" => {
            s.props = vec![Property {
                name: "z".into(),
                values: s.atoms.iter().map(|a| a.pos[2]).collect(),
                frames: Arc::new(
                    s.frames
                        .iter()
                        .map(|f| f.pos.iter().map(|p| p[2]).collect())
                        .collect(),
                ),
            }];
            s.prop_span(0);
        }
        other => panic!("no operation {other}"),
    }
    eprintln!("{op}: {:.2?}", t.elapsed());
}

/// The features that draw and read, at the sizes they meet:
/// `cargo test --release features_bench -- --ignored --nocapture`.
#[test]
#[ignore]
fn features_bench() {
    use crate::render::{build_scene, Appearance, Camera, Coloring, Palette, Repr, Style};
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/dump.lammpstrj");
    let dump =
        crate::formats::load_path(std::path::Path::new(path)).expect("dump.lammpstrj in the root");
    let hb = crate::formats::load_bytes(
        "4hhb.pdb.gz",
        include_bytes!("../samples/4hhb.pdb.gz").to_vec(),
    )
    .unwrap();
    let page = [0.97, 0.96, 0.94];
    let control = crate::trace::Control::default();

    // Ray tracing: the still view of a retina window, and saved images.
    for (name, s, repr) in [
        ("haemoglobin cartoon", &hb, Repr::Cartoon),
        ("dump spacefill", &dump, Repr::Spacefill),
    ] {
        let scene = build_scene(
            s,
            repr,
            Coloring::Chain,
            1.0,
            Style::Studio,
            Palette::Spectrum,
            true,
            true,
            &Appearance::default(),
        );
        let cam = Camera::frame(s);
        time(&format!("raster {name}, 1616×1436 (CPU)"), 3, || {
            crate::render::render_cpu(&scene, &cam, 1616, 1436, page, Style::Studio);
        });
        for (label, w, h) in [
            ("view 1616×1436", 1616, 1436),
            ("image 2× 3232×2872", 3232, 2872),
            ("image 3× 4848×4308", 4848, 4308),
            ("film 1080p", 1920, 1080),
        ] {
            time(&format!("trace {name}, {label}"), 1, || {
                crate::trace::render(&scene, &cam, w, h, page, Style::Studio, &control).unwrap();
            });
        }
    }

    // Labels drawn into a 4K frame.
    let font = crate::movie::font().unwrap();
    let cam = Camera::frame(&hb);
    let labels: Vec<([f32; 3], String)> = hb
        .atoms
        .iter()
        .step_by(hb.atoms.len() / 200)
        .map(|a| (a.pos, "HIS 87:A".to_string()))
        .collect();
    let labels: Vec<([f32; 3], &str)> = labels.iter().map(|(p, t)| (*p, t.as_str())).collect();
    let mut pic = crate::movie::Picture {
        width: 3840,
        height: 2160,
        rgba: vec![255; 3840 * 2160 * 4],
    };
    time(
        &format!("{} labels into a 4K frame", labels.len()),
        3,
        || {
            crate::movie::draw_labels(&mut pic, &font, &cam, &labels, 26.0, [0.0; 3], page);
        },
    );

    // Volumes: reading a 200³ grid as cube text, as an MRC map, and as OpenDX.
    let k = 200usize;
    let value = |i: usize, j: usize, l: usize| {
        ((i as f32 * 0.11).sin() + (j as f32 * 0.07).cos()) * (l as f32 - 100.0) / 100.0
    };
    let mut cube = format!(
        "big\ntest\n    1 0 0 0\n   {k} 0.2 0 0\n   {k} 0 0.2 0\n   {k} 0 0 0.2\n    6 6.0 0 0 0\n"
    );
    let mut dx = format!("object 1 class gridpositions counts {k} {k} {k}\norigin 0 0 0\ndelta 0.1 0 0\ndelta 0 0.1 0\ndelta 0 0 0.1\nobject 3 class array type double rank 0 items {} data follows\n", k * k * k);
    let mut mrc = vec![0u8; 1024];
    for (o, v) in [
        (0, k as i32),
        (4, k as i32),
        (8, k as i32),
        (12, 2),
        (28, k as i32),
        (32, k as i32),
        (36, k as i32),
        (64, 1),
        (68, 2),
        (72, 3),
    ] {
        mrc[o..o + 4].copy_from_slice(&v.to_le_bytes());
    }
    for (o, v) in [(40, 20.0f32), (44, 20.0), (48, 20.0)] {
        mrc[o..o + 4].copy_from_slice(&v.to_le_bytes());
    }
    for i in 0..k {
        for j in 0..k {
            for l in 0..k {
                let v = value(i, j, l);
                cube.push_str(&format!("{v:.5e} "));
                dx.push_str(&format!("{v:.5e} "));
            }
            cube.push('\n');
            dx.push('\n');
        }
    }
    for l in 0..k {
        for j in 0..k {
            for i in 0..k {
                mrc.extend(value(i, j, l).to_le_bytes());
            }
        }
    }
    eprintln!(
        "  (cube text {} MB, MRC {} MB)",
        cube.len() >> 20,
        mrc.len() >> 20
    );
    time("read 200³ cube file", 1, || {
        crate::formats::load_bytes("big.cube", cube.clone().into_bytes()).unwrap();
    });
    time("read 200³ OpenDX map", 1, || {
        crate::volume::read_map("big.dx", dx.clone().into_bytes()).unwrap();
    });
    time("read 200³ MRC map", 3, || {
        crate::volume::read_map("big.mrc", mrc.clone()).unwrap();
    });

    // Replicating the dump, every frame.
    time("replicate dump 2×2×2 (801 frames)", 1, || {
        crate::pbc::replicate(&dump, [2, 2, 2]).unwrap();
    });
}

/// Ray tracing alone, for A/B tests of the tracer:
/// `cargo test --release trace_bench -- --ignored --nocapture`. Writes the
/// pictures to `PROFOLD_TRACE_OUT`, when set, to compare before and after.
#[test]
#[ignore]
fn trace_bench() {
    use crate::render::{build_scene, Appearance, Camera, Coloring, Palette, Repr, Style};
    let hb = crate::formats::load_bytes(
        "4hhb.pdb.gz",
        include_bytes!("../samples/4hhb.pdb.gz").to_vec(),
    )
    .unwrap();
    let denge =
        std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/data.denge")).expect("data.denge");
    let sim = crate::formats::load_bytes("data.denge", denge).unwrap();
    let out = std::env::var("PROFOLD_TRACE_OUT").ok();
    for (name, s, repr) in [
        ("haemoglobin-cartoon", &hb, Repr::Cartoon),
        ("haemoglobin-spacefill", &hb, Repr::Spacefill),
        ("haemoglobin-surface", &hb, Repr::Surface),
        ("denge-spacefill", &sim, Repr::Spacefill),
    ] {
        let scene = build_scene(
            s,
            repr,
            Coloring::Chain,
            1.0,
            Style::Studio,
            Palette::Spectrum,
            true,
            true,
            &Appearance::default(),
        );
        let mut cam = Camera::frame(s);
        cam.rotate(0.3, 0.2);
        let mut frame = None;
        time(
            &format!("trace {name}, 1616×1436 ({} prims)", scene.len()),
            2,
            || {
                frame = crate::trace::render(
                    &scene,
                    &cam,
                    1616,
                    1436,
                    [0.97, 0.96, 0.94],
                    Style::Studio,
                    &Default::default(),
                );
            },
        );
        if let (Some(dir), Some(f)) = (&out, frame) {
            let mut rgba = f.bgra;
            for p in rgba.chunks_exact_mut(4) {
                p.swap(0, 2);
            }
            image::save_buffer(
                format!("{dir}/{name}.png"),
                &rgba,
                1616,
                1436,
                image::ColorType::Rgba8,
            )
            .unwrap();
        }
    }
}
