//! Structures built into the app, so there is something to look at before
//! any data of your own: one of each kind of thing the viewer is for.

pub struct Sample {
    /// The heading it is listed under. Samples of a group sit together.
    pub group: &'static str,
    /// The file name, which is how its format is told.
    pub file: &'static str,
    pub name: &'static str,
    pub note: &'static str,
    pub bytes: &'static [u8],
}

pub const SAMPLES: [Sample; 9] = [
    Sample {
        group: "Molecules",
        file: "caffeine.sdf",
        name: "Caffeine",
        note: "a small molecule, SDF",
        bytes: include_bytes!("../samples/caffeine.sdf"),
    },
    Sample {
        group: "Molecules",
        file: "c60.xyz",
        name: "Buckminsterfullerene",
        note: "C₆₀, XYZ",
        bytes: include_bytes!("../samples/c60.xyz"),
    },
    Sample {
        group: "Proteins",
        file: "1crn.pdb.gz",
        name: "Crambin",
        note: "a protein, PDB",
        bytes: include_bytes!("../samples/1crn.pdb.gz"),
    },
    Sample {
        group: "Proteins",
        file: "1ubq.pdb.gz",
        name: "Ubiquitin",
        note: "a small protein fold, PDB",
        bytes: include_bytes!("../samples/1ubq.pdb.gz"),
    },
    Sample {
        group: "Proteins",
        file: "1ema.pdb.gz",
        name: "Green fluorescent protein",
        note: "a β-barrel, PDB",
        bytes: include_bytes!("../samples/1ema.pdb.gz"),
    },
    Sample {
        group: "Proteins",
        file: "4hhb.pdb.gz",
        name: "Haemoglobin",
        note: "four chains and their haems, PDB",
        bytes: include_bytes!("../samples/4hhb.pdb.gz"),
    },
    Sample {
        group: "Crystals",
        file: "rocksalt.cif",
        name: "Rock salt",
        note: "a crystal from its symmetry, CIF",
        bytes: include_bytes!("../samples/rocksalt.cif"),
    },
    Sample {
        group: "Crystals",
        file: "silicon.vasp",
        name: "Silicon",
        note: "a diamond-cubic supercell, VASP",
        bytes: include_bytes!("../samples/silicon.vasp"),
    },
    Sample {
        group: "Simulations",
        file: "lj-melt.lammpstrj.gz",
        name: "Lennard-Jones melt",
        note: "a crystal melting, LAMMPS trajectory",
        bytes: include_bytes!("../samples/lj-melt.lammpstrj.gz"),
    },
];

/// The samples under their headings, in the order they are listed, each
/// with its index in [`SAMPLES`].
pub fn grouped() -> Vec<(&'static str, Vec<(usize, &'static Sample)>)> {
    let mut groups: Vec<(&'static str, Vec<(usize, &'static Sample)>)> = Vec::new();
    for (i, s) in SAMPLES.iter().enumerate() {
        match groups.last_mut() {
            Some((g, list)) if *g == s.group => list.push((i, s)),
            _ => groups.push((s.group, vec![(i, s)])),
        }
    }
    groups
}

#[cfg(test)]
mod tests {
    #[test]
    fn every_sample_opens() {
        for s in &super::SAMPLES {
            let structure = crate::formats::load_bytes(s.file, s.bytes.to_vec())
                .unwrap_or_else(|e| panic!("{}: {e:#}", s.file));
            assert!(!structure.atoms.is_empty(), "{}", s.file);
        }
    }

    #[test]
    fn groups_are_contiguous() {
        let groups = super::grouped();
        for (i, (g, _)) in groups.iter().enumerate() {
            assert!(!groups[..i].iter().any(|(h, _)| h == g), "{g} is split");
        }
    }
}
