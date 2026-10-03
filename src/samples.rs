//! Structures built into the app, so there is something to look at before
//! any data of your own: one of each kind of thing the viewer is for.

pub struct Sample {
    /// The file name, which is how its format is told.
    pub file: &'static str,
    pub name: &'static str,
    pub note: &'static str,
    pub bytes: &'static [u8],
}

pub const SAMPLES: [Sample; 6] = [
    Sample {
        file: "caffeine.sdf",
        name: "Caffeine",
        note: "a small molecule, SDF",
        bytes: include_bytes!("../samples/caffeine.sdf"),
    },
    Sample {
        file: "c60.xyz",
        name: "Buckminsterfullerene",
        note: "C₆₀, XYZ",
        bytes: include_bytes!("../samples/c60.xyz"),
    },
    Sample {
        file: "1crn.pdb.gz",
        name: "Crambin",
        note: "a protein, PDB",
        bytes: include_bytes!("../samples/1crn.pdb.gz"),
    },
    Sample {
        file: "rocksalt.cif",
        name: "Rock salt",
        note: "a crystal from its symmetry, CIF",
        bytes: include_bytes!("../samples/rocksalt.cif"),
    },
    Sample {
        file: "silicon.vasp",
        name: "Silicon",
        note: "a diamond-cubic supercell, VASP",
        bytes: include_bytes!("../samples/silicon.vasp"),
    },
    Sample {
        file: "lj-melt.lammpstrj.gz",
        name: "Lennard-Jones melt",
        note: "a crystal melting, LAMMPS trajectory",
        bytes: include_bytes!("../samples/lj-melt.lammpstrj.gz"),
    },
];

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
}
