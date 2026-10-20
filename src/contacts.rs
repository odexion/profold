//! Contacts between a ligand and the residues around it: hydrogen bonds,
//! salt bridges, metal coordination, and clashes.
//!
//! Most files carry no hydrogens, so a hydrogen bond is read from its heavy
//! atoms alone, as PyMOL's polar contacts are: a nitrogen or oxygen within
//! 3.5 Å of another.

use crate::protein::{Element, Structure};

#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub enum Kind {
    /// A metal ion held by a nitrogen, oxygen or sulphur.
    Metal,
    /// Opposite charges within reach of each other.
    SaltBridge,
    /// A nitrogen or oxygen within hydrogen-bonding distance of another.
    HBond,
    /// Two atoms closer than their radii allow.
    Clash,
}

impl Kind {
    pub const ALL: [Kind; 4] = [Kind::HBond, Kind::SaltBridge, Kind::Metal, Kind::Clash];

    pub fn label(self, n: usize) -> &'static str {
        match (self, n == 1) {
            (Kind::HBond, true) => "hydrogen bond",
            (Kind::HBond, false) => "hydrogen bonds",
            (Kind::SaltBridge, true) => "salt bridge",
            (Kind::SaltBridge, false) => "salt bridges",
            (Kind::Metal, true) => "metal contact",
            (Kind::Metal, false) => "metal contacts",
            (Kind::Clash, true) => "clash",
            (Kind::Clash, false) => "clashes",
        }
    }

    /// The colour its line is drawn in, as 0xRRGGBB.
    pub fn color(self) -> u32 {
        match self {
            Kind::HBond => 0x4fc3f7,
            Kind::SaltBridge => 0xe056c1,
            Kind::Metal => 0xb388ff,
            Kind::Clash => 0xff4d4d,
        }
    }
}

/// One contact, between atom `a` of the ligand and atom `b` beside it.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Contact {
    pub kind: Kind,
    pub a: usize,
    pub b: usize,
    pub distance: f32,
}

const HBOND: f32 = 3.5;
const SALT: f32 = 4.0;
const METAL: f32 = 2.8;
/// Closer than the two radii less this counts as a clash. Heavy atoms only,
/// so it is generous: a clash with hydrogens would be called sooner.
const CLASH_OVERLAP: f32 = 0.6;

/// Whether an element is a metal that binds as an ion.
fn metal(e: Element) -> bool {
    matches!(
        e.0,
        3 | 11 | 12 | 19 | 20 | 24..=30 | 37 | 38 | 42 | 47 | 48 | 55 | 56 | 74 | 78 | 79 | 80
    )
}

fn polar(e: Element) -> bool {
    matches!(e.0, 7 | 8)
}

/// The charge an atom carries at neutral pH, read from what it is: the side
/// chains of Asp, Glu, Lys and Arg, and in a ligand the oxygens of
/// phosphates, sulphates and carboxylates.
fn charge(s: &Structure, i: usize, bonded: &[Vec<usize>]) -> i8 {
    let a = &s.atoms[i];
    let r = &s.residues[a.residue as usize];
    if !r.outside_polymer() {
        return match (r.name.as_str(), a.name()) {
            ("ASP", "OD1" | "OD2") | ("GLU", "OE1" | "OE2") => -1,
            ("LYS", "NZ") | ("ARG", "NH1" | "NH2" | "NE") => 1,
            _ => 0,
        };
    }
    if a.element.0 != 8 {
        return 0;
    }
    // An oxygen hung on a phosphorus or sulphur, or one of a carboxylate's
    // two, with no other neighbour. A bond from oxygen to oxygen is a file's
    // slip, which old entries' CONECT records make, and does not count.
    let heavy: Vec<usize> = bonded[i]
        .iter()
        .copied()
        .filter(|&j| s.atoms[j].element.0 != 8)
        .collect();
    let [centre] = heavy[..] else { return 0 };
    let oxygens = bonded[centre]
        .iter()
        .filter(|&&j| s.atoms[j].element.0 == 8)
        .count();
    match s.atoms[centre].element.0 {
        15 | 16 => -1,
        6 if oxygens == 2 => -1,
        _ => 0,
    }
}

/// Every contact between residue `ligand` and residues `around`.
pub fn between(s: &Structure, ligand: usize, around: &[usize]) -> Vec<Contact> {
    let mut bonded = vec![Vec::new(); s.atoms.len()];
    for &[a, b] in &s.bonds {
        bonded[a as usize].push(b as usize);
        bonded[b as usize].push(a as usize);
    }
    let atoms_of = s.atoms_by_residue();
    let own: Vec<usize> = atoms_of.of(ligand).iter().map(|&i| i as usize).collect();
    let others: Vec<usize> = around
        .iter()
        .flat_map(|&r| atoms_of.of(r))
        .map(|&i| i as usize)
        .collect();
    let mut out = Vec::new();
    for &a in &own {
        let (pa, ea) = (s.atoms[a].pos, s.atoms[a].element);
        if ea == Element::H {
            continue;
        }
        for &b in &others {
            let (pb, eb) = (s.atoms[b].pos, s.atoms[b].element);
            // A metal's coordination is often taken for a bond by distance;
            // it is still the contact to show.
            if eb == Element::H || (bonded[a].contains(&b) && !metal(ea) && !metal(eb)) {
                continue;
            }
            let d = ((pa[0] - pb[0]).powi(2) + (pa[1] - pb[1]).powi(2) + (pa[2] - pb[2]).powi(2))
                .sqrt();
            if d > SALT {
                continue;
            }
            let kind = if (metal(ea) && (polar(eb) || eb.0 == 16))
                || (metal(eb) && (polar(ea) || ea.0 == 16))
            {
                (d <= METAL).then_some(Kind::Metal)
            } else if charge(s, a, &bonded) * charge(s, b, &bonded) < 0 {
                Some(Kind::SaltBridge)
            } else if polar(ea) && polar(eb) && d <= HBOND {
                Some(Kind::HBond)
            } else if d < ea.vdw() + eb.vdw() - CLASH_OVERLAP {
                Some(Kind::Clash)
            } else {
                None
            };
            if let Some(kind) = kind {
                out.push(Contact {
                    kind,
                    a,
                    b,
                    distance: d,
                });
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(name: &str) -> Structure {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("samples")
            .join(name);
        crate::formats::load_bytes(name, std::fs::read(path).unwrap()).unwrap()
    }

    #[test]
    fn a_haem_iron_sits_on_its_histidine() {
        let s = sample("4hhb.pdb.gz");
        let haem = s.copies("HEM")[0];
        let site = s.site(haem, 5.0);
        let contacts = between(&s, haem, &site);
        // The proximal histidine's NE2 holds the iron.
        let iron = contacts
            .iter()
            .find(|c| c.kind == Kind::Metal)
            .expect("the iron is held by nothing");
        let b = &s.atoms[iron.b];
        assert_eq!(s.atoms[iron.a].element.0, 26);
        assert_eq!(
            (s.residues[b.residue as usize].name.as_str(), b.name()),
            ("HIS", "NE2")
        );
        assert!(iron.distance < 2.5);
        // The propionates reach out to the solvent and to basic side chains.
        assert!(contacts
            .iter()
            .any(|c| c.kind == Kind::HBond || c.kind == Kind::SaltBridge));
        assert!(contacts.iter().all(|c| c.distance <= SALT));
    }

    #[test]
    fn a_carboxylate_oxygen_is_charged() {
        let s = sample("4hhb.pdb.gz");
        let mut bonded = vec![Vec::new(); s.atoms.len()];
        for &[a, b] in &s.bonds {
            bonded[a as usize].push(b as usize);
            bonded[b as usize].push(a as usize);
        }
        // A haem's propionates end in carboxylates; its ring has no oxygen.
        let haem = s.copies("HEM")[0];
        let oxygens: Vec<usize> = s
            .residue_atoms(haem)
            .into_iter()
            .filter(|&i| s.atoms[i].element.0 == 8)
            .collect();
        assert_eq!(oxygens.len(), 4);
        assert!(oxygens.iter().all(|&i| charge(&s, i, &bonded) == -1));
        // And aspartate's side chain, not its backbone.
        let asp = (0..s.residues.len())
            .find(|&r| s.residues[r].name == "ASP")
            .unwrap();
        for i in s.residue_atoms(asp) {
            let expected = if matches!(s.atoms[i].name(), "OD1" | "OD2") {
                -1
            } else {
                0
            };
            assert_eq!(charge(&s, i, &bonded), expected, "{}", s.atoms[i].name());
        }
    }
}
