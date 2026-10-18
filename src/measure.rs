//! Distances, angles and dihedrals between atoms the user picks.

use crate::protein::Structure;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    Distance,
    Angle,
    Dihedral,
}

impl Kind {
    pub const ALL: [Kind; 3] = [Kind::Distance, Kind::Angle, Kind::Dihedral];

    pub fn label(self) -> &'static str {
        match self {
            Kind::Distance => "Distance",
            Kind::Angle => "Angle",
            Kind::Dihedral => "Dihedral",
        }
    }

    /// How many atoms it takes.
    pub fn atoms(self) -> usize {
        match self {
            Kind::Distance => 2,
            Kind::Angle => 3,
            Kind::Dihedral => 4,
        }
    }

    /// The figure for these points, one per atom: ångström for a distance,
    /// degrees for the others.
    pub fn value(self, p: &[[f32; 3]]) -> f32 {
        match self {
            Kind::Distance => len(sub(p[1], p[0])),
            Kind::Angle => {
                let (a, b) = (sub(p[0], p[1]), sub(p[2], p[1]));
                let c = dot(a, b) / (len(a) * len(b)).max(1e-12);
                c.clamp(-1.0, 1.0).acos().to_degrees()
            }
            // IUPAC's sign: looking down the middle bond, a clockwise turn
            // from the first atom to the last is positive.
            Kind::Dihedral => {
                let (b1, b2, b3) = (sub(p[1], p[0]), sub(p[2], p[1]), sub(p[3], p[2]));
                let y = len(b2) * dot(b1, cross(b2, b3));
                let x = dot(cross(b1, b2), cross(b2, b3));
                y.atan2(x).to_degrees()
            }
        }
    }

    pub fn format(self, v: f32) -> String {
        match self {
            Kind::Distance => format!("{v:.2} Å"),
            _ => format!("{v:.1}°"),
        }
    }
}

/// One measurement, between atoms of one structure.
#[derive(Clone, PartialEq, Debug)]
pub struct Measurement {
    pub kind: Kind,
    /// Which structure: 0 is the one loaded, 1 the prediction laid over it.
    pub which: usize,
    pub atoms: Vec<usize>,
}

impl Measurement {
    /// Where its atoms are now, which a trajectory moves frame by frame.
    pub fn points(&self, s: &Structure) -> Option<Vec<[f32; 3]>> {
        self.atoms.iter().map(|&i| s.atoms.get(i).map(|a| a.pos)).collect()
    }

    pub fn value(&self, s: &Structure) -> Option<f32> {
        Some(self.kind.value(&self.points(s)?))
    }

    /// The segments drawn for it: atom to atom, in order.
    pub fn segments(&self, s: &Structure) -> Vec<([f32; 3], [f32; 3])> {
        let p = self.points(s).unwrap_or_default();
        p.windows(2).map(|w| (w[0], w[1])).collect()
    }

    /// Where its figure is written: halfway along a distance, at an angle's
    /// vertex, and halfway along a dihedral's middle bond.
    pub fn anchor(&self, s: &Structure) -> Option<[f32; 3]> {
        let p = self.points(s)?;
        Some(match self.kind {
            Kind::Distance => mid(p[0], p[1]),
            Kind::Angle => p[1],
            Kind::Dihedral => mid(p[1], p[2]),
        })
    }
}

fn sub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn len(a: [f32; 3]) -> f32 {
    dot(a, a).sqrt()
}

fn mid(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5, (a[2] + b[2]) * 0.5]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f32, b: f32) -> bool {
        (a - b).abs() < 1e-3
    }

    #[test]
    fn distances_and_angles() {
        assert!(close(Kind::Distance.value(&[[0.0; 3], [3.0, 4.0, 0.0]]), 5.0));
        let right = [[1.0, 0.0, 0.0], [0.0; 3], [0.0, 2.0, 0.0]];
        assert!(close(Kind::Angle.value(&right), 90.0));
        let straight = [[1.0, 0.0, 0.0], [0.0; 3], [-3.0, 0.0, 0.0]];
        assert!(close(Kind::Angle.value(&straight), 180.0));
    }

    #[test]
    fn dihedrals_turn_clockwise_positive() {
        let (b, c) = ([0.0; 3], [0.0, 0.0, 1.0]);
        let a = [1.0, 0.0, 0.0];
        // Looking down +z, from x round to y is clockwise.
        assert!(close(Kind::Dihedral.value(&[a, b, c, [0.0, 1.0, 1.0]]), 90.0));
        assert!(close(Kind::Dihedral.value(&[a, b, c, [0.0, -1.0, 1.0]]), -90.0));
        assert!(close(Kind::Dihedral.value(&[a, b, c, [1.0, 0.0, 1.0]]), 0.0));
        assert!(close(Kind::Dihedral.value(&[a, b, c, [-1.0, 0.0, 1.0]]).abs(), 180.0));
    }

    #[test]
    fn a_helix_has_helical_backbone_angles() {
        // Crambin's first helix runs from residue 7 to 19: its phi near -60°.
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("samples/1crn.pdb.gz");
        let s = crate::formats::load_bytes("1crn.pdb.gz", std::fs::read(path).unwrap()).unwrap();
        let r = |i: usize| &s.residues[i];
        let i = 10;
        let phi = Measurement {
            kind: Kind::Dihedral,
            which: 0,
            atoms: vec![
                r(i - 1).c.unwrap() as usize,
                r(i).n.unwrap() as usize,
                r(i).ca.unwrap() as usize,
                r(i).c.unwrap() as usize,
            ],
        };
        let v = phi.value(&s).unwrap();
        assert!((-80.0..-45.0).contains(&v), "phi {v}");
    }
}
