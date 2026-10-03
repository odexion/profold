//! The periodic table, as far as drawing needs it.
//!
//! An element is its atomic number, with zero standing for "unknown": a
//! LAMMPS type with no mass to go on, a coarse-grained bead, a label that
//! names nothing. Unknown atoms are still drawn, coloured by their type.

#[derive(Clone, Copy, PartialEq, Eq, Debug, Default, Hash)]
pub struct Element(pub u8);

/// Symbol, standard atomic weight, covalent radius (Cordero 2008), van der
/// Waals radius (Bondi where he gives one, Alvarez otherwise), Jmol colour.
const TABLE: [(&str, f32, f32, f32, u32); 104] = [
    ("X", 0.0, 0.0, 1.60, 0x9a8fb3),
    ("H", 1.008, 0.31, 1.10, 0xffffff),
    ("He", 4.0026, 0.28, 1.40, 0xd9ffff),
    ("Li", 6.94, 1.28, 1.82, 0xcc80ff),
    ("Be", 9.0122, 0.96, 1.53, 0xc2ff00),
    ("B", 10.81, 0.84, 1.92, 0xffb5b5),
    ("C", 12.011, 0.76, 1.70, 0x909090),
    ("N", 14.007, 0.71, 1.55, 0x3050f8),
    ("O", 15.999, 0.66, 1.52, 0xff0d0d),
    ("F", 18.998, 0.57, 1.47, 0x90e050),
    ("Ne", 20.180, 0.58, 1.54, 0xb3e3f5),
    ("Na", 22.990, 1.66, 2.27, 0xab5cf2),
    ("Mg", 24.305, 1.41, 1.73, 0x8aff00),
    ("Al", 26.982, 1.21, 1.84, 0xbfa6a6),
    ("Si", 28.085, 1.11, 2.10, 0xf0c8a0),
    ("P", 30.974, 1.07, 1.80, 0xff8000),
    ("S", 32.06, 1.05, 1.80, 0xffff30),
    ("Cl", 35.45, 1.02, 1.75, 0x1ff01f),
    ("Ar", 39.948, 1.06, 1.88, 0x80d1e3),
    ("K", 39.098, 2.03, 2.75, 0x8f40d4),
    ("Ca", 40.078, 1.76, 2.31, 0x3dff00),
    ("Sc", 44.956, 1.70, 2.15, 0xe6e6e6),
    ("Ti", 47.867, 1.60, 2.11, 0xbfc2c7),
    ("V", 50.942, 1.53, 2.07, 0xa6a6ab),
    ("Cr", 51.996, 1.39, 2.06, 0x8a99c7),
    ("Mn", 54.938, 1.39, 2.05, 0x9c7ac7),
    ("Fe", 55.845, 1.32, 2.04, 0xe06633),
    ("Co", 58.933, 1.26, 2.00, 0xf090a0),
    ("Ni", 58.693, 1.24, 1.63, 0x50d050),
    ("Cu", 63.546, 1.32, 1.40, 0xc88033),
    ("Zn", 65.38, 1.22, 1.39, 0x7d80b0),
    ("Ga", 69.723, 1.22, 1.87, 0xc28f8f),
    ("Ge", 72.630, 1.20, 2.11, 0x668f8f),
    ("As", 74.922, 1.19, 1.85, 0xbd80e3),
    ("Se", 78.971, 1.20, 1.90, 0xffa100),
    ("Br", 79.904, 1.20, 1.85, 0xa62929),
    ("Kr", 83.798, 1.16, 2.02, 0x5cb8d1),
    ("Rb", 85.468, 2.20, 3.03, 0x702eb0),
    ("Sr", 87.62, 1.95, 2.49, 0x00ff00),
    ("Y", 88.906, 1.90, 2.32, 0x94ffff),
    ("Zr", 91.224, 1.75, 2.23, 0x94e0e0),
    ("Nb", 92.906, 1.64, 2.18, 0x73c2c9),
    ("Mo", 95.95, 1.54, 2.17, 0x54b5b5),
    ("Tc", 98.0, 1.47, 2.16, 0x3b9e9e),
    ("Ru", 101.07, 1.46, 2.13, 0x248f8f),
    ("Rh", 102.91, 1.42, 2.10, 0x0a7d8c),
    ("Pd", 106.42, 1.39, 1.63, 0x006985),
    ("Ag", 107.87, 1.45, 1.72, 0xc0c0c0),
    ("Cd", 112.41, 1.44, 1.58, 0xffd98f),
    ("In", 114.82, 1.42, 1.93, 0xa67573),
    ("Sn", 118.71, 1.39, 2.17, 0x668080),
    ("Sb", 121.76, 1.39, 2.06, 0x9e63b5),
    ("Te", 127.60, 1.38, 2.06, 0xd47a00),
    ("I", 126.90, 1.39, 1.98, 0x940094),
    ("Xe", 131.29, 1.40, 2.16, 0x429eb0),
    ("Cs", 132.91, 2.44, 3.43, 0x57178f),
    ("Ba", 137.33, 2.15, 2.68, 0x00c900),
    ("La", 138.91, 2.07, 2.43, 0x70d4ff),
    ("Ce", 140.12, 2.04, 2.42, 0xffffc7),
    ("Pr", 140.91, 2.03, 2.40, 0xd9ffc7),
    ("Nd", 144.24, 2.01, 2.39, 0xc7ffc7),
    ("Pm", 145.0, 1.99, 2.38, 0xa3ffc7),
    ("Sm", 150.36, 1.98, 2.36, 0x8fffc7),
    ("Eu", 151.96, 1.98, 2.35, 0x61ffc7),
    ("Gd", 157.25, 1.96, 2.34, 0x45ffc7),
    ("Tb", 158.93, 1.94, 2.33, 0x30ffc7),
    ("Dy", 162.50, 1.92, 2.31, 0x1fffc7),
    ("Ho", 164.93, 1.92, 2.30, 0x00ff9c),
    ("Er", 167.26, 1.89, 2.29, 0x00e675),
    ("Tm", 168.93, 1.90, 2.27, 0x00d452),
    ("Yb", 173.05, 1.87, 2.26, 0x00bf38),
    ("Lu", 174.97, 1.87, 2.24, 0x00ab24),
    ("Hf", 178.49, 1.75, 2.23, 0x4dc2ff),
    ("Ta", 180.95, 1.70, 2.22, 0x4da6ff),
    ("W", 183.84, 1.62, 2.18, 0x2194d6),
    ("Re", 186.21, 1.51, 2.16, 0x267dab),
    ("Os", 190.23, 1.44, 2.16, 0x266696),
    ("Ir", 192.22, 1.41, 2.13, 0x175487),
    ("Pt", 195.08, 1.36, 1.75, 0xd0d0e0),
    ("Au", 196.97, 1.36, 1.66, 0xffd123),
    ("Hg", 200.59, 1.32, 1.55, 0xb8b8d0),
    ("Tl", 204.38, 1.45, 1.96, 0xa6544d),
    ("Pb", 207.2, 1.46, 2.02, 0x575961),
    ("Bi", 208.98, 1.48, 2.07, 0x9e4fb5),
    ("Po", 209.0, 1.40, 1.97, 0xab5c00),
    ("At", 210.0, 1.50, 2.02, 0x754f45),
    ("Rn", 222.0, 1.50, 2.20, 0x428296),
    ("Fr", 223.0, 2.60, 3.48, 0x420066),
    ("Ra", 226.0, 2.21, 2.83, 0x007d00),
    ("Ac", 227.0, 2.15, 2.47, 0x70abfa),
    ("Th", 232.04, 2.06, 2.45, 0x00baff),
    ("Pa", 231.04, 2.00, 2.43, 0x00a1ff),
    ("U", 238.03, 1.96, 2.41, 0x008fff),
    ("Np", 237.0, 1.90, 2.39, 0x0080ff),
    ("Pu", 244.0, 1.87, 2.43, 0x006bff),
    ("Am", 243.0, 1.80, 2.44, 0x545cf2),
    ("Cm", 247.0, 1.69, 2.45, 0x785ce3),
    ("Bk", 247.0, 1.68, 2.44, 0x8a4fe3),
    ("Cf", 251.0, 1.68, 2.45, 0xa136d4),
    ("Es", 252.0, 1.65, 2.45, 0xb31fd4),
    ("Fm", 257.0, 1.67, 2.45, 0xb31fba),
    ("Md", 258.0, 1.73, 2.46, 0xb30da6),
    ("No", 259.0, 1.76, 2.46, 0xbd0d87),
    ("Lr", 266.0, 1.61, 2.46, 0xc70066),
];

/// Elements that bond covalently. Two metals next to each other are a
/// lattice, not a bond, and drawing every contact in a crystal of iron turns
/// it into a scaffold.
const NONMETALS: [u8; 23] = [1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 51, 52, 53, 54, 85];

impl Element {
    pub const UNKNOWN: Element = Element(0);
    pub const H: Element = Element(1);
    pub const C: Element = Element(6);
    pub const N: Element = Element(7);
    pub const O: Element = Element(8);

    fn row(self) -> &'static (&'static str, f32, f32, f32, u32) {
        TABLE.get(self.0 as usize).unwrap_or(&TABLE[0])
    }

    pub fn is_known(self) -> bool {
        self.0 != 0
    }

    pub fn symbol(self) -> &'static str {
        self.row().0
    }

    pub fn from_number(z: u32) -> Self {
        if (z as usize) < TABLE.len() {
            Element(z as u8)
        } else {
            Element::UNKNOWN
        }
    }

    /// An element symbol in any case: "Fe", "FE", "fe". Deuterium and
    /// tritium are hydrogen.
    pub fn from_symbol(sym: &str) -> Self {
        let sym = sym.trim();
        if sym.eq_ignore_ascii_case("D") || sym.eq_ignore_ascii_case("T") {
            return Element::H;
        }
        TABLE
            .iter()
            .skip(1)
            .position(|row| row.0.eq_ignore_ascii_case(sym))
            .map(|i| Element(i as u8 + 1))
            .unwrap_or(Element::UNKNOWN)
    }

    /// A label from a file that has no element column: "C12", "Fe2", "OW",
    /// "HW1", "1HB". Two-letter symbols win when written in their own case,
    /// so "Ca" is calcium and "CA" an alpha carbon.
    pub fn from_label(label: &str) -> Self {
        let letters: String = label
            .trim()
            .trim_start_matches(|c: char| c.is_ascii_digit())
            .chars()
            .take_while(|c| c.is_ascii_alphabetic())
            .collect();
        let mut chars = letters.chars();
        let Some(first) = chars.next() else {
            return Element::UNKNOWN;
        };
        if let Some(second) = chars.next() {
            if first.is_ascii_uppercase() && second.is_ascii_lowercase() {
                let two = Element::from_symbol(&letters[..2]);
                if two.is_known() {
                    return two;
                }
            }
        }
        Element::from_symbol(&first.to_string())
    }

    /// The element whose standard weight is this mass, when one is close
    /// enough to be sure. A mass of exactly one is the reduced unit of a
    /// Lennard-Jones simulation far more often than it is hydrogen.
    pub fn from_mass(mass: f32) -> Self {
        if !(mass > 0.0) || mass == 1.0 {
            return Element::UNKNOWN;
        }
        TABLE
            .iter()
            .enumerate()
            .skip(1)
            .map(|(z, row)| (z, (row.1 - mass).abs()))
            .filter(|(_, d)| *d < 0.12)
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(z, _)| Element(z as u8))
            .unwrap_or(Element::UNKNOWN)
    }

    pub fn mass(self) -> f32 {
        self.row().1
    }

    /// Covalent radius in ångström; zero when unknown.
    pub fn covalent(self) -> f32 {
        self.row().2
    }

    /// Van der Waals radius in ångström.
    pub fn vdw(self) -> f32 {
        self.row().3
    }

    pub fn is_metal(self) -> bool {
        self.is_known() && !NONMETALS.contains(&self.0)
    }

    /// CPK-ish colours. The organic six keep their own muted shades so a
    /// protein sits well on paper white; everything else is Jmol's colour,
    /// taken a step toward grey for the same reason.
    pub fn color(self) -> [f32; 3] {
        match self.0 {
            6 => [0.42, 0.44, 0.47],
            7 => [0.24, 0.42, 0.78],
            8 => [0.83, 0.28, 0.26],
            16 => [0.85, 0.68, 0.20],
            15 => [0.90, 0.51, 0.20],
            1 => [0.85, 0.85, 0.87],
            _ => {
                let hex = self.row().4;
                let c = [
                    ((hex >> 16) & 0xff) as f32 / 255.0,
                    ((hex >> 8) & 0xff) as f32 / 255.0,
                    (hex & 0xff) as f32 / 255.0,
                ];
                let grey = 0.55;
                [
                    c[0] * 0.82 + grey * 0.18,
                    c[1] * 0.82 + grey * 0.18,
                    c[2] * 0.82 + grey * 0.18,
                ]
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symbols_and_labels() {
        assert_eq!(Element::from_symbol("FE"), Element(26));
        assert_eq!(Element::from_symbol("d"), Element::H);
        assert_eq!(Element::from_label("Ca"), Element(20));
        assert_eq!(Element::from_label("CA"), Element::C);
        assert_eq!(Element::from_label("OW"), Element::O);
        assert_eq!(Element::from_label("1HB"), Element::H);
        assert_eq!(Element::from_label("Fe2"), Element(26));
        assert_eq!(Element(26).symbol(), "Fe");
        assert_eq!(Element(103).symbol(), "Lr");
    }

    #[test]
    fn masses() {
        assert_eq!(Element::from_mass(12.011), Element::C);
        assert_eq!(Element::from_mass(63.55), Element(29));
        assert_eq!(Element::from_mass(1.0), Element::UNKNOWN);
        assert_eq!(Element::from_mass(3.0), Element::UNKNOWN);
    }
}
