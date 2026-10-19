//! A small software rasterizer for molecules.
//!
//! GPUI has no 3D pipeline, so we build the frame ourselves, on the graphics
//! card where there is one (see `gpu.rs`) and otherwise here: orthographic
//! camera, sphere impostors with a z-buffer, then two screen-space passes —
//! depth cueing and a depth-discontinuity outline — that give the result the
//! flat, printed look of a molecular illustration.

use std::collections::HashMap;

use crate::protein::{is_solvent, LigandPaint, Ss, Structure};

/// What an atom is, as far as its drawing settings go: its element, or for
/// an atom of no known element, the name of its type.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum Species {
    Element(u8),
    Kind(String),
}

/// How one species is drawn: a multiple of its usual size, how opaque, and
/// a colour of its own in place of the element's when colouring by element,
/// and on the ligands beside a polymer whatever the colouring.
/// Size or opacity at zero hides it, and its bonds with it.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Look {
    pub size: f32,
    pub alpha: f32,
    pub color: Option<[f32; 3]>,
}

impl Default for Look {
    fn default() -> Self {
        Look {
            size: 1.0,
            alpha: 1.0,
            color: None,
        }
    }
}

impl Look {
    pub fn hidden(self) -> bool {
        self.size <= 0.0 || self.alpha <= 0.0
    }
}

/// The size and opacity set for each species. Kept by element and by type
/// name rather than per structure, so hiding hydrogens stays hidden from one
/// file to the next.
#[derive(Clone, Default, PartialEq)]
pub struct Appearance {
    looks: HashMap<Species, Look>,
    /// Water and other solvent is hidden, as it is by default beside a
    /// protein, where a crystal's waters would only be clutter.
    pub hide_solvent: bool,
    /// Atoms hidden from every style, the backbone's included, by index.
    pub hidden: Option<std::sync::Arc<Vec<bool>>>,
}

impl Appearance {
    pub fn get(&self, species: &Species) -> Look {
        self.looks.get(species).copied().unwrap_or_default()
    }

    pub fn set(&mut self, species: Species, look: Look) {
        if look == Look::default() {
            self.looks.remove(&species);
        } else {
            self.looks.insert(species, look);
        }
    }

    pub fn is_default(&self) -> bool {
        self.looks.is_empty()
    }

    pub fn reset(&mut self) {
        self.looks.clear();
    }

    /// Give every species back its element's colour, keeping sizes and
    /// opacities.
    pub fn clear_colors(&mut self) {
        for look in self.looks.values_mut() {
            look.color = None;
        }
        self.looks.retain(|_, look| *look != Look::default());
    }

    /// Every species that differs from the default, for saving.
    pub fn entries(&self) -> impl Iterator<Item = (&Species, &Look)> {
        self.looks.iter()
    }

    /// One look per atom of `s`.
    fn resolve(&self, s: &Structure) -> Vec<Look> {
        let mut looks = self.resolve_species(s);
        if let Some(mask) = self.mask(s) {
            for (look, &h) in looks.iter_mut().zip(mask) {
                if h {
                    look.size = 0.0;
                }
            }
        }
        if self.hide_solvent {
            for (look, a) in looks.iter_mut().zip(&s.atoms) {
                if is_solvent(&s.residues[a.residue as usize].name) {
                    look.size = 0.0;
                }
            }
        }
        looks
    }

    /// The hidden atoms, when they are this structure's.
    fn mask(&self, s: &Structure) -> Option<&[bool]> {
        self.hidden.as_deref().map(Vec::as_slice).filter(|m| m.len() == s.atoms.len())
    }

    fn resolve_species(&self, s: &Structure) -> Vec<Look> {
        if self.looks.is_empty() {
            return vec![Look::default(); s.atoms.len()];
        }
        let elements: Vec<Look> = (0..=255u8)
            .map(|z| self.get(&Species::Element(z)))
            .collect();
        let kinds: Vec<Look> = s
            .kinds
            .iter()
            .map(|k| self.get(&Species::Kind(k.clone())))
            .collect();
        s.atoms
            .iter()
            .map(|a| {
                if a.element.is_known() {
                    elements[a.element.0 as usize]
                } else {
                    kinds.get(a.kind as usize).copied().unwrap_or_default()
                }
            })
            .collect()
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Repr {
    /// Ribbons for helices, arrows for strands, cord for the rest.
    Cartoon,
    /// Smooth tube through the alpha carbons.
    Tube,
    /// Small atoms joined by thin bonds.
    BallAndStick,
    /// Bonds as rods of one thickness, atoms only as their rounded ends.
    Licorice,
    /// Every atom at its van der Waals radius.
    Spacefill,
    /// The molecule's skin, smooth over its atoms; ligands stand inside it as
    /// atoms.
    Surface,
}

impl Repr {
    pub const ALL: [Repr; 6] = [
        Repr::Cartoon,
        Repr::Tube,
        Repr::BallAndStick,
        Repr::Licorice,
        Repr::Spacefill,
        Repr::Surface,
    ];
    pub fn label(self) -> &'static str {
        match self {
            Repr::Cartoon => "Cartoon",
            Repr::Tube => "Tube",
            Repr::BallAndStick => "Ball & stick",
            Repr::Licorice => "Licorice",
            Repr::Spacefill => "Spacefill",
            Repr::Surface => "Surface",
        }
    }

    /// Drawn from the protein backbone, so meaningless without one.
    pub fn needs_polymer(self) -> bool {
        matches!(self, Repr::Cartoon | Repr::Tube)
    }

    /// What a structure opens as: the cartoon for a protein, ball and stick
    /// for a molecule, and spheres for a simulation whose atoms are anonymous
    /// beads or too many to make out bonds between.
    pub fn default_for(s: &Structure) -> Repr {
        if s.has_polymer() {
            Repr::Cartoon
        } else if s.atoms.iter().all(|a| a.element.is_known()) && s.atoms.len() <= 50_000 {
            Repr::BallAndStick
        } else {
            Repr::Spacefill
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Coloring {
    Chain,
    /// Blue at the N terminus through to red at the C terminus.
    Rainbow,
    /// pLDDT, for predicted models.
    Confidence,
    /// Low to high B-factor along the palette's ramp: how still each part
    /// of a crystal structure sat.
    BFactor,
    /// Helix, strand and coil, each its own colour.
    Structure,
    /// Each residue by how much it shuns water, on Kyte and Doolittle's scale.
    Hydrophobic,
    /// Acidic residues red and basic blue; over a surface, the electrostatic
    /// potential there.
    Charge,
    Element,
    /// One colour per structure, when two are superimposed.
    Pair,
    /// Distance between the superimposed counterparts, per residue.
    Deviation,
}

impl Coloring {
    pub const ALL: [Coloring; 8] = [
        Coloring::Chain,
        Coloring::Rainbow,
        Coloring::Confidence,
        Coloring::BFactor,
        Coloring::Structure,
        Coloring::Hydrophobic,
        Coloring::Charge,
        Coloring::Element,
    ];
    /// The options that mean something while two structures are overlaid.
    pub const COMPARE: [Coloring; 4] = [
        Coloring::Pair,
        Coloring::Deviation,
        Coloring::Rainbow,
        Coloring::Element,
    ];
    pub fn label(self) -> &'static str {
        match self {
            Coloring::Chain => "Chain",
            Coloring::Rainbow => "Rainbow",
            Coloring::Confidence => "pLDDT",
            Coloring::BFactor => "B-factor",
            Coloring::Structure => "Secondary",
            Coloring::Hydrophobic => "Hydrophobic",
            Coloring::Charge => "Charge",
            Coloring::Element => "Element",
            Coloring::Pair => "Pair",
            Coloring::Deviation => "Deviation",
        }
    }
}

/// The colour ramps a structure can be painted with.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Palette {
    /// The classic molecular rainbow.
    Spectrum,
    /// A rainbow without the muddy band in the middle.
    Turbo,
    /// Perceptually even and safe for colour vision deficiency.
    Viridis,
    /// Even and readable to every kind of colour vision.
    Cividis,
    Magma,
    Inferno,
    Plasma,
    Mako,
    Rocket,
    Ember,
    Ice,
    Forest,
    Dusk,
    /// Soft, chalky colour, the way a molecule is painted in an illustration.
    Pastel,
    /// For print, and for a figure that has to survive a photocopier.
    Greyscale,
}

/// The families the palettes are shown in.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PaletteGroup {
    /// Rainbows, which tell neighbours apart best.
    Classic,
    /// Even in lightness, so they read the same to every eye and in grey.
    Perceptual,
    /// A single mood each.
    Thematic,
    /// Quiet colour for illustration and print.
    Muted,
}

impl PaletteGroup {
    pub fn label(self) -> &'static str {
        match self {
            PaletteGroup::Classic => "Classic",
            PaletteGroup::Perceptual => "Perceptual",
            PaletteGroup::Thematic => "Thematic",
            PaletteGroup::Muted => "Muted",
        }
    }
}

impl Palette {
    /// In menu order, each group together.
    pub const ALL: [Palette; 15] = [
        Palette::Spectrum,
        Palette::Turbo,
        Palette::Viridis,
        Palette::Cividis,
        Palette::Magma,
        Palette::Inferno,
        Palette::Plasma,
        Palette::Mako,
        Palette::Rocket,
        Palette::Ember,
        Palette::Ice,
        Palette::Forest,
        Palette::Dusk,
        Palette::Pastel,
        Palette::Greyscale,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Palette::Spectrum => "Spectrum",
            Palette::Turbo => "Turbo",
            Palette::Viridis => "Viridis",
            Palette::Cividis => "Cividis",
            Palette::Magma => "Magma",
            Palette::Inferno => "Inferno",
            Palette::Plasma => "Plasma",
            Palette::Mako => "Mako",
            Palette::Rocket => "Rocket",
            Palette::Ember => "Ember",
            Palette::Ice => "Ice",
            Palette::Forest => "Forest",
            Palette::Dusk => "Dusk",
            Palette::Pastel => "Pastel",
            Palette::Greyscale => "Greyscale",
        }
    }

    pub fn group(self) -> PaletteGroup {
        match self {
            Palette::Spectrum | Palette::Turbo => PaletteGroup::Classic,
            Palette::Viridis
            | Palette::Cividis
            | Palette::Magma
            | Palette::Inferno
            | Palette::Plasma
            | Palette::Mako
            | Palette::Rocket => PaletteGroup::Perceptual,
            Palette::Ember | Palette::Ice | Palette::Forest | Palette::Dusk => {
                PaletteGroup::Thematic
            }
            Palette::Pastel | Palette::Greyscale => PaletteGroup::Muted,
        }
    }

    /// Position `t` in 0..=1 along the ramp.
    pub fn ramp(self, t: f32) -> [f32; 3] {
        let t = t.clamp(0.0, 1.0);
        match self {
            // Hue sweep rather than a table, so it stays smooth end to end.
            Palette::Spectrum => hsl_to_rgb((1.0 - t) * 250.0, 0.62, 0.52),
            // Trimmed of its darkest ends, which read as black on a molecule.
            Palette::Turbo => sample(
                &[
                    [0.239, 0.208, 0.547],
                    [0.275, 0.511, 0.973],
                    [0.122, 0.789, 0.866],
                    [0.234, 0.961, 0.556],
                    [0.633, 0.992, 0.239],
                    [0.906, 0.843, 0.222],
                    [0.995, 0.575, 0.164],
                    [0.887, 0.262, 0.038],
                    [0.631, 0.069, 0.004],
                ],
                t,
            ),
            Palette::Viridis => sample(
                &[
                    [0.267, 0.005, 0.329],
                    [0.229, 0.322, 0.545],
                    [0.128, 0.567, 0.551],
                    [0.369, 0.789, 0.383],
                    [0.993, 0.906, 0.144],
                ],
                t,
            ),
            Palette::Cividis => sample(
                &[
                    [0.000, 0.135, 0.305],
                    [0.253, 0.344, 0.454],
                    [0.478, 0.508, 0.472],
                    [0.719, 0.699, 0.406],
                    [0.995, 0.909, 0.218],
                ],
                t,
            ),
            Palette::Magma => sample(
                &[
                    [0.001, 0.000, 0.014],
                    [0.316, 0.072, 0.485],
                    [0.716, 0.215, 0.475],
                    [0.986, 0.536, 0.382],
                    [0.987, 0.991, 0.749],
                ],
                t,
            ),
            Palette::Inferno => sample(
                &[
                    [0.001, 0.000, 0.014],
                    [0.342, 0.062, 0.429],
                    [0.736, 0.216, 0.330],
                    [0.978, 0.558, 0.035],
                    [0.988, 0.998, 0.645],
                ],
                t,
            ),
            Palette::Plasma => sample(
                &[
                    [0.050, 0.030, 0.528],
                    [0.494, 0.012, 0.658],
                    [0.798, 0.280, 0.470],
                    [0.973, 0.585, 0.254],
                    [0.940, 0.975, 0.131],
                ],
                t,
            ),
            Palette::Mako => sample(
                &[
                    [0.043, 0.035, 0.052],
                    [0.238, 0.235, 0.443],
                    [0.263, 0.478, 0.635],
                    [0.362, 0.716, 0.663],
                    [0.867, 0.949, 0.792],
                ],
                t,
            ),
            Palette::Rocket => sample(
                &[
                    [0.011, 0.018, 0.100],
                    [0.381, 0.121, 0.325],
                    [0.797, 0.105, 0.311],
                    [0.959, 0.533, 0.375],
                    [0.981, 0.921, 0.865],
                ],
                t,
            ),
            Palette::Ember => sample(
                &[
                    [0.20, 0.07, 0.32],
                    [0.55, 0.15, 0.40],
                    [0.83, 0.33, 0.32],
                    [0.96, 0.60, 0.25],
                    [0.99, 0.87, 0.44],
                ],
                t,
            ),
            Palette::Ice => sample(
                &[
                    [0.07, 0.12, 0.34],
                    [0.11, 0.37, 0.62],
                    [0.24, 0.64, 0.75],
                    [0.72, 0.92, 0.92],
                ],
                t,
            ),
            Palette::Forest => sample(
                &[
                    [0.07, 0.19, 0.14],
                    [0.13, 0.38, 0.24],
                    [0.40, 0.58, 0.27],
                    [0.73, 0.76, 0.40],
                    [0.94, 0.92, 0.70],
                ],
                t,
            ),
            Palette::Dusk => sample(
                &[
                    [0.12, 0.13, 0.33],
                    [0.36, 0.22, 0.50],
                    [0.67, 0.32, 0.55],
                    [0.92, 0.52, 0.52],
                    [0.99, 0.80, 0.62],
                ],
                t,
            ),
            Palette::Pastel => sample(
                &[
                    [0.55, 0.66, 0.90],
                    [0.55, 0.82, 0.82],
                    [0.68, 0.86, 0.60],
                    [0.97, 0.87, 0.57],
                    [0.97, 0.69, 0.58],
                    [0.91, 0.63, 0.76],
                ],
                t,
            ),
            Palette::Greyscale => sample(&[[0.16, 0.16, 0.17], [0.90, 0.90, 0.89]], t),
        }
    }

    /// A distinct colour per chain. Palettes with a set of their own use it:
    /// Spectrum and Pastel keep hand-picked hues, Viridis borrows the
    /// Okabe–Ito set, built to stay separable under colour vision deficiency,
    /// and Forest and Dusk lift theirs off their ramps with the hues pulled
    /// apart. The rest step along their own ramp, interleaved so neighbouring
    /// chains do not land next to each other.
    pub fn chain(self, i: usize) -> [f32; 3] {
        const SPECTRUM: [[f32; 3]; 8] = [
            [0.29, 0.47, 0.78],
            [0.85, 0.44, 0.31],
            [0.36, 0.66, 0.48],
            [0.72, 0.52, 0.80],
            [0.88, 0.70, 0.26],
            [0.36, 0.70, 0.76],
            [0.83, 0.47, 0.62],
            [0.52, 0.55, 0.62],
        ];
        const OKABE_ITO: [[f32; 3]; 8] = [
            [0.000, 0.447, 0.698],
            [0.902, 0.624, 0.000],
            [0.000, 0.620, 0.451],
            [0.800, 0.475, 0.655],
            [0.337, 0.706, 0.914],
            [0.835, 0.369, 0.000],
            [0.941, 0.894, 0.259],
            [0.350, 0.350, 0.350],
        ];
        const PASTEL: [[f32; 3]; 8] = [
            [0.96, 0.62, 0.55],
            [0.55, 0.72, 0.92],
            [0.70, 0.86, 0.58],
            [0.98, 0.85, 0.55],
            [0.78, 0.66, 0.90],
            [0.56, 0.84, 0.82],
            [0.95, 0.70, 0.82],
            [0.80, 0.74, 0.64],
        ];
        const FOREST: [[f32; 3]; 8] = [
            [0.16, 0.42, 0.26],
            [0.80, 0.70, 0.36],
            [0.42, 0.62, 0.28],
            [0.55, 0.40, 0.24],
            [0.30, 0.56, 0.52],
            [0.72, 0.78, 0.48],
            [0.10, 0.26, 0.20],
            [0.62, 0.56, 0.42],
        ];
        const DUSK: [[f32; 3]; 8] = [
            [0.36, 0.24, 0.56],
            [0.95, 0.60, 0.50],
            [0.68, 0.30, 0.52],
            [0.98, 0.82, 0.62],
            [0.20, 0.22, 0.44],
            [0.86, 0.44, 0.52],
            [0.54, 0.46, 0.72],
            [0.76, 0.58, 0.62],
        ];
        const ORDER: [usize; 8] = [0, 4, 2, 6, 1, 5, 3, 7];
        match self {
            Palette::Spectrum => SPECTRUM[i % 8],
            Palette::Viridis => OKABE_ITO[i % 8],
            Palette::Pastel => PASTEL[i % 8],
            Palette::Forest => FOREST[i % 8],
            Palette::Dusk => DUSK[i % 8],
            _ => self.ramp(ORDER[i % 8] as f32 / 7.0),
        }
    }
}

/// Even stops along a ramp table.
fn sample(stops: &[[f32; 3]], t: f32) -> [f32; 3] {
    let last = stops.len() - 1;
    let x = t.clamp(0.0, 1.0) * last as f32;
    let i = (x.floor() as usize).min(last - 1);
    mix(stops[i], stops[i + 1], x - i as f32)
}

/// How the frame is shaded.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Style {
    /// Smooth shading on paper, with distance washing out.
    Soft,
    /// Ink outlines, vivid colour, banded light: the look of a drawn figure.
    Illustrative,
    /// Soft colour laid on almost flat, held by ink outlines: the look of a
    /// painted molecular illustration, after David Goodsell.
    Flat,
    /// Glossy light and shadow where atoms crowd together, without ink: the
    /// look of a rendered figure for a paper or a cover.
    Studio,
}

impl Style {
    /// The drawn figure leads: it is what the viewer opens with.
    pub const ALL: [Style; 4] = [Style::Illustrative, Style::Flat, Style::Soft, Style::Studio];
    pub fn label(self) -> &'static str {
        match self {
            Style::Soft => "Soft",
            Style::Illustrative => "Illustrative",
            Style::Flat => "Flat",
            Style::Studio => "Studio",
        }
    }

    /// Whether the edges are drawn in ink.
    pub fn inked(self) -> bool {
        matches!(self, Style::Illustrative | Style::Flat)
    }

    /// How far the farthest geometry washes toward the paper. A drawn figure
    /// does not fade with distance, so its outlines carry the depth instead;
    /// a painted one fades a little, as the back of a crowded scene does.
    ///
    /// On a dark page the wash is toward black, which the back of a dark
    /// palette is already close to, so it is held to half.
    pub fn fog(self) -> f32 {
        let fog = match self {
            Style::Soft => 0.26,
            Style::Illustrative => 0.0,
            Style::Flat => 0.2,
            Style::Studio => 0.1,
        };
        if crate::theme::dark() { fog * 0.5 } else { fog }
    }

    /// The style as the GPU shader reads it.
    pub fn code(self) -> f32 {
        match self {
            Style::Soft => 0.0,
            Style::Illustrative => 1.0,
            Style::Flat => 2.0,
            Style::Studio => 3.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Camera {
    /// World to view, as rows. Dragging turns it about the screen's own
    /// axes, trackball fashion, so the molecule can go all the way round in
    /// any direction with nothing to get stuck on.
    pub orientation: [[f32; 3]; 3],
    /// Multiplier on the auto-fit scale.
    pub zoom: f32,
    pub pan: [f32; 2],
    pub center: [f32; 3],
    pub radius: f32,
    /// Near things drawn larger than far ones, as a lens sees them, rather
    /// than all at one scale.
    pub perspective: bool,
    /// Only what lies between these two depths is drawn, near then far, in
    /// view space about the centre: a slice through the structure.
    pub slab: Option<[f32; 2]>,
}

impl Camera {
    /// How far the eye stands from the centre in perspective, in radii of
    /// the structure framed.
    pub const EYE: f32 = 3.0;

    /// A point in view space as drawn: in perspective, moved out from the
    /// middle of the view by how near it is, and with it the factor sizes
    /// there scale by. Each sphere or rod keeps its round section, scaled
    /// whole, as impostors drawn in perspective commonly are.
    pub fn warp(&self, v: [f32; 3]) -> ([f32; 3], f32) {
        if !self.perspective {
            return (v, 1.0);
        }
        let k = 1.0 / (1.0 - v[2] / (Self::EYE * self.radius)).max(0.2);
        ([v[0] * k, v[1] * k, v[2]], k)
    }

    /// The nearest and furthest depths drawn.
    pub fn depths(&self) -> (f32, f32) {
        match self.slab {
            Some([near, far]) => (near, far),
            None => (f32::INFINITY, f32::NEG_INFINITY),
        }
    }

    /// The three-quarter view every structure opens at, which reads as
    /// zero degrees on the angle readout.
    pub const HOME_YAW: f32 = 0.55;
    pub const HOME_PITCH: f32 = 0.25;

    /// The opening orientation.
    pub fn home() -> [[f32; 3]; 3] {
        mul(tilt(Self::HOME_PITCH), spin(Self::HOME_YAW))
    }

    /// Turn the view about the screen's vertical axis (`across`, as a drag
    /// from side to side does) and its horizontal axis (`down`).
    pub fn rotate(&mut self, across: f32, down: f32) {
        let turned = mul(mul(tilt(down), spin(across)), self.orientation);
        self.orientation = orthonormal(turned);
    }

    /// Turn the view in the screen's plane, clockwise as seen for positive
    /// `a`.
    pub fn roll(&mut self, a: f32) {
        let (s, c) = a.sin_cos();
        let roll = [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]];
        self.orientation = orthonormal(mul(roll, self.orientation));
    }

    /// Rotation away from the opening view, in degrees about the screen's
    /// x axis (tilting) and y axis (turning), each in -180..=180.
    pub fn angles(&self) -> [f32; 2] {
        let wrap = |r: f32| {
            let d = r.to_degrees().rem_euclid(360.0);
            if d > 180.0 {
                d - 360.0
            } else {
                d
            }
        };
        // The turn from home, read as a tilt after a spin: tilt(a) * spin(b)
        // has -sin b and cos b along its first row, and -sin a and cos a
        // down its middle column.
        let r = mul(self.orientation, transpose(Self::home()));
        let turn = (-r[0][2]).atan2(r[0][0]);
        let tilt = (-r[2][1]).atan2(r[1][1]);
        [wrap(tilt), wrap(turn)]
    }

    /// Set the view to read `[x, y]` degrees on the angle readout.
    pub fn set_angles(&mut self, [x, y]: [f32; 2]) {
        self.orientation = mul(mul(tilt(x.to_radians()), spin(y.to_radians())), Self::home());
    }

    pub fn frame(structure: &Structure) -> Self {
        Self::frame_all(&[structure])
    }

    /// Frame everything on screen, so a superimposed pair stays in view even
    /// where the prediction wanders away from the reference.
    pub fn frame_all(structures: &[&Structure]) -> Self {
        let mut lo = [f32::MAX; 3];
        let mut hi = [f32::MIN; 3];
        // The cell's corners count too, so a crystal is framed whole. A
        // protein's crystal cell is not drawn, and does not count.
        let points = |s: &Structure| {
            let polymer = s.has_polymer();
            s.atoms
                .iter()
                .map(|a| a.pos)
                .chain(s.cell.iter().filter(|_| !polymer).flat_map(|c| c.corners()))
                .collect::<Vec<_>>()
        };
        for s in structures {
            for p in points(s) {
                for i in 0..3 {
                    lo[i] = lo[i].min(p[i]);
                    hi[i] = hi[i].max(p[i]);
                }
            }
        }
        if lo[0] > hi[0] {
            return Camera {
                orientation: Camera::home(),
                zoom: 1.0,
                pan: [0.0, 0.0],
                center: [0.0; 3],
                radius: 1.0,
                perspective: false,
                slab: None,
            };
        }
        let center = [
            (lo[0] + hi[0]) * 0.5,
            (lo[1] + hi[1]) * 0.5,
            (lo[2] + hi[2]) * 0.5,
        ];
        let mut r2: f32 = 0.0;
        for s in structures {
            for p in points(s) {
                let d = [
                    p[0] - center[0],
                    p[1] - center[1],
                    p[2] - center[2],
                ];
                r2 = r2.max(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            }
        }
        let radius = r2.sqrt().max(1.0);
        Camera {
            orientation: Camera::home(),
            zoom: 1.0,
            pan: [0.0, 0.0],
            center,
            radius,
            perspective: false,
            slab: None,
        }
    }

    /// The view turned by `angle` radians about `axis`, an axis of the
    /// view itself rather than of the structure.
    pub fn turned(&self, axis: [f32; 3], angle: f32) -> [[f32; 3]; 3] {
        let a = unit(axis);
        if dot(a, a) < 0.5 || angle == 0.0 {
            return self.orientation;
        }
        // Rodrigues' rotation, as rows.
        let (s, c) = angle.sin_cos();
        let t = 1.0 - c;
        let r = [
            [t * a[0] * a[0] + c, t * a[0] * a[1] - s * a[2], t * a[0] * a[2] + s * a[1]],
            [t * a[0] * a[1] + s * a[2], t * a[1] * a[1] + c, t * a[1] * a[2] - s * a[0]],
            [t * a[0] * a[2] - s * a[1], t * a[1] * a[2] + s * a[0], t * a[2] * a[2] + c],
        ];
        orthonormal(mul(r, self.orientation))
    }

    /// The turn from this view to `target`'s: its axis in view space and
    /// its angle, such that `turned` by them gives `target`.
    pub fn turn_to(&self, target: [[f32; 3]; 3]) -> ([f32; 3], f32) {
        // R = target · orientationᵀ, as a rotation of view space.
        let o = self.orientation;
        let r: [[f32; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| (0..3).map(|k| target[i][k] * o[j][k]).sum()));
        let cos = ((r[0][0] + r[1][1] + r[2][2] - 1.0) * 0.5).clamp(-1.0, 1.0);
        let angle = cos.acos();
        if angle < 1e-5 {
            return ([0.0, 0.0, 1.0], 0.0);
        }
        let axis = [r[2][1] - r[1][2], r[0][2] - r[2][0], r[1][0] - r[0][1]];
        if dot(axis, axis) > 1e-8 {
            return (unit(axis), angle);
        }
        // Half a turn: the axis is the column of R + I that is not nothing.
        let k = (0..3).max_by(|&a, &b| r[a][a].total_cmp(&r[b][b])).unwrap();
        let col = [r[0][k] + if k == 0 { 1.0 } else { 0.0 }, r[1][k] + if k == 1 { 1.0 } else { 0.0 }, r[2][k] + if k == 2 { 1.0 } else { 0.0 }];
        (unit(col), angle)
    }

    /// The turn that brings world direction `d` round to point at the
    /// viewer: its axis in view space, and its angle.
    pub fn turn_to_face(&self, d: [f32; 3]) -> ([f32; 3], f32) {
        let m = self.orientation;
        let v = unit([dot(m[0], d), dot(m[1], d), dot(m[2], d)]);
        let toward = [0.0, 0.0, 1.0];
        (cross(v, toward), dot(v, toward).clamp(-1.0, 1.0).acos())
    }

    pub(crate) fn basis(&self) -> [[f32; 3]; 3] {
        self.orientation
    }
}

/// Turning about the view's vertical axis, as rows: world to view.
fn spin(a: f32) -> [[f32; 3]; 3] {
    let (s, c) = a.sin_cos();
    [[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]]
}

/// Tilting about the view's horizontal axis, as rows: world to view.
fn tilt(a: f32) -> [[f32; 3]; 3] {
    let (s, c) = a.sin_cos();
    [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]]
}

pub(crate) fn mul(a: [[f32; 3]; 3], b: [[f32; 3]; 3]) -> [[f32; 3]; 3] {
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    m
}

fn transpose(a: [[f32; 3]; 3]) -> [[f32; 3]; 3] {
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = a[j][i];
        }
    }
    m
}

/// Square the rows up again, so thousands of small turns do not let the
/// rotation drift into a skew.
pub(crate) fn orthonormal(m: [[f32; 3]; 3]) -> [[f32; 3]; 3] {
    let dot = |a: [f32; 3], b: [f32; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let unit = |a: [f32; 3]| {
        let n = dot(a, a).sqrt().max(1e-12);
        [a[0] / n, a[1] / n, a[2] / n]
    };
    let x = unit(m[0]);
    let d = dot(m[1], x);
    let y = unit([m[1][0] - d * x[0], m[1][1] - d * x[1], m[1][2] - d * x[2]]);
    let z = [
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    ];
    [x, y, z]
}

#[derive(Clone, Copy)]
pub struct Sphere {
    pub c: [f32; 3],
    pub r: f32,
    pub color: [f32; 3],
    /// Opacity; anything under one is blended after the opaque scene.
    pub alpha: f32,
}

/// A sphere-swept segment. The tube is built from these rather than from a
/// chain of spheres: consecutive capsules share an end cap exactly, so the
/// surface has no seams to shade.
#[derive(Clone, Copy)]
pub struct Capsule {
    pub a: [f32; 3],
    pub b: [f32; 3],
    pub r: f32,
    pub ca: [f32; 3],
    pub cb: [f32; 3],
    pub alpha: f32,
}

/// A shaded triangle with per-vertex normals, used for cartoon ribbons.
#[derive(Clone, Copy)]
pub struct Tri {
    pub p: [[f32; 3]; 3],
    pub n: [[f32; 3]; 3],
    pub c: [f32; 3],
}

#[derive(Clone, Default)]
pub struct Scene {
    pub spheres: Vec<Sphere>,
    pub capsules: Vec<Capsule>,
    pub tris: Vec<Tri>,
    /// The residue each ribbon triangle and tube capsule is part of, beside
    /// them rather than in them so the drawing never carries it. `NO_RESIDUE`
    /// for anything else; shorter than its list when the end has none.
    pub tri_residue: Vec<u32>,
    pub capsule_residue: Vec<u32>,
    /// Which build this is, for telling a frame of it from a frame of any
    /// other. Set by `stamp`; zero until then.
    pub version: u64,
}

/// Marks a triangle or capsule that is no residue's.
pub const NO_RESIDUE: u32 = u32::MAX;

impl Scene {
    /// Pad the residue lists out to the primitives they name.
    fn align(&mut self) {
        self.tri_residue.resize(self.tris.len(), NO_RESIDUE);
        self.capsule_residue.resize(self.capsules.len(), NO_RESIDUE);
    }

    /// Paint residue `residue`'s stretch of ribbon or tube, among triangles
    /// `tris` and capsules `capsules`, in the highlight's colour, or with
    /// `faint` most of the way to it.
    pub fn paint_residue(
        &mut self,
        tris: std::ops::Range<usize>,
        capsules: std::ops::Range<usize>,
        residue: u32,
        style: Style,
        faint: bool,
    ) {
        let white = highlight_color(style);
        let paint = |c: [f32; 3]| if faint { mix(c, white, 0.7) } else { white };
        for i in tris {
            if self.tri_residue.get(i) == Some(&residue) {
                self.tris[i].c = paint(self.tris[i].c);
            }
        }
        for i in capsules {
            if self.capsule_residue.get(i) == Some(&residue) {
                let c = &mut self.capsules[i];
                (c.ca, c.cb) = (paint(c.ca), paint(c.cb));
            }
        }
    }
    /// Give the scene a number no other scene has had, so a viewer caching
    /// frames never mistakes a recoloured scene for the one before it.
    pub fn stamp(mut self) -> Self {
        static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);
        self.version = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        self
    }

    pub fn len(&self) -> usize {
        self.spheres.len() + self.capsules.len() + self.tris.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn merge(&mut self, mut other: Scene) {
        self.align();
        other.align();
        self.tri_residue.append(&mut other.tri_residue);
        self.capsule_residue.append(&mut other.capsule_residue);
        self.spheres.append(&mut other.spheres);
        self.capsules.append(&mut other.capsules);
        self.tris.append(&mut other.tris);
    }
}

// ---------------------------------------------------------------------------
// Scene construction
// ---------------------------------------------------------------------------

/// The AlphaFold pLDDT bands, which people read at a glance.
fn plddt_color(v: f32) -> [f32; 3] {
    if v >= 90.0 {
        [0.05, 0.34, 0.77]
    } else if v >= 70.0 {
        [0.40, 0.79, 0.94]
    } else if v >= 50.0 {
        [1.00, 0.86, 0.36]
    } else {
        [0.98, 0.49, 0.24]
    }
}

fn rgb_to_hsl(c: [f32; 3]) -> (f32, f32, f32) {
    let (r, g, b) = (c[0], c[1], c[2]);
    let max = r.max(g).max(b);
    let min = r.min(g).min(b);
    let l = (max + min) * 0.5;
    if (max - min).abs() < 1e-6 {
        return (0.0, 0.0, l);
    }
    let d = max - min;
    let s = if l > 0.5 {
        d / (2.0 - max - min)
    } else {
        d / (max + min)
    };
    let h = if max == r {
        ((g - b) / d + if g < b { 6.0 } else { 0.0 }) * 60.0
    } else if max == g {
        ((b - r) / d + 2.0) * 60.0
    } else {
        ((r - g) / d + 4.0) * 60.0
    };
    (h, s, l)
}

/// Push a colour toward poster ink: saturated, and away from the extremes of
/// lightness so it still holds an outline.
fn vivid(c: [f32; 3]) -> [f32; 3] {
    let (h, s, l) = rgb_to_hsl(c);
    if s < 0.02 {
        return [l, l, l];
    }
    hsl_to_rgb(h, (s * 1.5 + 0.18).min(1.0), 0.5 + (l - 0.5) * 0.62)
}

/// Deepen a colour a little, so it holds its own once lit and fogged.
fn rich(c: [f32; 3]) -> [f32; 3] {
    let (h, s, l) = rgb_to_hsl(c);
    if s < 0.02 {
        return c;
    }
    hsl_to_rgb(h, (s * 1.2 + 0.05).min(0.9), l)
}

/// Settle a colour into paint: full-bodied but never neon, and drawn toward
/// a mid lightness so the dark ink outline always stands off it.
fn chalky(c: [f32; 3]) -> [f32; 3] {
    let (h, s, l) = rgb_to_hsl(c);
    let l = 0.6 + (l - 0.5) * 0.45;
    if s < 0.02 {
        return [l, l, l];
    }
    hsl_to_rgb(h, (s * 0.95 + 0.04).min(0.7), l)
}

fn hsl_to_rgb(h: f32, s: f32, l: f32) -> [f32; 3] {
    let c = (1.0 - (2.0 * l - 1.0).abs()) * s;
    let hp = h / 60.0;
    let x = c * (1.0 - (hp % 2.0 - 1.0).abs());
    let (r, g, b) = match hp as i32 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };
    let m = l - c * 0.5;
    [r + m, g + m, b + m]
}

/// Deviation in angstrom, along the chosen ramp. Compressed at the low end,
/// where a tenth of an angstrom still means something, and saturating by 8 A,
/// past which the two have simply gone their separate ways.
fn deviation_color(d: f32, palette: Palette) -> [f32; 3] {
    let t = (d.max(0.0) / 8.0).powf(0.6).clamp(0.0, 1.0);
    palette.ramp(t)
}

fn residue_color(
    s: &Structure,
    residue: usize,
    coloring: Coloring,
    palette: Palette,
) -> [f32; 3] {
    let r = &s.residues[residue];
    match coloring {
        Coloring::Pair => s.tint.unwrap_or(s.chain_color(r.chain as usize, palette)),
        Coloring::Deviation => match s.deviation.get(residue) {
            Some(d) if *d >= 0.0 => deviation_color(*d, palette),
            // No counterpart to compare against.
            _ => [0.74, 0.73, 0.71],
        },
        Coloring::Chain => s.chain_color(r.chain as usize, palette),
        Coloring::Rainbow => palette.ramp(s.along(residue)),
        // pLDDT keeps AlphaFold's own bands: they are a published convention
        // and a reader knows them by sight.
        Coloring::Confidence => {
            if s.has_confidence {
                plddt_color(r.b)
            } else {
                s.chain_color(r.chain as usize, palette)
            }
        }
        Coloring::BFactor => match s.b_span() {
            Some((lo, hi)) => palette.ramp(((r.b - lo) / (hi - lo)).clamp(0.0, 1.0)),
            None => s.chain_color(r.chain as usize, palette),
        },
        Coloring::Structure => ss_color(r.ss),
        Coloring::Hydrophobic => match kyte_doolittle(r.code) {
            Some(h) => diverging(h / 4.5, HYDROPHILIC, HYDROPHOBIC),
            None => NEUTRAL,
        },
        Coloring::Charge => match r.name.as_str() {
            "ASP" | "GLU" => NEGATIVE,
            "LYS" | "ARG" => POSITIVE,
            _ => NEUTRAL,
        },
        Coloring::Element => [0.42, 0.44, 0.47],
    }
}

/// The ends of the hydrophobicity and charge scales, and their middle:
/// fixed, as a reader knows them by sight, rather than the palette's.
pub const HYDROPHILIC: [f32; 3] = [0.25, 0.53, 0.80];
pub const HYDROPHOBIC: [f32; 3] = [0.88, 0.50, 0.16];
pub const NEGATIVE: [f32; 3] = [0.86, 0.22, 0.22];
pub const POSITIVE: [f32; 3] = [0.22, 0.40, 0.88];
pub const NEUTRAL: [f32; 3] = [0.94, 0.93, 0.91];

/// From `low` at -1 through white at 0 to `high` at 1.
pub fn diverging(t: f32, low: [f32; 3], high: [f32; 3]) -> [f32; 3] {
    let t = t.clamp(-1.0, 1.0);
    match t {
        -1.0 => low,
        1.0 => high,
        t if t < 0.0 => mix(NEUTRAL, low, -t),
        t => mix(NEUTRAL, high, t),
    }
}

/// Kyte and Doolittle's hydropathy, from -4.5 (arginine) to 4.5 (isoleucine).
pub fn kyte_doolittle(code: char) -> Option<f32> {
    Some(match code {
        'I' => 4.5,
        'V' => 4.2,
        'L' => 3.8,
        'F' => 2.8,
        'C' => 2.5,
        'M' => 1.9,
        'A' => 1.8,
        'G' => -0.4,
        'T' => -0.7,
        'S' => -0.8,
        'W' => -0.9,
        'Y' => -1.3,
        'P' => -1.6,
        'H' => -3.2,
        'E' | 'Q' | 'D' | 'N' => -3.5,
        'K' => -3.9,
        'R' => -4.5,
        _ => return None,
    })
}

/// The charge an atom carries at neutral pH: the charged groups of Asp,
/// Glu, Lys and Arg spread over their atoms, the chain's C terminus, a
/// nucleic acid's phosphates, and metal ions.
fn charge_of(s: &Structure, i: usize) -> f32 {
    let a = &s.atoms[i];
    let r = &s.residues[a.residue as usize];
    match (r.name.as_str(), a.name()) {
        ("ASP", "OD1" | "OD2") | ("GLU", "OE1" | "OE2") => -0.5,
        ("LYS", "NZ") => 1.0,
        ("ARG", "NH1" | "NH2" | "NE") => 1.0 / 3.0,
        (_, "OXT") => -1.0,
        // A strand's phosphates, one charge each.
        (_, "OP1" | "OP2" | "O1P" | "O2P") if r.nucleic.is_some() => -0.5,
        _ if r.outside_polymer() && matches!(a.element.0, 12 | 20 | 25 | 26 | 27 | 28 | 29 | 30) => 2.0,
        _ if r.outside_polymer() && matches!(a.element.0, 3 | 11 | 19) => 1.0,
        _ => 0.0,
    }
}

/// How far the potential reaches, in ångström: the Debye length of
/// physiological salt, about 8 Å, and twice that as the cut-off.
const DEBYE: f32 = 8.0;

/// The electrostatic potential at each point, from the charges of `s`, as
/// Coulomb's law screened by salt and damped by a dielectric rising with
/// distance: the usual quick picture, not a Poisson-Boltzmann solution.
fn potentials(s: &Structure, points: &[[f32; 3]]) -> Vec<f32> {
    use rayon::prelude::*;
    let charges: Vec<(usize, f32)> = (0..s.atoms.len()).map(|i| (i, charge_of(s, i))).filter(|(_, q)| *q != 0.0).collect();
    let cell = 2.0 * DEBYE;
    let key = |p: [f32; 3]| [0, 1, 2].map(|k| (p[k] / cell).floor() as i32);
    let mut grid: HashMap<[i32; 3], Vec<(usize, f32)>> = HashMap::new();
    for &(i, q) in &charges {
        grid.entry(key(s.atoms[i].pos)).or_default().push((i, q));
    }
    points
        .par_iter()
        .map(|&p| {
            let c = key(p);
            let mut phi = 0.0;
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        for &(i, q) in grid.get(&[c[0] + dx, c[1] + dy, c[2] + dz]).into_iter().flatten() {
                            let d = sub(p, s.atoms[i].pos);
                            let r = dot(d, d).sqrt().max(1.0);
                            if r < cell {
                                phi += q * (-r / DEBYE).exp() / (r * r);
                            }
                        }
                    }
                }
            }
            phi
        })
        .collect()
}

/// The potential at which the surface is fully red or blue.
const PHI_FULL: f32 = 0.06;

/// The colours of helix, strand and coil: kept apart from the palette, as
/// pLDDT's bands are, so they read the same whatever it is.
pub const SS_COLORS: [(Ss, &str, [f32; 3]); 3] = [
    (Ss::Helix, "Helix", [0.86, 0.26, 0.46]),
    (Ss::Sheet, "Strand", [0.98, 0.76, 0.20]),
    (Ss::Coil, "Coil", [0.78, 0.78, 0.76]),
];

fn ss_color(ss: Ss) -> [f32; 3] {
    SS_COLORS.iter().find(|(s, _, _)| *s == ss).map_or([0.78; 3], |(_, _, c)| *c)
}

/// pLDDT's bands, from the top down, as AlphaFold DB shows them.
pub const PLDDT_BANDS: [(&str, f32); 4] = [("Very high", 95.0), ("Confident", 80.0), ("Low", 60.0), ("Very low", 30.0)];

pub fn plddt_band_color(v: f32) -> [f32; 3] {
    plddt_color(v)
}

fn catmull_rom(p0: [f32; 3], p1: [f32; 3], p2: [f32; 3], p3: [f32; 3], t: f32) -> [f32; 3] {
    let t2 = t * t;
    let t3 = t2 * t;
    let mut out = [0.0; 3];
    for i in 0..3 {
        out[i] = 0.5
            * ((2.0 * p1[i])
                + (-p0[i] + p2[i]) * t
                + (2.0 * p0[i] - 5.0 * p1[i] + 4.0 * p2[i] - p3[i]) * t2
                + (-p0[i] + 3.0 * p1[i] - 3.0 * p2[i] + p3[i]) * t3);
    }
    out
}

/// Colours for one structure under one scheme, with what the scheme needs to
/// know about the structure worked out once.
/// The colour of each element: the chemists' convention, CPK by way of
/// Jmol, unless `standard` is off, when the elements present take the
/// palette's colours, most common first.
pub fn element_colors(s: &Structure, palette: Palette, standard: bool) -> Vec<[f32; 3]> {
    let mut out: Vec<[f32; 3]> = (0..=255u8)
        .map(|z| crate::elements::Element(z).color())
        .collect();
    if !standard {
        for (rank, (e, _)) in s
            .composition()
            .into_iter()
            .filter(|(e, _)| e.is_known())
            .enumerate()
        {
            out[e.0 as usize] = palette.chain(rank);
        }
    }
    out
}

struct Paint<'a> {
    s: &'a Structure,
    /// By atomic number.
    elements: Vec<[f32; 3]>,
    coloring: Coloring,
    palette: Palette,
    polymer: bool,
    /// Size and opacity, per atom.
    looks: Vec<Look>,
    /// How many atoms are the structure's own; past them come the boundary
    /// copies of a crystal, each repeating the atom `images` names.
    own: usize,
    images: Vec<u32>,
    /// Atoms hidden from every style, the backbone's included.
    hidden: Option<&'a [bool]>,
}

impl Paint<'_> {
    fn residue(&self, residue: usize) -> [f32; 3] {
        residue_color(self.s, residue, self.coloring, self.palette)
    }

    /// The structure's own atom that atom `i` is, or is a copy of.
    fn source(&self, i: usize) -> usize {
        if i < self.own {
            i
        } else {
            self.images[i - self.own] as usize
        }
    }

    fn atom(&self, i: usize) -> [f32; 3] {
        let a = &self.s.atoms[i];
        let r = &self.s.residues[a.residue as usize];
        let ligand = self.polymer && r.outside_polymer();
        let comparing = matches!(self.coloring, Coloring::Pair | Coloring::Deviation);
        // A ligand given a colour of its own wears it under every scheme
        // but the comparison's.
        if ligand && !comparing {
            match self.s.ligand_paint.get(&r.name) {
                Some(LigandPaint::Solid(c)) => return *c,
                Some(LigandPaint::Carbons(c)) if a.element.0 == 6 => return *c,
                _ => {}
            }
        }
        // Beside a polymer, ligands and ions keep their elements' colours
        // under any scheme that colours the polymer itself.
        let by_element = match self.coloring {
            Coloring::Element => true,
            _ => ligand && !comparing,
        };
        if by_element {
            return match self.looks[i].color {
                Some(c) => c,
                None if a.element.is_known() => self.elements[a.element.0 as usize],
                // A simulation's anonymous types get a colour each.
                None => self.palette.chain(a.kind as usize),
            };
        }
        match self.coloring {
            Coloring::Rainbow if !self.polymer => self
                .palette
                .ramp(self.source(i) as f32 / (self.own.max(2) - 1) as f32),
            // Without chains, each molecule is told apart instead.
            Coloring::Chain if !self.polymer => self.palette.chain(a.residue as usize),
            _ => self.residue(a.residue as usize),
        }
    }

    /// Drawing radius of an atom drawn at its full size, times the size set
    /// for its species.
    fn radius(&self, i: usize) -> f32 {
        let e = self.s.atoms[i].element;
        let r = if e.is_known() {
            e.vdw()
        } else {
            self.s.unknown_radius
        };
        r * self.looks[i].size
    }

    fn alpha(&self, i: usize) -> f32 {
        self.looks[i].alpha
    }
}


/// Build the drawable primitives for a structure. `quality` in 0..=1 trades
/// tube smoothness for speed on very long chains.
#[allow(clippy::too_many_arguments)]
pub fn build_scene(
    s: &Structure,
    repr: Repr,
    coloring: Coloring,
    quality: f32,
    style: Style,
    palette: Palette,
    standard_colors: bool,
    cell: bool,
    appearance: &Appearance,
) -> Scene {
    // Colours are settled on the structure's own atoms, so showing the cell
    // and its boundary copies does not reshuffle them.
    let elements = element_colors(s, palette, standard_colors);
    let own = s.atoms.len();
    let filled = if cell { with_boundary(s) } else { None };
    let (s, images) = match &filled {
        Some((filled, images)) => (filled, images.clone()),
        None => (s, Vec::new()),
    };
    let paint = Paint {
        s,
        elements,
        coloring,
        palette,
        polymer: s.has_polymer(),
        looks: appearance.resolve(s),
        own,
        images,
        hidden: appearance.mask(s),
    };
    let mut scene = build_primitives(&paint, repr, quality);
    if cell {
        if let Some(c) = &s.cell {
            scene.capsules.extend(cell_edges(c));
        }
    }
    settle(&mut scene, style);
    scene
}

/// Settle a scene's colours for the style and the page: each style's tone,
/// and lifted off a dark page.
fn settle(scene: &mut Scene, style: Style) {
    let tone: fn([f32; 3]) -> [f32; 3] = match style {
        Style::Soft | Style::Studio => rich,
        Style::Illustrative => vivid,
        Style::Flat => chalky,
    };
    let dark = crate::theme::dark();
    for c in scene
        .spheres
        .iter_mut()
        .map(|p| &mut p.color)
        .chain(scene.capsules.iter_mut().flat_map(|c| [&mut c.ca, &mut c.cb]))
        .chain(scene.tris.iter_mut().map(|t| &mut t.c))
    {
        *c = tone(*c);
        if dark {
            *c = lifted(*c);
        }
    }
}

/// Raise the floor of a colour's lightness for a dark page, so the dark end
/// of a ramp (magma's black, ice's navy) still stands off the paper. The
/// order of lightness is kept, so a perceptual palette still reads in order.
fn lifted(c: [f32; 3]) -> [f32; 3] {
    let (h, s, l) = rgb_to_hsl(c);
    hsl_to_rgb(h, s, 0.28 + l * 0.72)
}

fn build_primitives(paint: &Paint, repr: Repr, quality: f32) -> Scene {
    let s = paint.s;
    let (coloring, palette) = (paint.coloring, paint.palette);
    match repr {
        Repr::Spacefill => Scene {
            spheres: (0..s.atoms.len())
                .filter(|&i| !paint.looks[i].hidden())
                .map(|i| Sphere {
                    c: s.atoms[i].pos,
                    r: paint.radius(i),
                    color: paint.atom(i),
                    alpha: paint.alpha(i),
                })
                .collect(),
            ..Scene::default()
        },
        Repr::Surface => build_surface(paint, quality),
        Repr::BallAndStick | Repr::Licorice => {
            let all: Vec<bool> = vec![true; s.atoms.len()];
            build_sticks(paint, &all, repr == Repr::Licorice)
        }
        Repr::Tube | Repr::Cartoon => {
            let mut scene = Scene::default();
            if repr == Repr::Tube {
                (scene.capsules, scene.capsule_residue) = build_tube(s, coloring, quality, palette, paint.hidden);
            } else {
                (scene.tris, scene.tri_residue) = build_cartoon(s, coloring, quality, palette, paint.hidden);
            }
            let (capsules, owners) = build_nucleic(s, coloring, quality, palette, paint.hidden);
            scene.capsule_residue.resize(scene.capsules.len(), NO_RESIDUE);
            scene.capsules.extend(capsules);
            scene.capsule_residue.extend(owners);
            // Ligands, ions and nucleic acids have no alpha carbons to draw a
            // backbone through, and would vanish; they are drawn as atoms.
            // So is water, when it is shown at all.
            let ligand: Vec<bool> = s
                .atoms
                .iter()
                .map(|a| {
                    let r = &s.residues[a.residue as usize];
                    r.ca.is_none() && r.nucleic.is_none()
                })
                .collect();
            if ligand.iter().any(|l| *l) {
                scene.merge(build_sticks(paint, &ligand, false));
            }
            scene
        }
    }
}

/// Whether atom `i` is under the skin in the surface style: a polymer's own
/// atoms beside a polymer, every atom of anything else, and never solvent.
fn under_skin(s: &Structure, polymer: bool, i: usize) -> bool {
    let r = &s.residues[s.atoms[i].residue as usize];
    !is_solvent(&r.name) && (!polymer || !r.outside_polymer())
}

/// The surface over the atoms under the skin, each patch in its atom's
/// colour; ligands beside a polymer stand in it as ball and stick.
fn build_surface(paint: &Paint, quality: f32) -> Scene {
    let s = paint.s;
    let polymer = paint.polymer;
    let skin: Vec<usize> = (0..s.atoms.len())
        .filter(|&i| under_skin(s, polymer, i) && !paint.looks[i].hidden())
        .collect();
    let points: Vec<[f32; 3]> = skin.iter().map(|&i| s.atoms[i].pos).collect();
    let radii: Vec<f32> = skin.iter().map(|&i| paint.radius(i)).collect();
    // A fixed share of the structure's width, so a molecule framed close
    // up is as smooth as a protein framed whole; past the size of a large
    // protein, the grid's own cap takes over.
    let (lo, hi) = points.iter().fold(([f32::MAX; 3], [f32::MIN; 3]), |(lo, hi), p| {
        ([0, 1, 2].map(|k| lo[k].min(p[k])), [0, 1, 2].map(|k| hi[k].max(p[k])))
    });
    let width = (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f32::max);
    let spacing = (width / 100.0).clamp(0.15, 1.2) / quality.max(0.5);
    let facets = crate::surface::gaussian(&points, &radii, spacing);
    let mut scene = Scene {
        tris: Vec::with_capacity(facets.len()),
        tri_residue: Vec::with_capacity(facets.len()),
        ..Scene::default()
    };
    // In charge, the skin shows the potential over it, a little out from
    // each triangle, where a probe of solvent would feel it.
    let phi = (paint.coloring == Coloring::Charge).then(|| {
        let at: Vec<[f32; 3]> = facets
            .iter()
            .map(|f| {
                let m = [0, 1, 2].map(|k| (f.p[0][k] + f.p[1][k] + f.p[2][k]) / 3.0);
                let n = unit([0, 1, 2].map(|k| f.n[0][k] + f.n[1][k] + f.n[2][k]));
                [m[0] + n[0] * 1.4, m[1] + n[1] * 1.4, m[2] + n[2] * 1.4]
            })
            .collect();
        potentials(s, &at)
    });
    for (k, f) in facets.into_iter().enumerate() {
        let atom = skin[f.atom as usize];
        let c = match &phi {
            Some(phi) => diverging(phi[k] / PHI_FULL, NEGATIVE, POSITIVE),
            None => paint.atom(atom),
        };
        scene.tris.push(Tri { p: f.p, n: f.n, c });
        scene.tri_residue.push(s.atoms[atom].residue);
    }
    if polymer {
        let ligand: Vec<bool> = (0..s.atoms.len()).map(|i| !under_skin(s, polymer, i)).collect();
        if ligand.iter().any(|l| *l) {
            scene.merge(build_sticks(paint, &ligand, false));
        }
    }
    scene
}

/// Atoms and bonds among the atoms marked in `keep`. Each bond is two
/// half-rods, each in its own atom's colour.
fn build_sticks(paint: &Paint, keep: &[bool], licorice: bool) -> Scene {
    let s = paint.s;
    // Rods are sized for ångström; a simulation in other units scales them
    // to its own atoms.
    let unit = |i: usize| {
        if s.atoms[i].element.is_known() {
            1.0
        } else {
            (s.unknown_radius / 0.8).min(1.0)
        }
    };
    let rod = if licorice { 0.22 } else { 0.14 };

    let mut bonded = vec![false; s.atoms.len()];
    let mut capsules = Vec::new();
    for &[a, b] in &s.bonds {
        let (a, b) = (a as usize, b as usize);
        if !keep[a] || !keep[b] || paint.looks[a].hidden() || paint.looks[b].hidden() {
            continue;
        }
        bonded[a] = true;
        bonded[b] = true;
        let pa = s.atoms[a].pos;
        let pb = s.atoms[b].pos;
        let mid = mix(pa, pb, 0.5);
        let r = rod * unit(a).min(unit(b));
        let (ca, cb) = (paint.atom(a), paint.atom(b));
        capsules.push(Capsule {
            a: pa,
            b: mid,
            r,
            ca,
            cb: ca,
            alpha: paint.alpha(a),
        });
        capsules.push(Capsule {
            a: mid,
            b: pb,
            r,
            ca: cb,
            cb,
            alpha: paint.alpha(b),
        });
    }

    let spheres = (0..s.atoms.len())
        .filter(|&i| keep[i] && !paint.looks[i].hidden())
        .filter_map(|i| {
            let r = if licorice {
                // The rod ends round off bonded atoms; only the lone ones,
                // ions and the like, need a sphere of their own.
                // Grown past one, a bonded atom swells out of its rods.
                if bonded[i] {
                    let size = paint.looks[i].size;
                    if size <= 1.0 {
                        return None;
                    }
                    return Some(Sphere {
                        c: s.atoms[i].pos,
                        r: rod * unit(i) * size,
                        color: paint.atom(i),
                        alpha: paint.alpha(i),
                    });
                }
                paint.radius(i) * 0.3
            } else {
                paint.radius(i) * 0.25
            };
            Some(Sphere {
                c: s.atoms[i].pos,
                r,
                color: paint.atom(i),
                alpha: paint.alpha(i),
            })
        })
        .collect();
    Scene {
        spheres,
        capsules,
        ..Scene::default()
    }
}

// ---------------------------------------------------------------------------
// Picking
// ---------------------------------------------------------------------------

/// What is under a point of the viewport.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Picked {
    /// The atom, or a residue's alpha carbon when `residue` is set.
    pub atom: usize,
    /// Whether a whole residue was hit, as in the backbone styles, where
    /// the ribbon stands for the residue rather than any one atom.
    pub residue: bool,
    /// View-space depth of the surface hit; larger is nearer.
    pub depth: f32,
}

/// Pointers this close to a thin rod or a small atom, in pixels, still
/// count as on it.
const PICK_SLACK: f32 = 3.0;

/// Backbone styles are picked against a rod of this radius through the
/// alpha carbons: wide enough to cover a ribbon seen flat.
const BACKBONE_PICK_RADIUS: f32 = 1.3;

/// The atom drawn nearest the viewer at pixel `(x, y)` of a `width` by
/// `height` frame, as `render` would draw `s` in `repr` from `cam`. Worked
/// out from the atoms rather than read back from the frame, so it costs
/// nothing while nobody points at anything.
#[allow(clippy::too_many_arguments)]
pub fn pick(
    s: &Structure,
    repr: Repr,
    appearance: &Appearance,
    cam: &Camera,
    width: usize,
    height: usize,
    x: f32,
    y: f32,
) -> Option<Picked> {
    if s.atoms.is_empty() || width == 0 || height == 0 {
        return None;
    }
    let m = cam.basis();
    let scale = (width.min(height) as f32) / (2.2 * cam.radius) * cam.zoom;
    let p = [
        (x - (width as f32 * 0.5 + cam.pan[0])) / scale,
        -(y - (height as f32 * 0.5 + cam.pan[1])) / scale,
    ];
    let slack = PICK_SLACK / scale;
    let to_view = |q: [f32; 3]| {
        let q = sub(q, cam.center);
        [dot(m[0], q), dot(m[1], q), dot(m[2], q)]
    };
    let paint = Paint {
        s,
        elements: Vec::new(),
        coloring: Coloring::Element,
        palette: Palette::Spectrum,
        polymer: s.has_polymer(),
        looks: appearance.resolve(s),
        own: s.atoms.len(),
        images: Vec::new(),
        hidden: appearance.mask(s),
    };

    let mut best: Option<Picked> = None;
    let mut offer = |depth: Option<f32>, atom: usize, residue: bool| {
        if let Some(depth) = depth {
            if best.is_none_or(|b| depth > b.depth) {
                best = Some(Picked { atom, residue, depth });
            }
        }
    };
    // As drawn: in perspective, and cut by the slab.
    let (near, far) = cam.depths();
    let within = |z: Option<f32>| z.filter(|z| *z <= near && *z >= far);
    let ball = |c: [f32; 3], r: f32| {
        let (v, k) = cam.warp(to_view(c));
        within(sphere_depth(v, (r * k).max(slack), p))
    };
    let rod = |a: [f32; 3], b: [f32; 3], r: f32| {
        let ((a, ka), (b, kb)) = (cam.warp(to_view(a)), cam.warp(to_view(b)));
        within(capsule_depth(a, b, (r * (ka + kb) * 0.5).max(slack), p))
    };

    // The atoms drawn as atoms, and how: as in `build_primitives`.
    let backbone = matches!(repr, Repr::Tube | Repr::Cartoon);
    let surface = repr == Repr::Surface;
    let polymer = paint.polymer;
    let drawn = |i: usize| {
        let r = &s.residues[s.atoms[i].residue as usize];
        !paint.looks[i].hidden()
            && if backbone {
                r.ca.is_none() && r.nucleic.is_none()
            } else if surface {
                !under_skin(s, polymer, i)
            } else {
                true
            }
    };
    // Under a skin, each atom stands for the patch over it, which is its
    // residue's in a polymer.
    if surface {
        for i in (0..s.atoms.len()).filter(|&i| under_skin(s, polymer, i) && !paint.looks[i].hidden()) {
            offer(ball(s.atoms[i].pos, paint.radius(i) + crate::surface::SWELL), i, polymer);
        }
    }
    match repr {
        Repr::Spacefill => {
            for i in (0..s.atoms.len()).filter(|&i| drawn(i)) {
                offer(ball(s.atoms[i].pos, paint.radius(i)), i, false);
            }
        }
        _ => {
            let licorice = repr == Repr::Licorice;
            let rod_r = if licorice { 0.22 } else { 0.14 };
            for i in (0..s.atoms.len()).filter(|&i| drawn(i)) {
                let r = paint.radius(i) * if licorice { 0.3 } else { 0.25 };
                offer(ball(s.atoms[i].pos, r.max(rod_r)), i, false);
            }
            // Each half of a bond belongs to the atom at its end.
            for &[a, b] in &s.bonds {
                let (a, b) = (a as usize, b as usize);
                if !drawn(a) || !drawn(b) {
                    continue;
                }
                let (pa, pb) = (s.atoms[a].pos, s.atoms[b].pos);
                let mid = mix(pa, pb, 0.5);
                offer(rod(pa, mid, rod_r), a, false);
                offer(rod(mid, pb, rod_r), b, false);
            }
        }
    }
    if backbone {
        for chain in &s.chains {
            // A strand is picked by its tube and its rungs.
            for run in nucleic_runs(s, chain, appearance.mask(s)) {
                for &ri in &run {
                    let at = s.residues[ri].nucleic.unwrap() as usize;
                    offer(ball(s.atoms[at].pos, NUCLEIC_RADIUS + 0.3), at, true);
                    if let Some(tip) = base_tip(s, ri) {
                        offer(rod(s.atoms[at].pos, s.atoms[tip].pos, RUNG_RADIUS + 0.3), at, true);
                    }
                }
                for w in run.windows(2) {
                    let (a, b) = (s.residues[w[0]].nucleic.unwrap() as usize, s.residues[w[1]].nucleic.unwrap() as usize);
                    let mid = mix(s.atoms[a].pos, s.atoms[b].pos, 0.5);
                    offer(rod(s.atoms[a].pos, mid, NUCLEIC_RADIUS + 0.3), a, true);
                    offer(rod(mid, s.atoms[b].pos, NUCLEIC_RADIUS + 0.3), b, true);
                }
            }
            for run in ca_runs(s, chain, appearance.mask(s)) {
                let ca = |k: usize| s.residues[run[k]].ca.unwrap() as usize;
                for k in 0..run.len() {
                    let here = s.atoms[ca(k)].pos;
                    offer(ball(here, BACKBONE_PICK_RADIUS), ca(k), true);
                    if k + 1 < run.len() {
                        let mid = mix(here, s.atoms[ca(k + 1)].pos, 0.5);
                        offer(rod(here, mid, BACKBONE_PICK_RADIUS), ca(k), true);
                        offer(rod(mid, s.atoms[ca(k + 1)].pos, BACKBONE_PICK_RADIUS), ca(k + 1), true);
                    }
                }
            }
        }
    }
    best
}

/// Depth of the near surface of a view-space sphere over point `p`.
fn sphere_depth(c: [f32; 3], r: f32, p: [f32; 2]) -> Option<f32> {
    sphere_hit(c, r, p[0], p[1]).map(|(z, _)| z)
}

/// Depth of the near surface of a view-space capsule over point `p`,
/// taking the closest point of its axis as seen on the screen.
fn capsule_depth(a: [f32; 3], b: [f32; 3], r: f32, p: [f32; 2]) -> Option<f32> {
    let d = [b[0] - a[0], b[1] - a[1]];
    let len2 = d[0] * d[0] + d[1] * d[1];
    let t = if len2 > 1e-12 {
        (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    sphere_depth(mix(a, b, t), r, p)
}

/// The colour a picked atom or residue glows with.
const HIGHLIGHT: [f32; 3] = [1.0, 1.0, 1.0];

/// The highlight's colour as a scene in `style` shows it.
fn highlight_color(style: Style) -> [f32; 3] {
    let mut one = Scene {
        tris: vec![Tri { p: [[0.0; 3]; 3], n: [[0.0; 3]; 3], c: HIGHLIGHT }],
        ..Scene::default()
    };
    settle(&mut one, style);
    one.tris[0].c
}

/// What marks a picked atom: a glow around it. A residue picked in an
/// atomistic style glows on each of its atoms; in a backbone style it is
/// painted on the ribbon instead, by `Scene::paint_residue`, and nothing is
/// built for it here. Nothing is drawn that the style does not draw already.
pub fn build_highlight(
    s: &Structure,
    picked: Picked,
    repr: Repr,
    style: Style,
    appearance: &Appearance,
) -> Scene {
    let Some(atom) = s.atoms.get(picked.atom) else {
        return Scene::default();
    };
    let paint = Paint {
        s,
        elements: Vec::new(),
        coloring: Coloring::Element,
        palette: Palette::Spectrum,
        polymer: s.has_polymer(),
        looks: appearance.resolve(s),
        own: s.atoms.len(),
        images: Vec::new(),
        hidden: appearance.mask(s),
    };
    // The radius an atom is drawn at in this style.
    let drawn = |i: usize| match repr {
        Repr::Spacefill => paint.radius(i),
        Repr::Surface => paint.radius(i) + crate::surface::SWELL,
        Repr::Licorice => (paint.radius(i) * 0.3).max(0.22),
        _ => paint.radius(i) * 0.25,
    };
    let glow = |c: [f32; 3], r: f32| Sphere {
        c,
        r,
        color: HIGHLIGHT,
        alpha: 0.5,
    };
    let backbone = matches!(repr, Repr::Tube | Repr::Cartoon);
    let residue = &s.residues[atom.residue as usize];
    let mut scene = Scene::default();
    let skin = repr == Repr::Surface && under_skin(s, paint.polymer, picked.atom);
    match (picked.residue, backbone, residue.ca.or(residue.nucleic)) {
        // Painted on the ribbon or the skin instead.
        (true, true, Some(_)) => {}
        (true, _, _) if skin => {}
        (true, _, _) => {
            scene.spheres.extend(
                (0..s.atoms.len())
                    .filter(|&i| s.atoms[i].residue == atom.residue && !paint.looks[i].hidden())
                    .map(|i| glow(s.atoms[i].pos, drawn(i) + 0.5)),
            );
        }
        (false, _, _) => {
            let i = picked.atom;
            scene.spheres.push(Sphere {
                alpha: 0.6,
                ..glow(atom.pos, drawn(i) + 0.7)
            });
        }
    }
    settle(&mut scene, style);
    scene
}


/// A binding site's side chains, as sticks in their elements' colours, from
/// the alpha carbon out. Only a backbone style needs them: the atomistic
/// styles draw every atom already.
#[allow(clippy::too_many_arguments)]
pub fn build_site(
    s: &Structure,
    residues: &[usize],
    repr: Repr,
    style: Style,
    palette: Palette,
    standard_colors: bool,
    appearance: &Appearance,
) -> Scene {
    if !matches!(repr, Repr::Tube | Repr::Cartoon) {
        return Scene::default();
    }
    let mut keep = vec![false; s.atoms.len()];
    let atoms_of = s.atoms_by_residue();
    for &r in residues {
        let Some(res) = s.residues.get(r) else { continue };
        let backbone = [res.n, res.c, res.o];
        for &i in atoms_of.of(r) {
            keep[i as usize] = !backbone.contains(&Some(i));
        }
    }
    let paint = Paint {
        s,
        elements: element_colors(s, palette, standard_colors),
        coloring: Coloring::Element,
        palette,
        polymer: s.has_polymer(),
        looks: appearance.resolve(s),
        own: s.atoms.len(),
        images: Vec::new(),
        hidden: appearance.mask(s),
    };
    let mut scene = build_sticks(&paint, &keep, false);
    settle(&mut scene, style);
    scene
}

/// How a selection's atoms are drawn, over whatever the structure is drawn as.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LayerStyle {
    Sticks,
    Spheres,
    /// Taken out of every style.
    Hidden,
}

impl LayerStyle {
    pub const ALL: [LayerStyle; 3] = [LayerStyle::Sticks, LayerStyle::Spheres, LayerStyle::Hidden];

    pub fn label(self) -> &'static str {
        match self {
            LayerStyle::Sticks => "Sticks",
            LayerStyle::Spheres => "Spheres",
            LayerStyle::Hidden => "Hide",
        }
    }
}

/// A selection's atoms drawn as sticks or spheres, in their elements'
/// colours or in `color` throughout.
#[allow(clippy::too_many_arguments)]
pub fn build_layer(
    s: &Structure,
    atoms: &[bool],
    layer: LayerStyle,
    color: Option<[f32; 3]>,
    style: Style,
    palette: Palette,
    standard_colors: bool,
    appearance: &Appearance,
) -> Scene {
    if atoms.len() != s.atoms.len() || layer == LayerStyle::Hidden {
        return Scene::default();
    }
    let paint = Paint {
        s,
        elements: element_colors(s, palette, standard_colors),
        coloring: Coloring::Element,
        palette,
        polymer: s.has_polymer(),
        looks: appearance.resolve(s),
        own: s.atoms.len(),
        images: Vec::new(),
        hidden: appearance.mask(s),
    };
    let mut scene = match layer {
        LayerStyle::Sticks => build_sticks(&paint, atoms, true),
        _ => Scene {
            spheres: (0..s.atoms.len())
                .filter(|&i| atoms[i] && !paint.looks[i].hidden())
                .map(|i| Sphere { c: s.atoms[i].pos, r: paint.radius(i), color: paint.atom(i), alpha: paint.alpha(i) })
                .collect(),
            ..Scene::default()
        },
    };
    if let Some(c) = color {
        for p in &mut scene.spheres {
            p.color = c;
        }
        for k in &mut scene.capsules {
            (k.ca, k.cb) = (c, c);
        }
    }
    settle(&mut scene, style);
    scene
}

/// The edges of the simulation cell, as fine dark rods sized to the box.
/// A crystal with the periodic copies of its atoms that fall on the cell's far
/// faces, edges and corners, so the cell is drawn full the way crystallographers
/// draw it: rock salt's eight atoms, all in one corner of the cell, become the
/// familiar cube of twenty-seven. Also returned, the atom each copy repeats.
///
/// Only for crystal files read as one structure. A simulation box's atoms are
/// not on lattice sites, a protein's cell is not drawn, and the copies of a
/// trajectory would have to follow it frame by frame.
fn with_boundary(s: &Structure) -> Option<(Structure, Vec<u32>)> {
    /// Within this fraction of a cell edge counts as on the face.
    const ON: f32 = 1e-3;
    const LIMIT: usize = 100_000;
    let c = s.cell?;
    let crystal = matches!(s.format, Some("CIF" | "VASP" | "XSF"));
    if !crystal || s.has_polymer() || s.frame_count() > 1 || s.atoms.len() > LIMIT {
        return None;
    }
    let [a, b, cv] = c.vectors;
    let (bc, ca, ab) = (cross(b, cv), cross(cv, a), cross(a, b));
    let det = dot(a, bc);
    if det.abs() < 1e-9 {
        return None;
    }

    let mut out = s.clone();
    let mut images = Vec::new();
    for (i, atom) in s.atoms.iter().enumerate() {
        let d = sub(atom.pos, c.origin);
        let f = [dot(d, bc) / det, dot(d, ca) / det, dot(d, ab) / det];
        // The shifts along each edge that land on the cell: none, and one
        // across when the atom sits on a face.
        let shifts = f.map(|x| {
            if x.abs() < ON {
                [0.0, 1.0]
            } else if (x - 1.0).abs() < ON {
                [0.0, -1.0]
            } else {
                [0.0, 0.0]
            }
        });
        for n in 1..8 {
            let t = [
                shifts[0][n & 1],
                shifts[1][(n >> 1) & 1],
                shifts[2][(n >> 2) & 1],
            ];
            // A shift of zero on a bit that asked for one is a repeat of a
            // smaller combination, already made.
            if (0..3).any(|k| (n >> k) & 1 == 1 && t[k] == 0.0) {
                continue;
            }
            let mut copy = atom.clone();
            for k in 0..3 {
                copy.pos[k] += t[0] * a[k] + t[1] * b[k] + t[2] * cv[k];
            }
            out.atoms.push(copy);
            images.push(i as u32);
        }
    }
    if images.is_empty() {
        return None;
    }
    // The structure's own bonds stand; the copies take theirs from distance.
    let own = s.atoms.len() as u32;
    out.bonds
        .extend(crate::bonds::infer(&out.atoms).into_iter().filter(|&[x, y]| x >= own || y >= own));
    Some((out, images))
}

fn cell_edges(cell: &crate::protein::Cell) -> Vec<Capsule> {
    // Lifted again with every colour on a dark page, which takes it to a
    // pale grey there.
    const EDGE: [f32; 3] = [0.30, 0.31, 0.34];
    let longest = cell.lengths().into_iter().fold(0.0f32, f32::max);
    let r = (longest * 0.0015).clamp(0.01, 0.5);
    cell.edges()
        .into_iter()
        .map(|(a, b)| Capsule {
            a,
            b,
            r,
            ca: EDGE,
            cb: EDGE,
            alpha: 1.0,
        })
        .collect()
}

// -- small vector helpers ----------------------------------------------------

#[inline]
fn sub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
fn unit(v: [f32; 3]) -> [f32; 3] {
    let n = dot(v, v).sqrt();
    if n < 1e-6 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}
#[inline]
fn mix(a: [f32; 3], b: [f32; 3], t: f32) -> [f32; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

/// Alpha-carbon runs for one chain, split wherever the trace jumps.
fn ca_runs(s: &Structure, chain: &crate::protein::Chain, hidden: Option<&[bool]>) -> Vec<Vec<usize>> {
    let mut runs: Vec<Vec<usize>> = Vec::new();
    let mut run: Vec<usize> = Vec::new();
    let cut = |run: &mut Vec<usize>, runs: &mut Vec<Vec<usize>>| {
        if run.len() > 1 {
            runs.push(std::mem::take(run));
        } else {
            run.clear();
        }
    };
    for &ri in &chain.residues {
        let ri = ri as usize;
        // A residue hidden breaks the trace there, as a gap would.
        let Some(ca) = s.residues[ri].ca.filter(|&ca| !hidden.is_some_and(|h| h[ca as usize])) else {
            cut(&mut run, &mut runs);
            continue;
        };
        let p = s.atoms[ca as usize].pos;
        if let Some(&prev) = run.last() {
            let q = s.atoms[s.residues[prev].ca.unwrap() as usize].pos;
            // Consecutive alpha carbons sit 3.8 A apart; 5.5 A means a gap.
            if dot(sub(p, q), sub(p, q)) > 5.5 * 5.5 {
                cut(&mut run, &mut runs);
            }
        }
        run.push(ri);
    }
    cut(&mut run, &mut runs);
    runs
}

const TUBE_RADIUS: f32 = 0.75;

fn build_tube(
    s: &Structure,
    coloring: Coloring,
    quality: f32,
    palette: Palette,
    hidden: Option<&[bool]>,
) -> (Vec<Capsule>, Vec<u32>) {
    let mut out = (Vec::new(), Vec::new());
    for chain in &s.chains {
        for run in ca_runs(s, chain, hidden) {
            let run: Vec<(usize, [f32; 3])> = run
                .iter()
                .map(|&ri| (ri, s.atoms[s.residues[ri].ca.unwrap() as usize].pos))
                .collect();
            spline_tube(s, &run, TUBE_RADIUS, coloring, quality, palette, &mut out);
        }
    }
    out
}

/// A smooth tube through points of a trace, each its residue's, sampled
/// along a Catmull-Rom spline and joined with capsules, into `out`'s
/// capsules and the residue each belongs to.
#[allow(clippy::too_many_arguments)]
fn spline_tube(
    s: &Structure,
    run: &[(usize, [f32; 3])],
    radius: f32,
    coloring: Coloring,
    quality: f32,
    palette: Palette,
    out: &mut (Vec<Capsule>, Vec<u32>),
) {
    if run.len() < 2 {
        return;
    }
    let steps = (14.0 * quality).round().max(2.0) as usize;
    let at = |i: isize| -> [f32; 3] { run[i.clamp(0, run.len() as isize - 1) as usize].1 };
    let mut points: Vec<([f32; 3], [f32; 3])> = Vec::with_capacity(run.len() * steps + 1);
    // Each sample belongs to the nearer of the two residues it lies between.
    let mut belongs: Vec<u32> = Vec::with_capacity(run.len() * steps + 1);
    for i in 0..run.len() - 1 {
        let (p0, p1, p2, p3) = (at(i as isize - 1), at(i as isize), at(i as isize + 1), at(i as isize + 2));
        let c0 = residue_color(s, run[i].0, coloring, palette);
        let c1 = residue_color(s, run[i + 1].0, coloring, palette);
        for k in 0..steps {
            let t = k as f32 / steps as f32;
            belongs.push(run[if t < 0.5 { i } else { i + 1 }].0 as u32);
            points.push((catmull_rom(p0, p1, p2, p3, t), mix(c0, c1, t)));
        }
    }
    let last = run[run.len() - 1];
    points.push((last.1, residue_color(s, last.0, coloring, palette)));
    belongs.push(last.0 as u32);
    out.1.extend(&belongs[..belongs.len() - 1]);
    for pair in points.windows(2) {
        out.0.push(Capsule { a: pair[0].0, b: pair[1].0, r: radius, ca: pair[0].1, cb: pair[1].1, alpha: 1.0 });
    }
}

/// A nucleic acid strand's backbone is drawn wider than a protein's.
const NUCLEIC_RADIUS: f32 = 1.0;
const RUNG_RADIUS: f32 = 0.45;

/// The runs of a chain's nucleic acid, residue by residue, where each is
/// linked to the next and none is hidden.
fn nucleic_runs(s: &Structure, chain: &crate::protein::Chain, hidden: Option<&[bool]>) -> Vec<Vec<usize>> {
    let mut runs: Vec<Vec<usize>> = Vec::new();
    let mut run: Vec<usize> = Vec::new();
    for &ri in &chain.residues {
        let ri = ri as usize;
        let at = s.residues[ri].nucleic.filter(|&a| !hidden.is_some_and(|h| h[a as usize]));
        match at {
            Some(a) => {
                // Neighbouring phosphates sit 6 to 7 Å apart.
                let far = run.last().is_some_and(|&prev| {
                    let q = s.atoms[s.residues[prev].nucleic.unwrap() as usize].pos;
                    let p = s.atoms[a as usize].pos;
                    dot(sub(p, q), sub(p, q)) > 9.0 * 9.0
                });
                if far && run.len() > 1 {
                    runs.push(std::mem::take(&mut run));
                } else if far {
                    run.clear();
                }
                run.push(ri);
            }
            None => {
                if run.len() > 1 {
                    runs.push(std::mem::take(&mut run));
                }
                run.clear();
            }
        }
    }
    if run.len() > 1 {
        runs.push(run);
    }
    runs
}

/// Where a base's rung ends: N1 of a purine, N3 of a pyrimidine, the atoms
/// that pair across the helix. Found once, when the file is read.
fn base_tip(s: &Structure, residue: usize) -> Option<usize> {
    s.residues[residue].base.map(|a| a as usize)
}

/// Nucleic acids as a tube along the backbone, through each phosphorus,
/// with a rung out to each base, as the cartoons of DNA and RNA are drawn.
fn build_nucleic(
    s: &Structure,
    coloring: Coloring,
    quality: f32,
    palette: Palette,
    hidden: Option<&[bool]>,
) -> (Vec<Capsule>, Vec<u32>) {
    let mut out = (Vec::new(), Vec::new());
    for chain in &s.chains {
        for run in nucleic_runs(s, chain, hidden) {
            let trace: Vec<(usize, [f32; 3])> = run
                .iter()
                .map(|&ri| (ri, s.atoms[s.residues[ri].nucleic.unwrap() as usize].pos))
                .collect();
            spline_tube(s, &trace, NUCLEIC_RADIUS, coloring, quality, palette, &mut out);
            for &(ri, p) in &trace {
                if let Some(tip) = base_tip(s, ri) {
                    let c = residue_color(s, ri, coloring, palette);
                    out.0.push(Capsule { a: p, b: s.atoms[tip].pos, r: RUNG_RADIUS, ca: c, cb: c, alpha: 1.0 });
                    out.1.push(ri as u32);
                }
            }
        }
    }
    out
}


// ---------------------------------------------------------------------------
// Cartoon
// ---------------------------------------------------------------------------

/// Half-width and half-thickness of the swept cross-section, per structure.
const COIL_R: f32 = 0.30;
const HELIX: (f32, f32) = (1.15, 0.22);
const SHEET: (f32, f32) = (1.20, 0.22);
/// The arrowhead flares to this before tapering to a point.
const ARROW_HW: f32 = 1.95;
/// Sides of the swept cross-section. Enough that a coil reads as round.
const SIDES: usize = 10;

/// A ribbon that follows the raw alpha carbons of a strand shows the pleat, so
/// sheet residues get a light smoothing pass first. Helices keep their spiral.
fn smooth_sheets(s: &Structure, run: &[usize], p: &mut [[f32; 3]]) {
    for _ in 0..2 {
        let src = p.to_vec();
        for k in 1..run.len() - 1 {
            if s.residues[run[k]].ss != Ss::Sheet {
                continue;
            }
            p[k] = [
                (src[k - 1][0] + 2.0 * src[k][0] + src[k + 1][0]) * 0.25,
                (src[k - 1][1] + 2.0 * src[k][1] + src[k + 1][1]) * 0.25,
                (src[k - 1][2] + 2.0 * src[k][2] + src[k + 1][2]) * 0.25,
            ];
        }
    }
}

fn build_cartoon(
    s: &Structure,
    coloring: Coloring,
    quality: f32,
    palette: Palette,
    hidden: Option<&[bool]>,
) -> (Vec<Tri>, Vec<u32>) {
    let mut tris = Vec::new();
    let mut owners = Vec::new();
    let steps = (8.0 * quality).round().max(2.0) as usize;

    for chain in &s.chains {
        for run in ca_runs(s, chain, hidden) {
            if run.len() < 3 {
                continue;
            }
            let n = run.len();
            let mut p: Vec<[f32; 3]> = run
                .iter()
                .map(|&ri| s.atoms[s.residues[ri].ca.unwrap() as usize].pos)
                .collect();
            smooth_sheets(s, &run, &mut p);

            // Per-residue frame. The ribbon's flat face is set by the carbonyl
            // direction, which is what makes a helix read as a coiled ribbon
            // rather than a twisted one.
            let mut side: Vec<[f32; 3]> = Vec::with_capacity(n);
            for k in 0..n {
                let t = if k == 0 {
                    unit(sub(p[1], p[0]))
                } else if k == n - 1 {
                    unit(sub(p[n - 1], p[n - 2]))
                } else {
                    unit(sub(p[k + 1], p[k - 1]))
                };
                let r = &s.residues[run[k]];
                let raw = match (r.c, r.o) {
                    (Some(c), Some(o)) => sub(s.atoms[o as usize].pos, s.atoms[c as usize].pos),
                    // No carbonyl: fall back to the local bend.
                    _ if k > 0 && k + 1 < n => cross(sub(p[k + 1], p[k]), sub(p[k - 1], p[k])),
                    _ => [0.0, 0.0, 1.0],
                };
                let mut w = unit(sub(raw, [t[0] * dot(raw, t), t[1] * dot(raw, t), t[2] * dot(raw, t)]));
                if dot(w, w) < 0.5 {
                    // Degenerate: any vector perpendicular to the tangent.
                    let a = if t[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
                    w = unit(cross(t, a));
                }
                // A strand alternates its carbonyls; keep the face continuous.
                if let Some(&prev) = side.last() {
                    if dot(w, prev) < 0.0 {
                        w = [-w[0], -w[1], -w[2]];
                    }
                }
                side.push(w);
            }

            let ss = |k: usize| s.residues[run[k]].ss;
            let profile = |k: usize| match ss(k) {
                Ss::Helix => HELIX,
                Ss::Sheet => SHEET,
                Ss::Coil => (COIL_R, COIL_R),
            };
            let arrow = |k: usize| ss(k) == Ss::Sheet && (k + 1 >= n || ss(k + 1) != Ss::Sheet);

            // Sample the spline, carrying the frame and cross-section with it.
            struct Sample {
                pos: [f32; 3],
                side: [f32; 3],
                hw: f32,
                ht: f32,
                color: [f32; 3],
                /// The nearer residue of the two the sample lies between.
                residue: u32,
            }
            let at = |i: isize| p[i.clamp(0, n as isize - 1) as usize];
            let mut samples: Vec<Sample> = Vec::with_capacity(n * steps + 1);
            for k in 0..n - 1 {
                let c0 = residue_color(s, run[k], coloring, palette);
                let c1 = residue_color(s, run[k + 1], coloring, palette);
                let (hw0, ht0) = profile(k);
                let (hw1, ht1) = profile(k + 1);
                for step in 0..steps {
                    let f = step as f32 / steps as f32;
                    let (hw, ht) = if arrow(k) {
                        // Flare at the head, then taper to a point.
                        (ARROW_HW * (1.0 - f) + 0.12 * f, ht0)
                    } else if arrow(k + 1) {
                        // Hold the strand width; the flare is a deliberate step.
                        (hw0, ht0)
                    } else {
                        (hw0 + (hw1 - hw0) * f, ht0 + (ht1 - ht0) * f)
                    };
                    samples.push(Sample {
                        pos: catmull_rom(
                            at(k as isize - 1),
                            at(k as isize),
                            at(k as isize + 1),
                            at(k as isize + 2),
                            f,
                        ),
                        side: mix(side[k], side[k + 1], f),
                        hw,
                        ht,
                        color: mix(c0, c1, f),
                        residue: run[if f < 0.5 { k } else { k + 1 }] as u32,
                    });
                }
            }
            let (hw, ht) = profile(n - 1);
            samples.push(Sample {
                pos: p[n - 1],
                side: side[n - 1],
                hw: if arrow(n - 1) { 0.12 } else { hw },
                ht,
                color: residue_color(s, run[n - 1], coloring, palette),
                residue: run[n - 1] as u32,
            });
            if samples.len() < 2 {
                continue;
            }

            // Turn each sample into a ring of points on an elliptical section.
            let count = samples.len();
            let mut rings: Vec<[([f32; 3], [f32; 3]); SIDES]> = Vec::with_capacity(count);
            for (i, sm) in samples.iter().enumerate() {
                let tan = if i == 0 {
                    unit(sub(samples[1].pos, sm.pos))
                } else if i == count - 1 {
                    unit(sub(sm.pos, samples[count - 2].pos))
                } else {
                    unit(sub(samples[i + 1].pos, samples[i - 1].pos))
                };
                let w = unit(sub(
                    sm.side,
                    [
                        tan[0] * dot(sm.side, tan),
                        tan[1] * dot(sm.side, tan),
                        tan[2] * dot(sm.side, tan),
                    ],
                ));
                let up = cross(tan, w);
                rings.push(std::array::from_fn(|j| {
                    let a = std::f32::consts::TAU * j as f32 / SIDES as f32;
                    let (sa, ca) = a.sin_cos();
                    let pos = [
                        sm.pos[0] + w[0] * sm.hw * ca + up[0] * sm.ht * sa,
                        sm.pos[1] + w[1] * sm.hw * ca + up[1] * sm.ht * sa,
                        sm.pos[2] + w[2] * sm.hw * ca + up[2] * sm.ht * sa,
                    ];
                    // Outward normal of an ellipse swaps the two half-axes.
                    let nrm = unit([
                        w[0] * sm.ht * ca + up[0] * sm.hw * sa,
                        w[1] * sm.ht * ca + up[1] * sm.hw * sa,
                        w[2] * sm.ht * ca + up[2] * sm.hw * sa,
                    ]);
                    (pos, nrm)
                }));
            }

            for i in 0..count - 1 {
                let (a, b) = (&rings[i], &rings[i + 1]);
                let color = samples[i].color;
                owners.extend([samples[i].residue; 2 * SIDES]);
                for j in 0..SIDES {
                    let k = (j + 1) % SIDES;
                    tris.push(Tri {
                        p: [a[j].0, a[k].0, b[k].0],
                        n: [a[j].1, a[k].1, b[k].1],
                        c: color,
                    });
                    tris.push(Tri {
                        p: [a[j].0, b[k].0, b[j].0],
                        n: [a[j].1, b[k].1, b[j].1],
                        c: color,
                    });
                }
            }

            // Close both ends so the ribbon is not hollow.
            for (ring, color, residue) in [
                (&rings[0], samples[0].color, samples[0].residue),
                (&rings[count - 1], samples[count - 1].color, samples[count - 1].residue),
            ] {
                owners.extend(std::iter::repeat_n(residue, SIDES - 2));
                let nrm = unit(cross(
                    sub(ring[1].0, ring[0].0),
                    sub(ring[2].0, ring[0].0),
                ));
                for j in 1..SIDES - 1 {
                    tris.push(Tri {
                        p: [ring[0].0, ring[j].0, ring[j + 1].0],
                        n: [nrm; 3],
                        c: color,
                    });
                }
            }
        }
    }

    (tris, owners)
}

// ---------------------------------------------------------------------------
// Rasterization
// ---------------------------------------------------------------------------

pub struct Frame {
    pub width: usize,
    pub height: usize,
    /// BGRA, which is what GPUI uploads directly.
    pub bgra: Vec<u8>,
}

/// Light from the upper left, slightly in front of the viewer.
const LIGHT: [f32; 3] = [-0.4082, 0.5585, 0.7226];
/// Halfway between the light and the view direction, for the highlight.
const HALF: [f32; 3] = [-0.2225, 0.3044, 0.9262];

#[inline]
fn shade(n: [f32; 3], base: [f32; 3], style: Style) -> [f32; 3] {
    let diff = (n[0] * LIGHT[0] + n[1] * LIGHT[1] + n[2] * LIGHT[2]).max(0.0);
    let nh = (n[0] * HALF[0] + n[1] * HALF[1] + n[2] * HALF[2]).max(0.0);
    let n2 = nh * nh;
    let n8 = n2 * n2 * n2 * n2;

    let (k, spec) = match style {
        Style::Soft => (0.34 + 0.70 * diff, n8 * n8 * n8 * 0.28),
        Style::Illustrative => {
            // Quantise the light into a few steps, then ease most of the way
            // back: enough banding to read as drawn, not enough to look harsh.
            const STEPS: f32 = 4.0;
            let banded = (diff * STEPS).floor() / (STEPS - 1.0);
            let mixed = diff * 0.35 + banded.min(1.0) * 0.65;
            (0.62 + 0.42 * mixed, n8 * n8 * 0.10)
        }
        // Barely lit: the shape comes from the outline, with just enough
        // fall-off that a sphere does not read as a disc.
        Style::Flat => (0.76 + 0.26 * diff, 0.0),
        // A glossy surface: deep shade on the far side, a tight highlight
        // over a broad sheen.
        Style::Studio => {
            let n32 = n8 * n8 * n8 * n8;
            (0.26 + 0.80 * diff, n32 * 0.5 + n8 * 0.04)
        }
    };
    [base[0] * k + spec, base[1] * k + spec, base[2] * k + spec]
}

/// One horizontal band of the frame: its slices of the colour and depth
/// buffers, from row `top`. Bands are drawn in parallel, each with only the
/// primitives that reach into it.
/// The face a slab cuts, looking straight back at the viewer, and how much
/// darker it is than the outside, so a cut reads as one.
const CAP: [f32; 3] = [0.0, 0.0, 1.0];
const CAP_SHADE: f32 = 0.82;

struct Band<'a> {
    /// Index of the band's first pixel in the frame.
    offset: usize,
    color: &'a mut [f32],
    depth: &'a mut [f32],
    style: Style,
    /// The slab: nothing nearer than the first or further than the second.
    near: f32,
    far: f32,
}

impl Band<'_> {
    #[inline]
    fn plot(&mut self, idx: usize, z: f32, n: [f32; 3], base: [f32; 3]) {
        let idx = idx - self.offset;
        if z <= self.depth[idx] || z > self.near || z < self.far {
            return;
        }
        self.depth[idx] = z;
        let c = shade(n, base, self.style);
        let o = idx * 3;
        self.color[o] = c[0];
        self.color[o + 1] = c[1];
        self.color[o + 2] = c[2];
    }

    /// A ribbon triangle, in view space with its vertex normals.
    fn tri(&mut self, grid: &Grid, v: &[[f32; 3]; 3], nv: &[[f32; 3]; 3], color: [f32; 3]) {
        let sx: [f32; 3] = std::array::from_fn(|i| grid.cx + v[i][0] * grid.scale);
        let sy: [f32; 3] = std::array::from_fn(|i| grid.cy - v[i][1] * grid.scale);

        let d = (sy[1] - sy[2]) * (sx[0] - sx[2]) + (sx[2] - sx[1]) * (sy[0] - sy[2]);
        if d.abs() < 1e-9 {
            return;
        }
        let inv = 1.0 / d;

        let Some((x0, x1, y0, y1)) = grid.clip(
            sx[0].min(sx[1]).min(sx[2]),
            sx[0].max(sx[1]).max(sx[2]),
            sy[0].min(sy[1]).min(sy[2]),
            sy[0].max(sy[1]).max(sy[2]),
        ) else {
            return;
        };

        for py in y0..y1 {
            let fy = py as f32 + 0.5;
            let row = py * grid.width;
            for px in x0..x1 {
                let fx = px as f32 + 0.5;
                // Barycentric coordinates; dividing by the signed area keeps
                // them correct whichever way the triangle winds.
                let l0 = ((sy[1] - sy[2]) * (fx - sx[2]) + (sx[2] - sx[1]) * (fy - sy[2])) * inv;
                let l1 = ((sy[2] - sy[0]) * (fx - sx[2]) + (sx[0] - sx[2]) * (fy - sy[2])) * inv;
                let l2 = 1.0 - l0 - l1;
                if l0 < 0.0 || l1 < 0.0 || l2 < 0.0 {
                    continue;
                }

                // Orthographic projection is affine, so depth interpolates
                // linearly and exactly.
                let z = l0 * v[0][2] + l1 * v[1][2] + l2 * v[2][2];
                let idx = row + px;
                if z <= self.depth[idx - self.offset] {
                    continue;
                }

                let mut n = [
                    l0 * nv[0][0] + l1 * nv[1][0] + l2 * nv[2][0],
                    l0 * nv[0][1] + l1 * nv[1][1] + l2 * nv[2][1],
                    l0 * nv[0][2] + l1 * nv[1][2] + l2 * nv[2][2],
                ];
                let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                if len < 1e-6 {
                    continue;
                }
                n = [n[0] / len, n[1] / len, n[2] / len];
                // A ribbon is seen from both sides; light whichever faces us.
                if n[2] < 0.0 {
                    n = [-n[0], -n[1], -n[2]];
                }
                self.plot(idx, z, n, color);
            }
        }
    }
}

/// Ray-sphere in orthographic view space: returns the near-surface depth and
/// normal for the ray through `(x, y)`.
#[inline]
fn sphere_hit(c: [f32; 3], r: f32, x: f32, y: f32) -> Option<(f32, [f32; 3])> {
    let dx = x - c[0];
    let dy = y - c[1];
    let d2 = dx * dx + dy * dy;
    let rr = r * r;
    if d2 > rr {
        return None;
    }
    let dz = (rr - d2).sqrt();
    Some((c[2] + dz, [dx / r, dy / r, dz / r]))
}

/// Draw a frame: on the graphics card when there is one, else on the CPU.
pub fn render(
    scene: &Scene,
    cam: &Camera,
    width: usize,
    height: usize,
    page: [f32; 3],
    style: Style,
) -> Frame {
    crate::gpu::get()
        .and_then(|gpu| gpu.render(scene, cam, width, height, page, style))
        .unwrap_or_else(|| render_cpu(scene, cam, width, height, page, style))
}

/// A frame on no page at all, as straight RGBA: drawn once on white and
/// once on black, and the difference read as how much of the page shows
/// through each pixel. That holds for everything laid over the page,
/// edges, fog and glass alike, and needs nothing of either renderer.
pub fn render_transparent(scene: &Scene, cam: &Camera, width: usize, height: usize, style: Style) -> Vec<u8> {
    let white = render(scene, cam, width, height, [1.0; 3], style);
    let black = render(scene, cam, width, height, [0.0; 3], style);
    let mut out = vec![0u8; width * height * 4];
    for ((o, w), b) in out.chunks_exact_mut(4).zip(white.bgra.chunks_exact(4)).zip(black.bgra.chunks_exact(4)) {
        // BGRA in, RGBA out.
        let (w, b) = ([w[2], w[1], w[0]].map(|v| v as f32 / 255.0), [b[2], b[1], b[0]].map(|v| v as f32 / 255.0));
        let through = (0..3).map(|k| (w[k] - b[k]).clamp(0.0, 1.0)).sum::<f32>() / 3.0;
        let alpha = 1.0 - through;
        for k in 0..3 {
            let c = if alpha > 1e-3 { (b[k] / alpha).clamp(0.0, 1.0) } else { 0.0 };
            o[k] = (c * 255.0 + 0.5) as u8;
        }
        o[3] = (alpha * 255.0 + 0.5) as u8;
    }
    out
}

/// Whether frames are drawn on the graphics card.
pub fn accelerated() -> bool {
    crate::gpu::get().is_some()
}

pub fn render_cpu(
    scene: &Scene,
    cam: &Camera,
    width: usize,
    height: usize,
    page: [f32; 3],
    style: Style,
) -> Frame {
    use rayon::prelude::*;

    let n = width * height;
    let mut color = vec![0.0; n * 3];
    let mut depth = vec![f32::NEG_INFINITY; n];
    color.par_chunks_mut(3).for_each(|c| c.copy_from_slice(&page));

    if scene.is_empty() || width == 0 || height == 0 {
        return to_bgra(&color, width, height);
    }

    let m = cam.basis();
    let scale = (width.min(height) as f32) / (2.2 * cam.radius) * cam.zoom;
    let cx = width as f32 * 0.5 + cam.pan[0];
    let cy = height as f32 * 0.5 + cam.pan[1];
    let inv = 1.0 / scale;

    let to_view = |p: [f32; 3]| {
        let q = [p[0] - cam.center[0], p[1] - cam.center[1], p[2] - cam.center[2]];
        [
            m[0][0] * q[0] + m[0][1] * q[1] + m[0][2] * q[2],
            m[1][0] * q[0] + m[1][1] * q[1] + m[1][2] * q[2],
            m[2][0] * q[0] + m[2][1] * q[1] + m[2][2] * q[2],
        ]
    };
    let rot = |v: [f32; 3]| {
        [
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
        ]
    };

    let grid = Grid {
        cx,
        cy,
        scale,
        inv,
        width,
        height,
        top: 0,
        bottom: height,
    };
    let (near, far) = cam.depths();
    // Where the near plane cuts through a solid, it is drawn cut: a flat face
    // on the plane, so atoms sliced open read as solid rather than hollow.
    let cap = |front: f32, back: f32| (front > near && back <= near).then_some(near);

    // Everything in view space, with the rows it covers, opaque in the
    // order it was built and translucent held back to be laid over the
    // finished opaque frame back to front.
    let mut zmin = f32::MAX;
    let mut zmax = f32::MIN;
    let mut track = |lo: f32, hi: f32| {
        zmin = zmin.min(lo);
        zmax = zmax.max(hi);
    };
    let mut opaque: Vec<(Prim, (usize, usize))> = Vec::new();
    let mut translucent: Vec<(f32, Prim, (usize, usize))> = Vec::new();
    for (k, s) in scene.spheres.iter().enumerate() {
        if s.alpha <= 0.0 || s.r <= 0.0 {
            continue;
        }
        let (v, w) = cam.warp(to_view(s.c));
        let r = s.r * w;
        track(v[2] - s.r, v[2] + s.r);
        let Some(rows) = grid.sphere_rows(v, r) else {
            continue;
        };
        if s.alpha < OPAQUE {
            translucent.push((v[2], Prim::Sphere(k, v, r), rows));
        } else {
            opaque.push((Prim::Sphere(k, v, r), rows));
        }
    }
    for (k, c) in scene.capsules.iter().enumerate() {
        if c.alpha <= 0.0 {
            continue;
        }
        let (a, wa) = cam.warp(to_view(c.a));
        let (b, wb) = cam.warp(to_view(c.b));
        let r = c.r * (wa + wb) * 0.5;
        track(a[2].min(b[2]) - c.r, a[2].max(b[2]) + c.r);
        let Some(rows) = grid.capsule_rows(a, b, r) else {
            continue;
        };
        if c.alpha < OPAQUE {
            translucent.push(((a[2] + b[2]) * 0.5, Prim::Capsule(k, a, b, r), rows));
        } else {
            opaque.push((Prim::Capsule(k, a, b, r), rows));
        }
    }
    for (k, tri) in scene.tris.iter().enumerate() {
        let v = tri.p.map(|p| cam.warp(to_view(p)).0);
        track(
            v[0][2].min(v[1][2]).min(v[2][2]),
            v[0][2].max(v[1][2]).max(v[2][2]),
        );
        let top = cy - v[0][1].max(v[1][1]).max(v[2][1]) * scale;
        let bottom = cy - v[0][1].min(v[1][1]).min(v[2][1]) * scale;
        let y0 = top.floor().max(0.0) as usize;
        let y1 = (bottom.ceil().max(0.0) as usize).min(height);
        if y0 < y1 {
            let nv = [rot(tri.n[0]), rot(tri.n[1]), rot(tri.n[2])];
            opaque.push((Prim::Tri(k, v, nv), (y0, y1)));
        }
    }
    translucent.sort_by(|a, b| a.0.total_cmp(&b.0));

    // Bands a few to a thread, so a band crowded with atoms does not hold
    // the rest up.
    let bands = (rayon::current_num_threads() * 4).clamp(1, height);
    let rows = height.div_ceil(bands);
    let bin = |prims: &mut dyn Iterator<Item = (usize, (usize, usize))>| {
        let mut bins: Vec<Vec<u32>> = vec![Vec::new(); height.div_ceil(rows)];
        for (i, (y0, y1)) in prims {
            for b in bins.iter_mut().take((y1 - 1) / rows + 1).skip(y0 / rows) {
                b.push(i as u32);
            }
        }
        bins
    };
    let opaque_bins = bin(&mut opaque.iter().map(|p| p.1).enumerate());
    let translucent_bins = bin(&mut translucent.iter().map(|p| p.2).enumerate());

    color
        .par_chunks_mut(rows * width * 3)
        .zip(depth.par_chunks_mut(rows * width))
        .zip(&opaque_bins)
        .enumerate()
        .for_each(|(b, ((color, depth), bin))| {
            let top = b * rows;
            let grid = grid.rows(top, (top + rows).min(height));
            let mut t = Band {
                offset: top * width,
                color,
                depth,
                style,
                near,
                far,
            };
            for &i in bin {
                match &opaque[i as usize].0 {
                    Prim::Sphere(k, v, r) => {
                        let s = &scene.spheres[*k];
                        grid.sphere(*v, *r, |idx, z, n, _| match cap(z, 2.0 * v[2] - z) {
                            Some(z) => t.plot(idx, z, CAP, s.color.map(|c| c * CAP_SHADE)),
                            None => t.plot(idx, z, n, s.color),
                        });
                    }
                    Prim::Capsule(k, a, b, r) => {
                        let c = &scene.capsules[*k];
                        grid.capsule(*a, *b, *r, |idx, z, n, f| {
                            let base = mix(c.ca, c.cb, f);
                            let axis = a[2] + (b[2] - a[2]) * f;
                            match cap(z, 2.0 * axis - z) {
                                Some(z) => t.plot(idx, z, CAP, base.map(|c| c * CAP_SHADE)),
                                None => t.plot(idx, z, n, base),
                            }
                        });
                    }
                    Prim::Tri(k, v, nv) => t.tri(&grid, v, nv, scene.tris[*k].c),
                }
            }
        });

    let span = (zmax - zmin).max(1e-3);

    // Depth cue: wash distant geometry toward the paper. A drawn figure does
    // not fade with distance, so the outlines carry the depth instead.
    let fog_strength = style.fog();
    if fog_strength > 0.0 {
        color
            .par_chunks_mut(3)
            .zip(&depth)
            .filter(|(_, z)| **z != f32::NEG_INFINITY)
            .for_each(|(c, z)| {
                let d = ((z - zmin) / span).clamp(0.0, 1.0);
                let fog = (1.0 - d) * fog_strength;
                for k in 0..3 {
                    c[k] += (page[k] - c[k]) * fog;
                }
            });
    }

    // Shadow where atoms crowd together, for the studio style.
    if style == Style::Studio {
        let radius = ao_radius(width, height);
        let reach = radius / scale;
        color
            .par_chunks_mut(width * 3)
            .enumerate()
            .for_each(|(y, row)| {
                for x in 0..width {
                    let z = depth[y * width + x];
                    if z == f32::NEG_INFINITY {
                        continue;
                    }
                    let mut occ = 0.0;
                    for o in AO {
                        let qx = x as i32 + (o[0] * radius + 0.5).floor() as i32;
                        let qy = y as i32 + (o[1] * radius + 0.5).floor() as i32;
                        if qx < 0 || qy < 0 || qx >= width as i32 || qy >= height as i32 {
                            continue;
                        }
                        let nz = depth[qy as usize * width + qx as usize];
                        if nz != f32::NEG_INFINITY {
                            occ += occluding(nz - z, reach);
                        }
                    }
                    let k = 1.0 - AO_STRENGTH * occ / AO.len() as f32;
                    // The card keeps colour in eight bits, so a highlight is
                    // white before it is shaded; shaded here the same way.
                    for c in &mut row[x * 3..x * 3 + 3] {
                        *c = c.min(1.0) * k;
                    }
                }
            });
    }

    // Outline. Soft darkens what sits behind a depth jump and the rim of the
    // silhouette; illustrative lays down real ink, on both sides of the
    // silhouette, thickened so it holds at any resolution.
    let jump = (span * 0.012).max(0.45);
    match style {
        Style::Soft => {
            let src = color.clone();
            color
                .par_chunks_mut(width * 3)
                .enumerate()
                .for_each(|(y, row)| {
                    for x in 0..width {
                        let idx = y * width + x;
                        let z = depth[idx];
                        if z == f32::NEG_INFINITY {
                            continue;
                        }
                        let mut behind = false;
                        let mut open = false;
                        for nb in neighbours(idx, x, y, width, height).into_iter().flatten() {
                            let nz = depth[nb];
                            if nz == f32::NEG_INFINITY {
                                open = true;
                            } else if nz - z > jump {
                                behind = true;
                            }
                        }
                        let k = if behind {
                            0.62
                        } else if open {
                            0.80
                        } else {
                            continue;
                        };
                        for i in 0..3 {
                            row[x * 3 + i] = src[idx * 3 + i] * k;
                        }
                    }
                });
        }
        Style::Illustrative | Style::Flat => {
            const INK: [f32; 3] = [0.07, 0.07, 0.09];
            let mut edge = vec![false; n];
            edge.par_chunks_mut(width).enumerate().for_each(|(y, row)| {
                for (x, e) in row.iter_mut().enumerate() {
                    let idx = y * width + x;
                    let z = depth[idx];
                    *e = neighbours(idx, x, y, width, height)
                        .into_iter()
                        .flatten()
                        .any(|nb| {
                            let nz = depth[nb];
                            // The silhouette is inked from both sides, so the
                            // line reads the same against the page and
                            // against the body.
                            let silhouette =
                                (z == f32::NEG_INFINITY) != (nz == f32::NEG_INFINITY);
                            let overlap = z != f32::NEG_INFINITY
                                && nz != f32::NEG_INFINITY
                                && nz - z > jump;
                            silhouette || overlap
                        });
                }
            });

            // One pixel of line per ~700 of frame, so the weight is the same
            // whether or not the display is retina.
            let weight = (width.min(height) / 700).max(1);
            for _ in 1..weight {
                let src = edge.clone();
                edge.par_chunks_mut(width).enumerate().for_each(|(y, row)| {
                    for (x, e) in row.iter_mut().enumerate() {
                        let idx = y * width + x;
                        if !src[idx] {
                            *e = neighbours(idx, x, y, width, height)
                                .into_iter()
                                .flatten()
                                .any(|nb| src[nb]);
                        }
                    }
                });
            }

            color
                .par_chunks_mut(3)
                .zip(&edge)
                .filter(|(_, e)| **e)
                .for_each(|(c, _)| c.copy_from_slice(&INK));
        }
        Style::Studio => {}
    }

    // -- translucent atoms and bonds ----------------------------------------
    // Painter's order, farthest first, each blended where it is in front of
    // the opaque scene. They take the same fog, and a sphere keeps an ink
    // rim in the illustrative style so it still reads as drawn.
    if !translucent.is_empty() {
        let rim = (width.min(height) / 700).max(1) as f32 + 0.5;
        color
            .par_chunks_mut(rows * width * 3)
            .zip(depth.par_chunks(rows * width))
            .zip(&translucent_bins)
            .enumerate()
            .for_each(|(b, ((color, depth), bin))| {
                let top = b * rows;
                let offset = top * width;
                let grid = grid.rows(top, (top + rows).min(height));
                let mut blend = |idx: usize,
                                 z: f32,
                                 n: [f32; 3],
                                 base: [f32; 3],
                                 alpha: f32,
                                 edge: bool| {
                    let idx = idx - offset;
                    if z <= depth[idx] || z > near || z < far {
                        return;
                    }
                    let mut c = shade(n, base, style);
                    let d = ((z - zmin) / span).clamp(0.0, 1.0);
                    let fog = (1.0 - d) * fog_strength;
                    for k in 0..3 {
                        c[k] += (page[k] - c[k]) * fog;
                    }
                    let (c, alpha) = match (edge, style) {
                        (true, Style::Illustrative | Style::Flat) => ([0.07, 0.07, 0.09], alpha.max(0.45)),
                        (true, Style::Soft) => (c.map(|x| x * 0.8), alpha),
                        _ => (c, alpha),
                    };
                    // The card clamps what it blends in to eight bits.
                    let o = idx * 3;
                    for k in 0..3 {
                        color[o + k] = c[k].min(1.0) * alpha + color[o + k] * (1.0 - alpha);
                    }
                };
                for &i in bin {
                    match &translucent[i as usize].1 {
                        Prim::Sphere(k, v, r) => {
                            let s = &scene.spheres[*k];
                            let inner = (r * scale - rim).max(0.0);
                            grid.sphere(*v, *r, |idx, z, n, from_centre| {
                                blend(idx, z, n, s.color, s.alpha, from_centre > inner)
                            });
                        }
                        Prim::Capsule(k, a, b, r) => {
                            let c = &scene.capsules[*k];
                            grid.capsule(*a, *b, *r, |idx, z, n, f| {
                                blend(idx, z, n, mix(c.ca, c.cb, f), c.alpha, false)
                            });
                        }
                        Prim::Tri(..) => {}
                    }
                }
            });
    }

    to_bgra(&color, width, height)
}

/// Where the studio style looks for neighbours that shade a pixel: a spiral
/// of sixteen points over the unit disc. The GPU shader holds the same list.
const AO: [[f32; 2]; 16] = [
    [0.1768, 0.0000],
    [-0.2258, 0.2068],
    [0.0346, -0.3938],
    [0.2846, 0.3712],
    [-0.5222, -0.0924],
    [0.4947, -0.3147],
    [-0.1655, 0.6155],
    [-0.3156, -0.6076],
    [0.6846, 0.2500],
    [-0.7123, 0.2940],
    [0.3434, -0.7337],
    [0.2537, 0.8089],
    [-0.7647, -0.4432],
    [0.8971, -0.1972],
    [-0.5475, 0.7788],
    [-0.1265, -0.9761],
];

/// How dark a fully enclosed pixel goes.
const AO_STRENGTH: f32 = 0.7;

/// How far, in pixels, the studio style looks for neighbours: a fixed share
/// of the frame, so the shadows read the same at any size and in any units.
pub(crate) fn ao_radius(width: usize, height: usize) -> f32 {
    (width.min(height) as f32 * 0.02).max(2.0)
}

/// How much a neighbour `dz` nearer the viewer shades a pixel: more the
/// higher it stands, up to `reach`, then less again, so a body floating well
/// in front casts no halo on what is behind it. Smooth, so the CPU and GPU,
/// whose depths differ in the last bits, agree.
fn occluding(dz: f32, reach: f32) -> f32 {
    let up = (dz / reach).clamp(0.0, 1.0);
    let down = 1.0 - ((dz - reach) / reach).clamp(0.0, 1.0);
    up * down
}

/// Below this an atom is drawn as see-through.
pub(crate) const OPAQUE: f32 = 0.999;

/// A primitive in view space, with its index in the scene.
enum Prim {
    /// Which, where in view space as drawn, and how big there.
    Sphere(usize, [f32; 3], f32),
    Capsule(usize, [f32; 3], [f32; 3], f32),
    /// The corners and their normals.
    Tri(usize, [[f32; 3]; 3], [[f32; 3]; 3]),
}

/// The mapping between view space and the frame's pixels, and the two
/// impostor shapes resolved against it.
#[derive(Clone, Copy)]
struct Grid {
    cx: f32,
    cy: f32,
    scale: f32,
    inv: f32,
    width: usize,
    height: usize,
    /// The rows drawn to, for drawing one band.
    top: usize,
    bottom: usize,
}

impl Grid {
    fn rows(&self, top: usize, bottom: usize) -> Grid {
        Grid {
            top,
            bottom,
            ..*self
        }
    }

    /// A rectangle clamped to the frame and to the band's rows.
    #[inline]
    fn clip(
        &self,
        left: f32,
        right: f32,
        top: f32,
        bottom: f32,
    ) -> Option<(usize, usize, usize, usize)> {
        let (x0, x1, y0, y1) = bbox(left, right, top, bottom, self.width, self.height)?;
        let (y0, y1) = (y0.max(self.top), y1.min(self.bottom));
        (y0 < y1).then_some((x0, x1, y0, y1))
    }

    fn sphere_bounds(&self, v: [f32; 3], r: f32) -> [f32; 4] {
        let sx = self.cx + v[0] * self.scale;
        let sy = self.cy - v[1] * self.scale;
        let rp = r * self.scale;
        [sx - rp, sx + rp, sy - rp, sy + rp]
    }

    fn capsule_bounds(&self, a: [f32; 3], b: [f32; 3], r: f32) -> [f32; 4] {
        let sax = self.cx + a[0] * self.scale;
        let say = self.cy - a[1] * self.scale;
        let sbx = self.cx + b[0] * self.scale;
        let sby = self.cy - b[1] * self.scale;
        let rp = r * self.scale;
        [
            sax.min(sbx) - rp,
            sax.max(sbx) + rp,
            say.min(sby) - rp,
            say.max(sby) + rp,
        ]
    }

    /// The rows a sphere covers, or `None` if it is off the frame.
    fn sphere_rows(&self, v: [f32; 3], r: f32) -> Option<(usize, usize)> {
        let [l, r, t, b] = self.sphere_bounds(v, r);
        self.clip(l, r, t, b).map(|(_, _, y0, y1)| (y0, y1))
    }

    fn capsule_rows(&self, a: [f32; 3], b: [f32; 3], r: f32) -> Option<(usize, usize)> {
        let [l, r, t, b] = self.capsule_bounds(a, b, r);
        self.clip(l, r, t, b).map(|(_, _, y0, y1)| (y0, y1))
    }

    #[inline]
    fn view_x(&self, px: usize) -> f32 {
        (px as f32 + 0.5 - self.cx) * self.inv
    }

    #[inline]
    fn view_y(&self, py: usize) -> f32 {
        -(py as f32 + 0.5 - self.cy) * self.inv
    }

    /// Every pixel the sphere covers, with its depth, normal, and distance
    /// from the centre in pixels.
    fn sphere(&self, v: [f32; 3], r: f32, mut f: impl FnMut(usize, f32, [f32; 3], f32)) {
        let [l, rt, t, b] = self.sphere_bounds(v, r);
        let Some((x0, x1, y0, y1)) = self.clip(l, rt, t, b) else {
            return;
        };
        for py in y0..y1 {
            let y = self.view_y(py);
            let row = py * self.width;
            for px in x0..x1 {
                let x = self.view_x(px);
                if let Some((z, n)) = sphere_hit(v, r, x, y) {
                    let d = ((x - v[0]).powi(2) + (y - v[1]).powi(2)).sqrt() * self.scale;
                    f(row + px, z, n, d);
                }
            }
        }
    }

    /// Every pixel the capsule covers, with its depth, normal, and how far
    /// along from `a` to `b` the surface there is, for blending the colour.
    fn capsule(
        &self,
        a: [f32; 3],
        b: [f32; 3],
        r: f32,
        mut f: impl FnMut(usize, f32, [f32; 3], f32),
    ) {
        let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let u = if len > 1e-6 {
            [d[0] / len, d[1] / len, d[2] / len]
        } else {
            [0.0, 0.0, 1.0]
        };

        let [l, rt, t, bt] = self.capsule_bounds(a, b, r);
        let Some((x0, x1, y0, y1)) = self.clip(l, rt, t, bt) else {
            return;
        };

        // Coefficients of the ray-cylinder quadratic that do not vary per
        // pixel. `q` collapses to zero when the axis points at the viewer, in
        // which case the end caps alone describe the shape.
        let q = 1.0 - u[2] * u[2];
        let rr = r * r;

        for py in y0..y1 {
            let y = self.view_y(py);
            let row = py * self.width;
            for px in x0..x1 {
                let x = self.view_x(px);
                let wx = x - a[0];
                let wy = y - a[1];

                if q > 1e-6 && len > 1e-6 {
                    let k = wx * u[0] + wy * u[1];
                    let bq = -2.0 * k * u[2];
                    let cq = wx * wx + wy * wy - k * k - rr;
                    let disc = bq * bq - 4.0 * q * cq;
                    if disc >= 0.0 {
                        let cz = (-bq + disc.sqrt()) / (2.0 * q);
                        let along = k + cz * u[2];
                        let s = along / len;
                        if (0.0..=1.0).contains(&s) {
                            let nrm = [
                                (wx - u[0] * along) / r,
                                (wy - u[1] * along) / r,
                                (cz - u[2] * along) / r,
                            ];
                            f(row + px, a[2] + cz, nrm, s);
                            continue;
                        }
                    }
                }

                // Outside the cylinder's span, the caps take over.
                let mut best: Option<(f32, [f32; 3], f32)> = None;
                for (centre, s) in [(a, 0.0), (b, 1.0)] {
                    if let Some((z, nrm)) = sphere_hit(centre, r, x, y) {
                        if best.as_ref().is_none_or(|(bz, _, _)| z > *bz) {
                            best = Some((z, nrm, s));
                        }
                    }
                }
                if let Some((z, nrm, s)) = best {
                    f(row + px, z, nrm, s);
                }
            }
        }
    }
}


/// The four pixels sharing an edge with this one.
#[inline]
fn neighbours(
    idx: usize,
    x: usize,
    y: usize,
    width: usize,
    height: usize,
) -> [Option<usize>; 4] {
    [
        (x > 0).then(|| idx - 1),
        (x + 1 < width).then_some(idx + 1),
        (y > 0).then(|| idx - width),
        (y + 1 < height).then_some(idx + width),
    ]
}

/// Clamp a floating-point rectangle to the frame, or `None` if it falls
/// entirely outside.
#[inline]
fn bbox(
    left: f32,
    right: f32,
    top: f32,
    bottom: f32,
    width: usize,
    height: usize,
) -> Option<(usize, usize, usize, usize)> {
    let x0 = left.floor().max(0.0) as usize;
    let x1 = (right.ceil().max(0.0) as usize).min(width);
    let y0 = top.floor().max(0.0) as usize;
    let y1 = (bottom.ceil().max(0.0) as usize).min(height);
    if x0 >= x1 || y0 >= y1 {
        None
    } else {
        Some((x0, x1, y0, y1))
    }
}

fn to_bgra(color: &[f32], width: usize, height: usize) -> Frame {
    use rayon::prelude::*;
    let mut bgra = vec![0u8; width * height * 4];
    bgra.par_chunks_mut(4).zip(color.par_chunks(3)).for_each(|(q, c)| {
        q[0] = (c[2].clamp(0.0, 1.0) * 255.0 + 0.5) as u8;
        q[1] = (c[1].clamp(0.0, 1.0) * 255.0 + 0.5) as u8;
        q[2] = (c[0].clamp(0.0, 1.0) * 255.0 + 0.5) as u8;
        q[3] = 255;
    });
    Frame {
        width,
        height,
        bgra,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pixel(f: &Frame, x: usize, y: usize) -> [u8; 3] {
        let o = (y * f.width + x) * 4;
        [f.bgra[o + 2], f.bgra[o + 1], f.bgra[o]]
    }

    fn sample(name: &str) -> Structure {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("samples").join(name);
        crate::formats::load_bytes(name, std::fs::read(path).unwrap()).unwrap()
    }

    fn scene(s: &Structure, repr: Repr, coloring: Coloring, cell: bool) -> Scene {
        build_scene(
            s, repr, coloring, 1.0, Style::Soft, Palette::Spectrum, false, cell,
            &Appearance::default(),
        )
    }

    #[test]
    fn a_positive_roll_turns_the_view_clockwise() {
        let mut cam = Camera { orientation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], zoom: 1.0, pan: [0.0; 2], center: [0.0; 3], radius: 10.0, perspective: false, slab: None };
        cam.roll(std::f32::consts::FRAC_PI_2);
        // What was to the right now points down the screen (view y is up).
        let m = cam.basis();
        let v = [m[0][0], m[1][0], m[2][0]];
        assert!(v[0].abs() < 1e-5 && (v[1] + 1.0).abs() < 1e-5 && v[2].abs() < 1e-5, "{v:?}");
    }

    #[test]
    fn a_crystal_fills_its_cell() {
        let s = sample("rocksalt.cif");
        assert_eq!(s.atoms.len(), 8);
        let bare = scene(&s, Repr::BallAndStick, Coloring::Element, false);
        let full = scene(&s, Repr::BallAndStick, Coloring::Element, true);
        assert_eq!(bare.spheres.len(), 8);
        // Fourteen sodium on the corners and faces, thirteen chloride on the
        // edges and at the centre.
        assert_eq!(full.spheres.len(), 27);
        // Every Na-Cl contact in the cube: nine lines of three atoms along
        // each axis, two bonds to a line, each bond drawn as two halves.
        assert_eq!(full.capsules.len() - 12, 54 * 2);
        // Every atom lies within the cell.
        let c = s.cell.unwrap();
        let edge = c.lengths()[0];
        for p in &full.spheres {
            for k in 0..3 {
                let x = p.c[k] - c.origin[k];
                assert!((-1e-3..=edge + 1e-3).contains(&x), "an atom left the cell: {:?}", p.c);
            }
        }
        // The copies wear their originals' colours, and the colours do not
        // change when the cell is shown.
        let colours = |sc: &Scene| {
            let mut v: Vec<[u32; 3]> = sc.spheres.iter().map(|p| p.color.map(f32::to_bits)).collect();
            v.sort();
            v.dedup();
            v
        };
        assert_eq!(colours(&bare), colours(&full));
        // The structure itself is untouched.
        assert_eq!(s.atoms.len(), 8);
    }

    #[test]
    fn only_crystals_take_boundary_copies() {
        // Proteins, molecules and simulation boxes draw the same atoms with
        // the cell shown as without; at most the cell's twelve edges are added.
        for name in ["1crn.pdb.gz", "caffeine.sdf", "c60.xyz", "lj-melt.lammpstrj.gz"] {
            let s = sample(name);
            for repr in [Repr::Spacefill, Repr::BallAndStick, Repr::Licorice, Repr::Cartoon] {
                for coloring in [Coloring::Element, Coloring::Rainbow, Coloring::Chain] {
                    let bare = scene(&s, repr, coloring, false);
                    let full = scene(&s, repr, coloring, true);
                    assert_eq!(bare.spheres.len(), full.spheres.len(), "{name} {repr:?}");
                    assert!(bare.spheres.iter().zip(&full.spheres).all(|(a, b)| a.c == b.c && a.color == b.color));
                    let edges = if s.cell.is_some() { 12 } else { 0 };
                    assert_eq!(full.capsules.len(), bare.capsules.len() + edges, "{name} {repr:?}");
                }
            }
        }
    }

    #[test]
    fn palettes_are_listed_by_group() {
        let mut seen = Vec::new();
        for p in Palette::ALL {
            if seen.last() != Some(&p.group()) {
                assert!(!seen.contains(&p.group()), "{p:?} is apart from its group");
                seen.push(p.group());
            }
        }
        assert_eq!(seen.len(), 4);
        for (i, p) in Palette::ALL.iter().enumerate() {
            assert!(
                Palette::ALL[..i].iter().all(|q| q.label() != p.label()),
                "two palettes are called {}",
                p.label()
            );
        }
    }

    #[test]
    fn chains_are_told_apart_in_every_palette() {
        // Chains beside each other in a structure are the ones that touch, so
        // each must stand clear of the next, and the first four, which most
        // complexes stop at, of each other.
        let apart = |a: [f32; 3], b: [f32; 3]| {
            ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
        };
        for p in Palette::ALL {
            for i in 0..8 {
                for j in i + 1..8 {
                    if j == i + 1 || j < 4 {
                        let d = apart(p.chain(i), p.chain(j));
                        assert!(d > 0.2, "{p:?}: chains {i} and {j} are {d:.3} apart");
                    }
                }
            }
        }
    }

    #[test]
    fn palettes_recolour_elements_unless_fixed() {
        let s = crate::formats::load_bytes(
            "caffeine.sdf",
            include_bytes!("../samples/caffeine.sdf").to_vec(),
        )
        .unwrap();
        let colors = |palette, fixed| {
            build_scene(
                &s,
                Repr::BallAndStick,
                Coloring::Element,
                1.0,
                Style::Soft,
                palette,
                fixed,
                false,
                &Appearance::default(),
            )
            .spheres
            .iter()
            .map(|p| p.color)
            .collect::<Vec<_>>()
        };
        // Fixed, elements keep their standard colours whatever the palette,
        // deepened as the style deepens every colour.
        let standard = colors(Palette::Spectrum, true);
        assert_eq!(standard[0], rich(crate::elements::Element::O.color()));
        assert_eq!(colors(Palette::Viridis, true), standard);
        // Otherwise each palette paints them its own way.
        assert_ne!(colors(Palette::Viridis, false), standard);
        assert_ne!(colors(Palette::Magma, false), colors(Palette::Viridis, false));
    }

    #[test]
    fn ligands_keep_element_colours_beside_chains() {
        let s = crate::formats::load_bytes(
            "4hhb.pdb.gz",
            include_bytes!("../samples/4hhb.pdb.gz").to_vec(),
        )
        .unwrap();
        assert!(s.has_ligands());
        let pink = [1.0, 0.0, 1.0];
        let mut looks = Appearance::default();
        looks.set(Species::Element(6), Look { color: Some(pink), ..Look::default() });
        let colors = build_scene(
            &s,
            Repr::Spacefill,
            Coloring::Chain,
            1.0,
            Style::Soft,
            Palette::Spectrum,
            false,
            false,
            &looks,
        )
        .spheres
        .iter()
        .map(|p| p.color)
        .collect::<Vec<_>>();
        // The haems' carbons take the chosen colour; the chains' keep theirs.
        for (a, c) in s.atoms.iter().zip(&colors) {
            let r = &s.residues[a.residue as usize];
            if r.name == "HEM" && a.element.0 == 6 {
                assert_eq!(*c, rich(pink));
            } else if r.ca.is_some() {
                assert_eq!(*c, rich(Palette::Spectrum.chain(r.chain as usize)));
            }
        }
    }

    #[test]
    fn a_ligand_takes_its_own_colour() {
        let mut s = crate::formats::load_bytes(
            "4hhb.pdb.gz",
            include_bytes!("../samples/4hhb.pdb.gz").to_vec(),
        )
        .unwrap();
        assert_eq!(s.ligands().iter().find(|(n, _)| n == "HEM").map(|l| l.1), Some(4));
        assert!(s.ligands().iter().all(|(n, _)| n != "HOH"), "water is no ligand");
        let gold = [0.98, 0.78, 0.20];
        let paint = |s: &Structure| {
            build_scene(
                s,
                Repr::Spacefill,
                Coloring::Chain,
                1.0,
                Style::Soft,
                Palette::Spectrum,
                true,
                false,
                &Appearance::default(),
            )
            .spheres
            .iter()
            .map(|p| p.color)
            .collect::<Vec<_>>()
        };
        let haem = |s: &Structure, i: usize, z: u8| {
            let a = &s.atoms[i];
            s.residues[a.residue as usize].name == "HEM" && a.element.0 == z
        };

        s.ligand_paint.insert("HEM".into(), LigandPaint::Carbons(gold));
        let colors = paint(&s);
        let nitrogen = rich(crate::elements::Element(7).color());
        for (i, c) in colors.iter().enumerate() {
            if haem(&s, i, 6) {
                assert_eq!(*c, rich(gold));
            } else if haem(&s, i, 7) {
                assert_eq!(*c, nitrogen, "a haem's nitrogens keep their own colour");
            }
        }

        s.ligand_paint.insert("HEM".into(), LigandPaint::Solid(gold));
        let colors = paint(&s);
        for (i, c) in colors.iter().enumerate() {
            let r = &s.residues[s.atoms[i].residue as usize];
            if r.name == "HEM" {
                assert_eq!(*c, rich(gold));
            } else if r.ca.is_some() {
                assert_ne!(*c, rich(gold), "the paint spilt onto the protein");
            }
        }
    }

    /// Where an atom's centre lands in a `size` square frame, and its depth.
    fn on_screen(cam: &Camera, size: usize, p: [f32; 3]) -> (f32, f32, f32) {
        let m = cam.basis();
        let q = sub(p, cam.center);
        let v = [dot(m[0], q), dot(m[1], q), dot(m[2], q)];
        let scale = size as f32 / (2.2 * cam.radius) * cam.zoom;
        let half = size as f32 * 0.5;
        (half + cam.pan[0] + v[0] * scale, half + cam.pan[1] - v[1] * scale, v[2])
    }

    #[test]
    fn picking_finds_the_nearest_atom_under_the_pointer() {
        let s = sample("1crn.pdb.gz");
        let cam = Camera::frame(&s);
        let look = Appearance::default();
        // The atom nearest the viewer of all: nothing can be in front of it.
        let front = (0..s.atoms.len())
            .max_by(|&a, &b| on_screen(&cam, 400, s.atoms[a].pos).2.total_cmp(&on_screen(&cam, 400, s.atoms[b].pos).2))
            .unwrap();
        let (x, y, _) = on_screen(&cam, 400, s.atoms[front].pos);
        let hit = pick(&s, Repr::Spacefill, &look, &cam, 400, 400, x, y).unwrap();
        assert_eq!(hit.atom, front);
        assert!(!hit.residue);
        // Panned, the same atom is wherever it was panned to.
        let mut moved = cam;
        moved.pan = [37.0, -21.0];
        let hit = pick(&s, Repr::Spacefill, &look, &moved, 400, 400, x + 37.0, y - 21.0).unwrap();
        assert_eq!(hit.atom, front);
        // A corner of the frame is empty.
        assert_eq!(pick(&s, Repr::Spacefill, &look, &cam, 400, 400, 2.0, 2.0), None);
    }

    #[test]
    fn a_ribbon_picks_its_residue() {
        let s = sample("1crn.pdb.gz");
        let cam = Camera::frame(&s);
        let look = Appearance::default();
        let r = &s.residues[20];
        let ca = r.ca.unwrap() as usize;
        let (x, y, _) = on_screen(&cam, 400, s.atoms[ca].pos);
        let hit = pick(&s, Repr::Cartoon, &look, &cam, 400, 400, x, y).unwrap();
        assert!(hit.residue);
        // A residue in front may hide this one, but whatever is hit is an
        // alpha carbon standing for its residue.
        assert_eq!(s.residues[s.atoms[hit.atom].residue as usize].ca, Some(hit.atom as u32));
        // Hidden atoms are not there to be picked.
        let mut gone = Appearance::default();
        for z in 0..=255u8 {
            gone.set(Species::Element(z), Look { size: 0.0, ..Look::default() });
        }
        assert_eq!(pick(&s, Repr::Spacefill, &gone, &cam, 400, 400, x, y), None);
    }

    #[test]
    fn a_picked_residue_glows_without_new_atoms() {
        let s = sample("1crn.pdb.gz");
        let ca = s.residues[20].ca.unwrap() as usize;
        let picked = Picked { atom: ca, residue: true, depth: 0.0 };
        let look = Appearance::default();
        // On a ribbon nothing is added: the ribbon itself is painted.
        let ribbon = build_highlight(&s, picked, Repr::Cartoon, Style::Soft, &look);
        assert!(ribbon.is_empty());
        // In atoms, a glow on each of its atoms.
        let atoms = s.atoms.iter().filter(|a| a.residue == 20).count();
        let sticks = build_highlight(&s, picked, Repr::BallAndStick, Style::Soft, &look);
        assert_eq!(sticks.spheres.len(), atoms);
        assert!(sticks.capsules.is_empty());
        let one = build_highlight(&s, Picked { residue: false, ..picked }, Repr::Spacefill, Style::Soft, &look);
        assert_eq!(one.spheres.len(), 1);
        assert!(one.spheres[0].alpha < 1.0);
    }

    #[test]
    fn a_residue_is_painted_on_its_own_stretch_of_ribbon() {
        let s = sample("1crn.pdb.gz");
        for repr in [Repr::Cartoon, Repr::Tube] {
            let base = scene(&s, repr, Coloring::Chain, false);
            assert_eq!(base.tri_residue.len(), base.tris.len());
            assert_eq!(base.capsule_residue.len(), base.capsules.len());
            let mut painted = base.clone();
            painted.paint_residue(0..base.tris.len(), 0..base.capsules.len(), 20, Style::Soft, false);
            let white = highlight_color(Style::Soft);
            let changed: Vec<u32> = (0..base.tris.len())
                .filter(|&i| painted.tris[i].c != base.tris[i].c)
                .map(|i| base.tri_residue[i])
                .chain(
                    (0..base.capsules.len())
                        .filter(|&i| painted.capsules[i].ca != base.capsules[i].ca)
                        .map(|i| base.capsule_residue[i]),
                )
                .collect();
            assert!(!changed.is_empty(), "{repr:?}: nothing painted");
            assert!(changed.iter().all(|&r| r == 20), "{repr:?}: another residue painted");
            assert!(painted.tris.iter().filter(|t| t.c == white).count() + painted.capsules.iter().filter(|c| c.ca == white).count() > 0);
            // Most of the way is neither the ribbon's colour nor the highlight's.
            let mut faint = base.clone();
            faint.paint_residue(0..base.tris.len(), 0..base.capsules.len(), 20, Style::Soft, true);
            assert!(faint.tris.iter().all(|t| t.c != white));
            assert!(faint.capsules.iter().all(|c| c.ca != white));
        }
        // Merged with atoms beside it, the lists still line up.
        let mut both = scene(&s, Repr::BallAndStick, Coloring::Chain, false);
        both.merge(scene(&s, Repr::Tube, Coloring::Chain, false));
        assert_eq!(both.capsule_residue.len(), both.capsules.len());
        assert_eq!(both.capsule_residue.iter().filter(|&&r| r != NO_RESIDUE).count(), scene(&s, Repr::Tube, Coloring::Chain, false).capsules.len());
    }

    #[test]
    fn a_turn_brings_a_direction_to_face_the_viewer() {
        let cam = Camera { orientation: Camera::home(), zoom: 1.0, pan: [0.0; 2], center: [0.0; 3], radius: 10.0, perspective: false, slab: None };
        for d in [[1.0, 0.0, 0.0], [0.3, -0.8, 0.2], [0.0, 0.0, -1.0], [-0.5, 0.5, 0.7]] {
            let (axis, angle) = cam.turn_to_face(d);
            let m = Camera { orientation: cam.turned(axis, angle), ..cam }.basis();
            let n = unit(d);
            let z = dot(m[2], n);
            // Straight behind has no one axis to turn about; anything else faces.
            if dot(axis, axis) > 1e-6 {
                assert!(z > 0.999, "{d:?} ends {z}");
            }
        }
    }

    #[test]
    fn a_surface_is_picked_and_painted_by_residue() {
        let s = sample("1crn.pdb.gz");
        let cam = Camera::frame(&s);
        let look = Appearance::default();
        let base = scene(&s, Repr::Surface, Coloring::Chain, false);
        assert!(base.tris.len() > 10_000);
        assert_eq!(base.tri_residue.len(), base.tris.len());
        // Water is left out of the skin.
        assert!(base.tri_residue.iter().all(|&r| !is_solvent(&s.residues[r as usize].name)));

        // A point on the skin picks the residue under it.
        let hit = (100..300)
            .step_by(10)
            .find_map(|x| pick(&s, Repr::Surface, &look, &cam, 400, 400, x as f32, 200.0))
            .expect("nothing picked across the middle");
        assert!(hit.residue);
        let residue = s.atoms[hit.atom].residue;
        // Painted on the skin rather than glowing.
        let glow = build_highlight(&s, hit, Repr::Surface, Style::Soft, &look);
        assert!(glow.is_empty());
        let mut painted = base.clone();
        painted.paint_residue(0..base.tris.len(), 0..0, residue, Style::Soft, false);
        let changed = (0..base.tris.len()).filter(|&i| painted.tris[i].c != base.tris[i].c).count();
        assert!(changed > 0);
    }

    #[test]
    fn colours_follow_b_factors_and_secondary_structure() {
        let s = sample("1crn.pdb.gz");
        let (lo, hi) = s.b_span().expect("crambin has B-factors");
        assert!(hi > lo);
        let polymer: Vec<usize> = (0..s.residues.len()).filter(|&r| !s.residues[r].outside_polymer()).collect();
        let stillest = *polymer.iter().min_by(|&&a, &&b| s.residues[a].b.total_cmp(&s.residues[b].b)).unwrap();
        let floppiest = *polymer.iter().max_by(|&&a, &&b| s.residues[a].b.total_cmp(&s.residues[b].b)).unwrap();
        let p = Palette::Viridis;
        assert_eq!(residue_color(&s, stillest, Coloring::BFactor, p), p.ramp(0.0));
        assert_eq!(residue_color(&s, floppiest, Coloring::BFactor, p), p.ramp(1.0));

        let helix = polymer.iter().find(|&&r| s.residues[r].ss == Ss::Helix).unwrap();
        let strand = polymer.iter().find(|&&r| s.residues[r].ss == Ss::Sheet).unwrap();
        let (h, e) = (residue_color(&s, *helix, Coloring::Structure, p), residue_color(&s, *strand, Coloring::Structure, p));
        assert_ne!(h, e);
        // The palette does not move them.
        assert_eq!(h, residue_color(&s, *helix, Coloring::Structure, Palette::Magma));
    }

    #[test]
    fn a_flat_or_predicted_column_is_no_b_factor() {
        let mut predicted = sample("1crn.pdb.gz");
        predicted.has_confidence = true;
        assert_eq!(predicted.b_span(), None, "pLDDT taken for B-factors");
        let mut flat = sample("1crn.pdb.gz");
        for r in &mut flat.residues {
            r.b = 100.0;
        }
        assert_eq!(flat.b_span(), None);
    }

    #[test]
    fn a_transparent_frame_shows_the_page_only_round_the_molecule() {
        let s = sample("caffeine.sdf");
        let cam = Camera::frame(&s);
        let sc = scene(&s, Repr::Spacefill, Coloring::Element, false);
        let px = |rgba: &[u8], x: usize, y: usize| -> [u8; 4] { rgba[(y * 200 + x) * 4..][..4].try_into().unwrap() };
        // Without fog, the molecule is solid and the page clear.
        let rgba = render_transparent(&sc, &cam, 200, 200, Style::Illustrative);
        assert_eq!(px(&rgba, 2, 2)[3], 0, "the corner is not clear");
        assert!(px(&rgba, 100, 100)[3] > 250, "the middle is not solid");
        // With it, the far side fades into what is behind; laid on the page
        // it was drawn on, every pixel is as drawn there.
        let page = [0.97, 0.96, 0.94];
        let rgba = render_transparent(&sc, &cam, 200, 200, Style::Soft);
        let drawn = render(&sc, &cam, 200, 200, page, Style::Soft);
        for (x, y) in [(100, 100), (80, 120), (130, 90), (2, 2)] {
            let m = px(&rgba, x, y);
            let a = m[3] as f32 / 255.0;
            let d = &drawn.bgra[(y * 200 + x) * 4..][..3];
            for (k, want) in [d[2], d[1], d[0]].into_iter().enumerate() {
                let over = m[k] as f32 * a + page[k] * 255.0 * (1.0 - a);
                assert!((over - want as f32).abs() <= 3.0, "({x}, {y}) {over} for {want}");
            }
        }
    }

    #[test]
    fn a_slab_cuts_and_perspective_widens_the_near_side() {
        let s = sample("1crn.pdb.gz");
        let base = Camera::frame(&s);
        let sc = scene(&s, Repr::Spacefill, Coloring::Element, false);
        let page = [0.97, 0.96, 0.94];
        let covered = |cam: &Camera| {
            let f = render_cpu(&sc, cam, 300, 300, page, Style::Flat);
            f.bgra.chunks_exact(4).filter(|p| p[2] != (page[0] * 255.0 + 0.5) as u8 || p[1] != (page[1] * 255.0 + 0.5) as u8).count()
        };
        let whole = covered(&base);
        let r = base.radius;
        // A thin slab through the middle leaves less of it, and one out past
        // the structure leaves nothing.
        let thin = covered(&Camera { slab: Some([0.1 * r, -0.1 * r]), ..base });
        assert!(thin < whole * 9 / 10, "{thin} of {whole}");
        assert_eq!(covered(&Camera { slab: Some([3.0 * r, 2.5 * r]), ..base }), 0);
        // In perspective the near half of the structure, spread wider, covers
        // more than it does flat.
        let near = |perspective| covered(&Camera { perspective, slab: Some([2.0 * r, 0.0]), ..base });
        assert!(near(true) > near(false), "{} vs {}", near(true), near(false));
        // And a pick finds what is drawn: the nearest atom, at its warped place.
        let cam = Camera { perspective: true, ..base };
        let m = cam.basis();
        let view = |p: [f32; 3]| {
            let q = sub(p, cam.center);
            [dot(m[0], q), dot(m[1], q), dot(m[2], q)]
        };
        let front = (0..s.atoms.len()).max_by(|&a, &b| view(s.atoms[a].pos)[2].total_cmp(&view(s.atoms[b].pos)[2])).unwrap();
        let (v, _) = cam.warp(view(s.atoms[front].pos));
        let scale = 400.0 / (2.2 * cam.radius) * cam.zoom;
        let hit = pick(&s, Repr::Spacefill, &Appearance::default(), &cam, 400, 400, 200.0 + v[0] * scale, 200.0 - v[1] * scale);
        assert_eq!(hit.map(|h| h.atom), Some(front));
        // A slab in front of it hides it from the pick too.
        let cut = Camera { slab: Some([v[2] - 2.0, v[2] - 6.0]), ..cam };
        let hit = pick(&s, Repr::Spacefill, &Appearance::default(), &cut, 400, 400, 200.0 + v[0] * scale, 200.0 - v[1] * scale);
        assert_ne!(hit.map(|h| h.atom), Some(front));
    }

    #[test]
    fn a_layer_draws_its_atoms_and_hiding_takes_them_from_the_cartoon() {
        let s = sample("1crn.pdb.gz");
        let mask: Vec<bool> = s.atoms.iter().map(|a| a.residue < 5).collect();
        let n = mask.iter().filter(|x| **x).count();
        let look = Appearance::default();
        let spheres = build_layer(&s, &mask, LayerStyle::Spheres, None, Style::Soft, Palette::Spectrum, false, &look);
        assert_eq!(spheres.spheres.len(), n);
        let sticks = build_layer(&s, &mask, LayerStyle::Sticks, Some([1.0, 0.0, 0.0]), Style::Soft, Palette::Spectrum, false, &look);
        assert!(!sticks.capsules.is_empty());
        assert!(sticks.capsules.iter().all(|c| c.ca == sticks.capsules[0].ca), "not all in the one colour");
        // Hidden, those residues leave the ribbon.
        let whole = scene(&s, Repr::Cartoon, Coloring::Chain, false);
        let hidden = Appearance { hidden: Some(std::sync::Arc::new(mask.clone())), ..Appearance::default() };
        let cut = build_scene(&s, Repr::Cartoon, Coloring::Chain, 1.0, Style::Soft, Palette::Spectrum, false, false, &hidden);
        assert!(cut.tris.len() < whole.tris.len());
        assert!(cut.tri_residue.iter().all(|&r| r >= 5));
        let atoms = build_scene(&s, Repr::Spacefill, Coloring::Chain, 1.0, Style::Soft, Palette::Spectrum, false, false, &hidden);
        assert_eq!(atoms.spheres.len(), s.atoms.len() - n);
    }

    #[test]
    fn the_rainbow_runs_over_the_polymer_whatever_else_is_there() {
        // Haemoglobin's waters come after its chains in the file.
        let s = sample("4hhb.pdb.gz");
        assert!(s.residues.iter().any(|r| r.name == "HOH"));
        let polymer: Vec<usize> = (0..s.residues.len()).filter(|&r| !s.residues[r].outside_polymer()).collect();
        let p = Palette::Viridis;
        assert_eq!(residue_color(&s, polymer[0], Coloring::Rainbow, p), p.ramp(0.0));
        assert_eq!(residue_color(&s, *polymer.last().unwrap(), Coloring::Rainbow, p), p.ramp(1.0));
    }

    #[test]
    fn hydrophobic_and_charged_residues_and_the_potential_over_a_skin() {
        let s = sample("1crn.pdb.gz");
        let find = |name: &str| (0..s.residues.len()).find(|&r| s.residues[r].name == name).unwrap();
        let p = Palette::Spectrum;
        assert_eq!(residue_color(&s, find("ILE"), Coloring::Hydrophobic, p), HYDROPHOBIC);
        assert_eq!(residue_color(&s, find("ARG"), Coloring::Hydrophobic, p), HYDROPHILIC);
        assert_eq!(residue_color(&s, find("ARG"), Coloring::Charge, p), POSITIVE);
        assert_eq!(residue_color(&s, find("GLU"), Coloring::Charge, p), NEGATIVE);
        assert_eq!(residue_color(&s, find("THR"), Coloring::Charge, p), NEUTRAL);

        // Just off an arginine's guanidinium it is positive, off a
        // glutamate's carboxylate negative.
        let atom = |r: usize, name: &str| s.atoms[s.residue_atoms(r).into_iter().find(|&i| s.atoms[i].name() == name).unwrap()].pos;
        let off = |p: [f32; 3]| [p[0] + 1.5, p[1], p[2]];
        let phi = potentials(&s, &[off(atom(find("ARG"), "NH1")), off(atom(find("GLU"), "OE1"))]);
        assert!(phi[0] > 0.0 && phi[1] < 0.0, "{phi:?}");

        // Over the skin the colours run both ways.
        let sc = scene(&s, Repr::Surface, Coloring::Charge, false);
        let red = sc.tris.iter().filter(|t| t.c[0] > t.c[2] + 0.1).count();
        let blue = sc.tris.iter().filter(|t| t.c[2] > t.c[0] + 0.1).count();
        assert!(red > 100 && blue > 100, "{red} red, {blue} blue");
    }

    #[test]
    fn a_turn_between_two_views_lands_on_the_second() {
        let mut a = Camera { orientation: Camera::home(), zoom: 1.0, pan: [0.0; 2], center: [0.0; 3], radius: 10.0, perspective: false, slab: None };
        let mut b = a;
        b.rotate(1.2, -0.7);
        for (from, to) in [(a, b), (b, a)] {
            let (axis, angle) = from.turn_to(to.orientation);
            let m = from.turned(axis, angle);
            assert!((0..3).all(|i| (0..3).all(|j| (m[i][j] - to.orientation[i][j]).abs() < 1e-4)));
        }
        // Half a turn has no axis to read off the skew part; it still lands.
        a.orientation = Camera::home();
        let mut half = a;
        half.rotate(std::f32::consts::PI, 0.0);
        let (axis, angle) = a.turn_to(half.orientation);
        let m = a.turned(axis, angle);
        assert!((0..3).all(|i| (0..3).all(|j| (m[i][j] - half.orientation[i][j]).abs() < 1e-3)), "{m:?}");
    }

    #[test]
    fn angles_start_at_zero() {
        let mut cam = Camera::frame_all(&[]);
        assert_eq!(cam.angles(), [0.0, 0.0]);
        cam.rotate(std::f32::consts::FRAC_PI_2, 0.0);
        cam.rotate(0.0, -10f32.to_radians());
        let [x, y] = cam.angles();
        assert!((y - 90.0).abs() < 1e-3 && (x + 10.0).abs() < 1e-3, "{x} {y}");
        // A full turn comes back round.
        cam.rotate(std::f32::consts::TAU, 0.0);
        assert!((cam.angles()[1] - 90.0).abs() < 1e-2);
        cam.set_angles([30.0, -45.0]);
        let [x, y] = cam.angles();
        assert!((x - 30.0).abs() < 1e-3 && (y + 45.0).abs() < 1e-3, "{x} {y}");
    }

    #[test]
    fn the_view_tumbles_over_the_top_and_keeps_its_bearings() {
        let mut cam = Camera::frame_all(&[]);
        let step = 0.008 * 7.0;
        // Drag down, a long way: right over the top and back round, with
        // nowhere to stick. Every step turns the view by the same amount.
        let mut last = cam.basis();
        for _ in 0..200 {
            cam.rotate(0.0, step);
            let now = cam.basis();
            let turned = mul(now, transpose(last));
            let angle = ((turned[0][0] + turned[1][1] + turned[2][2] - 1.0) / 2.0).clamp(-1.0, 1.0).acos();
            assert!((angle - step).abs() < 1e-3, "a step turned {angle}");
            last = now;
        }
        // Upside down, a drag to the right still turns the near face right:
        // about the screen's vertical, whichever way the molecule is up.
        cam.set_angles([180.0, 0.0]);
        let before = cam.basis();
        cam.rotate(0.1, 0.0);
        assert_eq!(mul(cam.basis(), transpose(before))[1][1].round(), 1.0);
        let turned = mul(cam.basis(), transpose(before));
        assert!((turned[0][0] - 0.1f32.cos()).abs() < 1e-4 && turned[0][2] < 0.0);
        // And it stays a rotation.
        for _ in 0..10_000 {
            cam.rotate(0.013, -0.007);
        }
        let m = cam.basis();
        let i = mul(m, transpose(m));
        for (r, row) in i.iter().enumerate() {
            for (c, v) in row.iter().enumerate() {
                assert!((v - if r == c { 1.0 } else { 0.0 }).abs() < 1e-4);
            }
        }
    }

    #[test]
    fn translucent_spheres_blend_over_opaque_ones() {
        let behind = Sphere {
            c: [0.0, 0.0, -2.0],
            r: 1.0,
            color: [1.0, 0.0, 0.0],
            alpha: 1.0,
        };
        let front = Sphere {
            c: [0.0, 0.0, 2.0],
            r: 1.0,
            color: [0.0, 0.0, 1.0],
            alpha: 0.5,
        };
        let cam = Camera {
            orientation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            zoom: 1.0,
            pan: [0.0, 0.0],
            center: [0.0; 3],
            radius: 2.0,
            perspective: false,
            slab: None,
        };
        let draw = |spheres: Vec<Sphere>| {
            let scene = Scene {
                spheres,
                ..Scene::default()
            };
            render(&scene, &cam, 64, 64, [1.0; 3], Style::Soft)
        };
        let both = pixel(&draw(vec![behind, front]), 32, 32);
        // Red shows through the blue in front of it.
        assert!(both[0] > 40 && both[2] > 40, "{both:?}");
        // Fully transparent is not drawn at all.
        let mut gone = front;
        gone.alpha = 0.0;
        let alone = pixel(&draw(vec![behind]), 32, 32);
        assert_eq!(pixel(&draw(vec![behind, gone]), 32, 32), alone);
    }
}

/// Frame times for a structure, run with
/// `PROFOLD_BENCH=path cargo test --release bench_render -- --ignored --nocapture`.
#[cfg(test)]
mod bench {
    use super::*;
    use std::time::Instant;

    /// The card alone, long enough for its clocks to settle, at the size of
    /// the app's viewport on a retina screen.
    #[test]
    #[ignore]
    fn bench_gpu() {
        let path = std::env::var("PROFOLD_BENCH").expect("PROFOLD_BENCH names a structure");
        let s = crate::formats::load_path(std::path::Path::new(&path)).unwrap();
        let mut glass = Appearance::default();
        if let Some(kind) = s.kinds.iter().find(|k| k.as_str() == "6") {
            glass.set(Species::Kind(kind.clone()), Look { alpha: 0.02, ..Look::default() });
        }
        for (label, repr, appearance) in [
            ("spacefill", Repr::Spacefill, &Appearance::default()),
            ("spacefill+glass", Repr::Spacefill, &glass),
            ("ball&stick", Repr::BallAndStick, &Appearance::default()),
        ] {
            let scene = build_scene(
                &s, repr, Coloring::Element, 1.0, Style::Illustrative, Palette::Spectrum,
                false, true, appearance,
            )
            .stamp();
            let mut cam = Camera::frame(&s);
            let n = 300;
            for _ in 0..30 {
                render(&scene, &cam, 1616, 1436, [1.0; 3], Style::Illustrative);
            }
            let t = Instant::now();
            for _ in 0..n {
                cam.rotate(0.02, 0.0);
                let frame = render(&scene, &cam, 1616, 1436, [1.0; 3], Style::Illustrative);
                std::hint::black_box(frame);
            }
            println!("{label}: {:.2?}/frame", t.elapsed() / n);
        }
    }

    #[test]
    #[ignore]
    fn bench_render() {
        let path = std::env::var("PROFOLD_BENCH").expect("PROFOLD_BENCH names a structure");
        let t = Instant::now();
        let s = crate::formats::load_path(std::path::Path::new(&path)).unwrap();
        println!("load: {:?}, {} atoms, {} bonds", t.elapsed(), s.atoms.len(), s.bonds.len());
        // As well as everything opaque, the solvent turned almost to glass,
        // which is how a simulation box is usually looked into.
        let mut glass = Appearance::default();
        if let Some(kind) = s.kinds.iter().find(|k| k.as_str() == "6") {
            glass.set(Species::Kind(kind.clone()), Look { alpha: 0.02, ..Look::default() });
        }
        for (repr, appearance) in [
            (Repr::Spacefill, &Appearance::default()),
            (Repr::Spacefill, &glass),
            (Repr::BallAndStick, &Appearance::default()),
        ] {
            for style in Style::ALL {
                let t = Instant::now();
                let scene = build_scene(
                    &s, repr, Coloring::Element, 1.0, style, Palette::Spectrum, false, true,
                    appearance,
                )
                .stamp();
                let built = t.elapsed();
                let mut cam = Camera::frame(&s);
                for (w, h) in [(900, 700), (1800, 1400)] {
                    let n = 8;
                    let time = |draw: &dyn Fn(&Camera) -> Frame, cam: &mut Camera| {
                        draw(cam);
                        let t = Instant::now();
                        for i in 0..n {
                            cam.rotate(0.01 * i as f32, 0.0);
                            std::hint::black_box(draw(cam));
                        }
                        t.elapsed() / n
                    };
                    let cpu = time(&|c| render_cpu(&scene, c, w, h, [1.0; 3], style), &mut cam);
                    let gpu = time(&|c| render(&scene, c, w, h, [1.0; 3], style), &mut cam);
                    println!(
                        "{:?}{} {:?} build {:?}  {w}x{h}: cpu {:?}  gpu {:?}",
                        repr,
                        if appearance.is_default() { "" } else { "+glass" },
                        style, built, cpu, gpu
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod gpu_tests {
    use super::*;

    /// The card draws what the CPU draws: past the odd pixel along an edge,
    /// where the two round a pixel centre differently, the frames agree.
    #[test]
    fn gpu_matches_cpu() {
        let Some(gpu) = crate::gpu::get() else {
            eprintln!("no GPU here; skipped");
            return;
        };
        let samples = [
            ("caffeine.sdf", include_bytes!("../samples/caffeine.sdf").to_vec()),
            ("1crn.pdb.gz", include_bytes!("../samples/1crn.pdb.gz").to_vec()),
            ("rocksalt.cif", include_bytes!("../samples/rocksalt.cif").to_vec()),
        ];
        let out = std::env::var("PROFOLD_DIFF_DIR").ok();
        for (name, bytes) in samples {
            let s = crate::formats::load_bytes(name, bytes).unwrap();
            for repr in Repr::ALL {
                if repr.needs_polymer() && !s.has_polymer() {
                    continue;
                }
                for style in Style::ALL {
                    for see_through in [false, true] {
                        let mut appearance = Appearance::default();
                        if see_through {
                            // Half the elements turned to glass.
                            for (i, (e, _)) in s.composition().into_iter().enumerate() {
                                if i % 2 == 0 {
                                    appearance.set(
                                        Species::Element(e.0),
                                        Look { alpha: 0.5, ..Look::default() },
                                    );
                                }
                            }
                        }
                        let scene = build_scene(
                            &s, repr, Coloring::Element, 1.0, style, Palette::Spectrum,
                            false, true, &appearance,
                        )
                        .stamp();
                        let mut cam = Camera::frame(&s);
                        cam.rotate(0.3, 0.0);
                        // Lines one, two and three pixels wide, which the card
                        // finds by three routes.
                        let sizes: &[(usize, usize)] = if name == "caffeine.sdf" {
                            &[(640, 480), (1500, 1450), (2200, 2150)]
                        } else {
                            &[(640, 480)]
                        };
                        for &(w, h) in sizes {
                        let cpu = render_cpu(&scene, &cam, w, h, [0.97, 0.96, 0.94], style);
                        let gpu = gpu
                            .render(&scene, &cam, w, h, [0.97, 0.96, 0.94], style)
                            .unwrap();
                        // Glass is blended into eight bits a layer at a time on
                        // the card, in floating point on the CPU, so stacked
                        // layers drift a few levels apart.
                        let tolerance = if see_through { 16 } else { 6 };
                        let off = cpu
                            .bgra
                            .chunks_exact(4)
                            .zip(gpu.bgra.chunks_exact(4))
                            .filter(|(a, b)| a.iter().zip(*b).any(|(x, y)| x.abs_diff(*y) > tolerance))
                            .count();
                        let frac = off as f32 / (w * h) as f32;
                        let label = format!("{name}-{repr:?}-{style:?}-{see_through}-{w}");
                        println!("{label}: {:.3}% differ", frac * 100.0);
                        if let Some(dir) = &out {
                            for (tag, f) in [("cpu", &cpu), ("gpu", &gpu)] {
                                let mut rgba = f.bgra.clone();
                                rgba.chunks_exact_mut(4).for_each(|p| p.swap(0, 2));
                                image::save_buffer(
                                    format!("{dir}/{label}-{tag}.png"),
                                    &rgba, w as u32, h as u32, image::ColorType::Rgba8,
                                )
                                .unwrap();
                            }
                        }
                        assert!(frac < 0.01, "{label}: {:.2}% of pixels differ", frac * 100.0);
                        }
                    }
                }
            }
        }
    }
}

/// The pictures in the README, redrawn with
/// `cargo test --release readme_gallery -- --ignored`; they land in `docs/`.
#[cfg(test)]
mod gallery {
    use super::*;

    #[test]
    fn gpu_matches_cpu_in_perspective_through_a_slab() {
        let Some(gpu) = crate::gpu::get() else {
            eprintln!("no GPU here; skipped");
            return;
        };
        let s = crate::formats::load_bytes("1crn.pdb.gz", include_bytes!("../samples/1crn.pdb.gz").to_vec()).unwrap();
        for repr in [Repr::Cartoon, Repr::BallAndStick, Repr::Spacefill, Repr::Surface] {
            for style in [Style::Illustrative, Style::Soft] {
                let scene = build_scene(&s, repr, Coloring::Chain, 1.0, style, Palette::Spectrum, false, false, &Appearance::default()).stamp();
                let mut cam = Camera::frame(&s);
                cam.rotate(0.3, 0.2);
                cam.perspective = true;
                cam.slab = Some([0.4 * cam.radius, -0.3 * cam.radius]);
                let page = [0.97, 0.96, 0.94];
                let cpu = render_cpu(&scene, &cam, 640, 480, page, style);
                let gpu = gpu.render(&scene, &cam, 640, 480, page, style).unwrap();
                let off = cpu
                    .bgra
                    .chunks_exact(4)
                    .zip(gpu.bgra.chunks_exact(4))
                    .filter(|(a, b)| a.iter().zip(*b).any(|(x, y)| x.abs_diff(*y) > 6))
                    .count();
                let frac = off as f32 / (640.0 * 480.0);
                assert!(frac < 0.01, "{repr:?} {style:?}: {:.2}% of pixels differ", frac * 100.0);
            }
        }
    }

    #[test]
    #[ignore]
    fn readme_gallery() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("docs");
        std::fs::create_dir_all(&dir).unwrap();
        let shots: [(&str, &str, Repr, Coloring, Style, Palette, f32, bool); 14] = [
            ("haemoglobin", "4hhb.pdb.gz", Repr::Cartoon, Coloring::Chain, Style::Illustrative, Palette::Spectrum, 0.0, false),
            ("gfp", "1ema.pdb.gz", Repr::Cartoon, Coloring::Rainbow, Style::Illustrative, Palette::Spectrum, 0.0, false),
            ("ubiquitin-soft", "1ubq.pdb.gz", Repr::Cartoon, Coloring::Rainbow, Style::Soft, Palette::Spectrum, 0.0, false),
            ("crambin-tube", "1crn.pdb.gz", Repr::Tube, Coloring::Rainbow, Style::Illustrative, Palette::Viridis, 0.0, false),
            ("crambin-ballstick", "1crn.pdb.gz", Repr::BallAndStick, Coloring::Element, Style::Illustrative, Palette::Spectrum, 0.0, false),
            ("haemoglobin-spacefill", "4hhb.pdb.gz", Repr::Spacefill, Coloring::Chain, Style::Soft, Palette::Spectrum, 0.0, false),
            ("caffeine", "caffeine.sdf", Repr::BallAndStick, Coloring::Element, Style::Illustrative, Palette::Spectrum, 0.0, false),
            ("c60", "c60.xyz", Repr::Licorice, Coloring::Element, Style::Soft, Palette::Plasma, 0.0, false),
            ("rocksalt", "rocksalt.cif", Repr::BallAndStick, Coloring::Element, Style::Illustrative, Palette::Spectrum, 0.0, true),
            ("silicon", "silicon.vasp", Repr::BallAndStick, Coloring::Element, Style::Soft, Palette::Ice, 0.0, true),
            ("melt", "lj-melt.lammpstrj.gz", Repr::Spacefill, Coloring::Element, Style::Illustrative, Palette::Spectrum, 1.0, true),
            ("ubiquitin-viridis", "1ubq.pdb.gz", Repr::Cartoon, Coloring::Rainbow, Style::Illustrative, Palette::Viridis, 0.0, false),
            ("ubiquitin-magma", "1ubq.pdb.gz", Repr::Cartoon, Coloring::Rainbow, Style::Illustrative, Palette::Magma, 0.0, false),
            ("ubiquitin-mako", "1ubq.pdb.gz", Repr::Cartoon, Coloring::Rainbow, Style::Illustrative, Palette::Mako, 0.0, false),
        ];
        for (name, file, repr, coloring, style, palette, at, cell) in shots {
            let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("samples").join(file);
            let mut s = crate::formats::load_path(&path).unwrap();
            // `at` is how far through a trajectory to draw it, 0..=1.
            let last = s.frame_count() - 1;
            s.set_frame((at * last as f32).round() as usize);
            let scene = build_scene(
                &s, repr, coloring, 1.0, style, palette, false, cell, &Appearance::default(),
            )
            .stamp();
            let cam = Camera::frame(&s);
            // Drawn large, then cut to what was drawn and centred on a 4:3
            // page of its own, so every picture sits the same way.
            let (w, h) = (2000, 2000);
            let frame = render(&scene, &cam, w, h, crate::theme::view(), style);
            let mut rgba = frame.bgra;
            for px in rgba.chunks_exact_mut(4) {
                px.swap(0, 2);
            }
            let img = image::RgbaImage::from_raw(w as u32, h as u32, rgba).unwrap();
            let paper = *img.get_pixel(0, 0);
            let (mut x0, mut y0, mut x1, mut y1) = (w as u32, h as u32, 0, 0);
            for (x, y, p) in img.enumerate_pixels() {
                if p.0.iter().zip(paper.0).any(|(a, b)| a.abs_diff(b) > 4) {
                    (x0, y0, x1, y1) = (x0.min(x), y0.min(y), x1.max(x), y1.max(y));
                }
            }
            let (cw, ch) = (x1 - x0 + 1, y1 - y0 + 1);
            let side = ((cw as f32 * 1.12).max(ch as f32 * 1.12 * 4.0 / 3.0)) as u32;
            let (pw, ph) = (side, side * 3 / 4);
            let mut page = image::RgbaImage::from_pixel(pw, ph, paper);
            let cut = image::imageops::crop_imm(&img, x0, y0, cw, ch).to_image();
            image::imageops::overlay(&mut page, &cut, ((pw - cw) / 2) as i64, ((ph - ch) / 2) as i64);
            image::imageops::resize(&page, 1200, 900, image::imageops::FilterType::Lanczos3)
                .save(dir.join(format!("{name}.png")))
                .unwrap();
        }
    }

    /// Every style over a spread of structures, on both pages, side by side:
    /// `PROFOLD_SHEET=dir cargo test --release contact_sheet -- --ignored`.
    #[test]
    #[ignore]
    fn contact_sheet() {
        let Some(dir) = std::env::var_os("PROFOLD_SHEET") else { return };
        let dir = std::path::PathBuf::from(dir);
        std::fs::create_dir_all(&dir).unwrap();
        let rows: [(&str, &str, Repr, Coloring, Palette, bool); 9] = [
            ("haemoglobin", "4hhb.pdb.gz", Repr::Cartoon, Coloring::Chain, Palette::Spectrum, false),
            ("gfp", "1ema.pdb.gz", Repr::Cartoon, Coloring::Rainbow, Palette::Turbo, false),
            ("ubq-magma", "1ubq.pdb.gz", Repr::Cartoon, Coloring::Rainbow, Palette::Magma, false),
            ("ubq-mako", "1ubq.pdb.gz", Repr::Tube, Coloring::Rainbow, Palette::Mako, false),
            ("crambin", "1crn.pdb.gz", Repr::BallAndStick, Coloring::Element, Palette::Spectrum, false),
            ("hb-fill", "4hhb.pdb.gz", Repr::Spacefill, Coloring::Chain, Palette::Greyscale, false),
            ("caffeine", "caffeine.sdf", Repr::BallAndStick, Coloring::Element, Palette::Spectrum, false),
            ("rocksalt", "rocksalt.cif", Repr::BallAndStick, Coloring::Element, Palette::Spectrum, true),
            ("silicon", "silicon.vasp", Repr::BallAndStick, Coloring::Element, Palette::Ice, true),
        ];
        let tile = 360u32;
        for dark in [false, true] {
            crate::theme::set_dark(dark);
            let mut sheet =
                image::RgbaImage::new(tile * Style::ALL.len() as u32, tile * rows.len() as u32);
            for (r, (_, file, repr, coloring, palette, cell)) in rows.iter().enumerate() {
                let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("samples").join(file);
                let s = crate::formats::load_path(&path).unwrap();
                for (c, style) in Style::ALL.iter().enumerate() {
                    let scene = build_scene(
                        &s, *repr, *coloring, 1.0, *style, *palette, false, *cell,
                        &Appearance::default(),
                    )
                    .stamp();
                    let cam = Camera::frame(&s);
                    let side = tile as usize * 2;
                    let frame = render(&scene, &cam, side, side, crate::theme::view(), *style);
                    let mut rgba = frame.bgra;
                    for px in rgba.chunks_exact_mut(4) {
                        px.swap(0, 2);
                    }
                    let img = image::RgbaImage::from_raw(side as u32, side as u32, rgba).unwrap();
                    let small = image::imageops::resize(&img, tile, tile, image::imageops::FilterType::Lanczos3);
                    image::imageops::overlay(&mut sheet, &small, (c as u32 * tile) as i64, (r as u32 * tile) as i64);
                }
            }
            let name = if dark { "dark.png" } else { "light.png" };
            sheet.save(dir.join(name)).unwrap();
        }
        crate::theme::set_dark(false);
    }
}
