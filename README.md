<p align="center">
  <img src="assets/logo.svg" width="120" alt="profold logo">
</p>

<h1 align="center">profold</h1>

<p align="center">
  A small, fast molecular viewer for proteins, molecules, crystals and simulations,<br>
  written in Rust with <a href="https://www.gpui.rs">GPUI</a> and drawn on the GPU.
</p>

<p align="center">
  <b><a href="#install">Install</a></b> · macOS · Linux · Windows · Apple silicon, Intel/AMD and ARM
</p>

<p align="center">
  <img src="docs/app-haemoglobin.png" alt="profold showing human deoxyhaemoglobin (4HHB) as a cartoon coloured by chain">
</p>

Open a structure or a trajectory and turn it around: a protein from the PDB, a
ligand, a crystal, or a LAMMPS or GROMACS simulation box. Search for a protein
or a molecule by name, or paste an identifier, and it fetches the structure.
Each tab holds its own, so a protein, a trajectory and caffeine can be open side
by side. Pick atoms and measure between them, look into a ligand's binding site,
label what matters, analyse a trajectory frame by frame, and save a ray-traced
picture or a film of it. With the SimpleFold plugin, a fold tab folds a sequence with
[Apple's SimpleFold](https://github.com/apple/ml-simplefold), then lays the
prediction over the experiment and tells you how close it came.

## Install

**macOS and Linux**, in a terminal:

```bash
curl -LsSf https://github.com/odexion/profold/releases/latest/download/install.sh | sh
```

**Windows**, in PowerShell:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://github.com/odexion/profold/releases/latest/download/install.ps1 | iex"
```

That's all. The installer picks the right build for your machine, shows a
progress bar while it downloads, checks the download against its SHA-256 and
installs profold for your user only, with no `sudo` and no administrator rights:

| | Installs to | Start it from |
| --- | --- | --- |
| macOS | `~/Applications/profold.app`, plus `profold` in `~/.local/bin` | Applications, Spotlight, or `profold` |
| Linux | `~/.local/bin/profold`, plus a launcher and icon | the applications menu, or `profold` |
| Windows | `%LOCALAPPDATA%\Programs\profold`, on your `PATH` | the Start menu, or `profold` |

Builds are made for **Apple silicon and Intel Macs**, **x86-64 and ARM64
Linux** (glibc 2.35 or newer, such as Ubuntu 22.04+, Debian 12+, Fedora 36+)
and **x86-64 and ARM64 Windows 10/11**. On a Mac, installing this way never
shows the "unidentified developer" warning, because nothing downloaded by
`curl` is quarantined.

**Updates.** profold checks for a new release when it starts (at most once a
day), downloads it in the background and shows **Restart to update** in the tab
bar. The version is at the right of the tab bar; click it to check for updates
yourself or to turn the automatic check off. To install a particular version,
set `PROFOLD_VERSION=0.2.0` before running the installer.

<details>
<summary>Uninstalling</summary>

- **macOS:** delete `~/Applications/profold.app` and `~/.local/bin/profold`.
- **Linux:** delete `~/.local/bin/profold`,
  `~/.local/share/applications/dev.profold.desktop` and
  `~/.local/share/icons/hicolor/256x256/apps/dev.profold.png`.
- **Windows:** delete `%LOCALAPPDATA%\Programs\profold` and the *profold*
  shortcut in the Start menu, and remove that folder from your user `PATH`.

Settings, the folding environment and its weights live in `~/.profold` (and
`~/.cache/torch`); delete those too to remove everything.
</details>

## Highlights

- **Search by name.** Type *green fluorescent protein* or *caffeine* and pick
  from structures in the PDB, AlphaFold models and molecules in PubChem, best
  match first. A gene symbol works too: *TP53* brings up human p53. No
  identifier needed.
- **One field for everything.** A name, a path, a PDB ID, a UniProt accession,
  a PubChem CID or a sequence. It works out which one you gave it.
- **Tabs.** Each tab has its own structure, camera and job, so a fold or a
  download in one never holds up the others.
- **Many file formats.** PDB, mmCIF and BinaryCIF, crystallographic CIF, XYZ,
  GROMACS, MOL2, SDF, LAMMPS data and dump, VASP, XCrySDen, DCD and XTC
  trajectories, and volumes from Gaussian cube, CCP4/MRC and OpenDX. Gzipped
  files open as they are.
- **Cartoons without annotations.** Helices and sheets are assigned from the
  backbone geometry with a DSSP-style hydrogen-bond analysis, so even a bare
  prediction gets ribbons and arrows.
- **Six representations.** Cartoon, Tube, Ball & stick, Licorice, Spacefill and
  a molecular Surface.
- **Four looks.** *Illustrative* has ink outlines and banded light, like a drawn
  figure. *Flat* lays soft colour almost flat inside the ink, after David
  Goodsell. *Soft* has smooth shading and depth fog. *Studio* has glossy light
  and shadow, for a paper or a cover. Light and dark mode for the window.
- **Colour by** chain, position, B-factor, secondary structure, hydrophobicity,
  charge (the electrostatic potential over a surface), amino acid, pLDDT,
  element or any per-atom property from the file, with fifteen palettes.
- **Look closer.** Hover to name an atom, click to pick it, measure distances,
  angles and dihedrals, and open a ligand's binding site with its hydrogen
  bonds, salt bridges, metal contacts and clashes. Select with a small
  language (`resn HIS and within 5 of ligand`, `c_pe < -3`), and pin labels to
  atoms. A Ramachandran plot, and AlphaFold's predicted aligned error.
- **Trajectories.** Every multi-frame format plays, steps and scrubs, and
  periodic cells are drawn as a fine box. Analysis plots RMSD, radius of
  gyration, SASA, polar contacts, RMSF and secondary structure through time.
- **Simulation boxes.** Wrap atoms into the cell, mend molecules the walls cut,
  unwrap a trajectory in time, or replicate the cell. Colour and select by any
  column of a LAMMPS dump or extended XYZ.
- **Isosurfaces** of densities and orbitals, solid or as a mesh, with a second
  lobe for negative values.
- **Per-element and per-type styling.** Size, opacity and colour sliders for
  each species. Turn the solvent to glass and look inside the box.
- **Fast.** Impostor spheres and capsules on the GPU through wgpu, off the main
  thread. Without a GPU it falls back to a parallel software rasterizer that
  draws the same image.
- **Folding in its own tab** (plugin). Opt in from the + menu.
- **Prediction vs. experiment** (plugin). Fold a PDB entry's own sequence,
  superimpose it with Kabsch, and read off the Cα RMSD, coloured by deviation.
- **Publication images.** Save the viewport as a PNG at one to three times the screen size,
  ray traced if you like, with soft shadows and ambient occlusion, on the
  graphics card where it can.
- **Films.** Play a trajectory while the camera turns, rocks, or follows a path
  you recorded by hand, with smoothing between frames. Export as GIF or PNG
  frames, or with ffmpeg as MP4 or transparent WebM or ProRes, at up to 4K, with
  a frame counter, simulation time, a scale bar, the colour legend and a caption.
- **A console** (⌘K) for every command, `profold --run script.txt` to run a
  script of them, undo and redo, and sessions that reopen as they were left.

## Gallery

### Proteins

| | | |
| :---: | :---: | :---: |
| <img src="docs/haemoglobin.png" alt="Haemoglobin cartoon by chain"> | <img src="docs/gfp.png" alt="GFP cartoon in rainbow"> | <img src="docs/ubiquitin-soft.png" alt="Ubiquitin cartoon, soft style"> |
| **Haemoglobin** (4HHB): cartoon by chain, haems as atoms | **Green fluorescent protein** (1EMA): the β-barrel and its chromophore | **Ubiquitin** (1UBQ): the β-grasp fold in the *Soft* style |
| <img src="docs/haemoglobin-spacefill.png" alt="Haemoglobin spacefill"> | <img src="docs/crambin-tube.png" alt="Crambin tube"> | <img src="docs/crambin-ballstick.png" alt="Crambin ball and stick"> |
| **Haemoglobin**: spacefill, four chains packed together | **Crambin** (1CRN): a backbone tube | **Crambin**: every atom in ball & stick |
| <img src="docs/haemoglobin-surface.png" alt="Haemoglobin surface coloured by electrostatic potential"> | <img src="docs/gfp-flat.png" alt="GFP cartoon in the Flat style"> | <img src="docs/ubiquitin-amino.png" alt="Ubiquitin licorice coloured by amino acid"> |
| **Haemoglobin**: its surface, coloured by the electrostatic potential | **Green fluorescent protein** in the *Flat* style | **Ubiquitin**: licorice, each amino acid in RasMol's colour |

### Ray traced

| | |
| :---: | :---: |
| <img src="docs/haemoglobin-traced.png" alt="Haemoglobin spacefill, ray traced in the Studio style"> | <img src="docs/gfp-traced.png" alt="GFP cartoon, ray traced in the Studio style"> |
| **Haemoglobin**: spacefill in the *Studio* style, with shadows and ambient occlusion | **Green fluorescent protein**: the chromophore inside its β-barrel, in the *Studio* style |

### Molecules, crystals and simulations

| | | |
| :---: | :---: | :---: |
| <img src="docs/caffeine.png" alt="Caffeine ball and stick"> | <img src="docs/c60.png" alt="C60 licorice"> | <img src="docs/rocksalt.png" alt="Rock salt unit cell"> |
| **Caffeine**: ball & stick, from an SDF file | **Buckminsterfullerene**: licorice, from XYZ | **Rock salt**: grown from its symmetry operations in a CIF |
| <img src="docs/silicon.png" alt="Silicon supercell"> | <img src="docs/melt.png" alt="Lennard-Jones melt"> | <img src="docs/app-melt.png" alt="Trajectory player in the app"> |
| **Silicon**: a diamond-cubic 3×3×3 supercell from VASP | **Lennard-Jones melt**: the last frame of a LAMMPS dump | **Trajectories** play, step and scrub along the foot of the viewport |

### Palettes

| Viridis | Magma | Mako |
| :---: | :---: | :---: |
| <img src="docs/ubiquitin-viridis.png" alt="Ubiquitin in Viridis"> | <img src="docs/ubiquitin-magma.png" alt="Ubiquitin in Magma"> | <img src="docs/ubiquitin-mako.png" alt="Ubiquitin in Mako"> |

Spectrum, Turbo, Cividis, Inferno, Plasma, Rocket, Ember, Ice, Forest, Dusk,
Pastel and Greyscale are there too.

## Using it

<img src="docs/app-start.png" alt="The empty window, listing samples and example identifiers" align="right" width="420">

One field takes everything, and the app works out what you gave it:

| Input | What happens |
| --- | --- |
| `green fluorescent protein`: a name | Searches the PDB, AlphaFold DB and PubChem |
| `TP53`: a gene symbol | Finds the protein it codes for, human first, as structures and a model |
| `~/sims/dump.melt`: a path | Opens the file |
| `1CRN`: a PDB ID | Downloads the experimental structure from RCSB |
| `P69905`: a UniProt accession | Downloads the predicted model from AlphaFold DB |
| `CID 2519`: a PubChem compound | Downloads its 3D conformer from PubChem |
| `MVLSPADKTNVK…`: residue letters, FASTA and all | Finds PDB entries with that sequence |

<img src="docs/app-search.png" alt="Searching for green fluorescent protein: PDB structures, AlphaFold models and PubChem molecules, best match first">

A search lists what it finds as you type, grouped by database: experimental
structures, predicted models, and molecules. Return opens the best match, the
arrow keys or a click another one. The databases rank by their own relevance,
which for a protein favours the latest of a thousand variants, so profold ranks
the hits again by how well their names match, with a nudge towards the first
structure ever deposited (1UBQ for ubiquitin, 1MBN for myoglobin) and towards
human proteins. A name that is a molecule's, like *caffeine* or *aspirin*, puts
PubChem first. A compound without a 3D conformer, a peptide drug say, opens as
PubChem's flat drawing.

A gene symbol, *TP53*, *BRCA1* or *EGFR*, is looked up in UniProt's gene names
too. The protein it codes for comes up as its AlphaFold model, and the PDB
entries are the ones whose chains are that protein, ranked by how much of it
they hold, so p53's core domain comes before a peptide of it and KRAS's own
structures before its promoter's DNA. The human gene wins when several species
share the symbol; a word that only some bacterium's gene answers to, like
*heme*, stays the molecule.

**Tabs.** ⌘T opens a tab, ⇧⌘T a fold tab, ⌘W closes one, and ⌘1 to ⌘9 and
⌃Tab move between them. Drag a tab along the bar to reorder it. The tabs you
had open come back at the next launch, reopened from their files and
identifiers (a prediction is not kept, so a fold tab comes back empty).
Whatever is opened, searched or dropped on the window goes into the tab in front
and replaces what was there, as in a browser. Each tab keeps its own camera and
trajectory, which pauses when the tab goes behind another. A job in one tab
carries on while you use another, and its tab shows the progress; closing the
tab stops it.

In the viewport: drag to orbit, shift-drag to pan, scroll to zoom, double-click
to reframe. The molecule follows the pointer, so dragging right brings the near
face right. The angle and zoom readouts in the title bar can be clicked and
typed into. ⌘-drag rolls the view (Ctrl-drag off the Mac), ⌘-click flies to
an atom, and small axes in the corner say which way x, y and z point. The bar
along the top is **Open · Samples · Camera · Image · Film · Save**: **Camera**
centres, focuses, spins, rocks, and switches perspective and the clipping slab;
**Image** saves the molecule alone as a PNG at up to three times the viewport
size, transparent if you like. **Image › Ray traced** lights the still view with
rays, for soft shadows, ambient occlusion and smooth edges, and saved images
with it. **Film › Export film…** (⇧⌘E) makes an animation of it (see
[Pictures and films](#pictures-and-films)). **Save** writes the structure, a
chain or a binding site as PDB or mmCIF, or the whole session as a `.profold`
file. ⌘K opens the console, where every command can be typed, and ⌘/ lists the
shortcuts.

<br clear="right">

### Samples

No data to hand? Ten structures are built into the app, listed on the empty
screen and under **Samples** at any time:

- **Molecules:** caffeine and buckminsterfullerene
- **Proteins:** crambin, ubiquitin, green fluorescent protein, haemoglobin and
  SAM synthetase
- **Crystals and simulations:** rock salt grown from its symmetry, a silicon
  supercell, and a Lennard-Jones crystal melting over 41 frames

Where they come from is in [`samples/`](samples/README.md).

## Viewing

**Representations.** *Cartoon* is the default for a protein: ribbons for
helices, arrows for strands, cord for everything else. Ligands, ions and nucleic
acids are drawn as atoms beside it rather than left out. *Tube* is a plain
backbone trace. *Ball & stick* is the default for anything else, and *Licorice*
draws bonds as rods of one thickness. *Spacefill* puts every atom at its van der
Waals radius, and is how a box of anonymous beads opens. *Surface* is a smooth
Gaussian skin over the molecule, with its ligands standing inside it as atoms.
DNA and RNA are drawn as a tube along the backbone with a rung to each base.

**Styles.** *Illustrative* is the default and has the drawn look: ink outlines
on both sides of the silhouette and across every overlap, saturated colour,
light quantised into a few bands, and no depth fade, because a drawn figure does
not have one. *Flat* lays soft colour almost flat inside the ink, the look of a
painted molecular illustration. *Soft* has smooth shading on paper, with
distance washing out. *Studio* has glossy light and shadow where atoms crowd
together, without ink. The line weight scales with the frame, so it holds on a
retina display.

**Colour** by chain (by molecule, without chains), by position along the chain,
by B-factor, by secondary structure, by hydrophobicity (Kyte and Doolittle), by
charge, by amino acid (RasMol's colours, or the palette's), by pLDDT, by
element, or by a per-atom property. Over a surface, *Charge* shows the
electrostatic potential. Predicted models open coloured by pLDDT. The palette
colours the chains and the ramps; elements wear the standard colours (carbon
grey, oxygen red, nitrogen blue) unless *Element colours* is set to *Palette*.
Ligands beside a protein keep their element colours whatever the protein is
coloured by, and any chain or element can be given a colour of its own. Where
atoms are drawn, *Carbons only* lets the colouring take the carbons while every
other element keeps its own colour.

<img src="docs/app-dark-surface.png" alt="Green fluorescent protein as a surface coloured by electrostatic potential, in dark mode">

**Dark mode** is the ☾ in the tab bar, and is remembered between launches. The
layout button beside it moves the side panel under the view, for tall or
phone-shaped windows.

**Atoms.** The side panel lists every element, or every atom type for a
simulation, with its own size, opacity and colour. Preferences outlive the
session in `~/.profold/settings.json`.

Anything slow, like a download, an install or a fold, takes over the viewport
with its own progress. Where the underlying tool reports a percentage the bar is
real. Where it does not, the bar says so by travelling instead.

## Looking closer

Hover over the molecule and a card names what is under the pointer: a
residue's code and full name, its fold, B-factor or pLDDT; water, ions, common
ligands and modified bases by their names. A click picks it out, painted white
on the ribbon, and its card has **Focus** and **Label**. The sequence, chain by
chain, is linked both ways with the view, and shift-click picks out a run in
either. Pick atoms in turn to measure distances, angles and dihedrals, drawn as
dashes over the molecule.

<img src="docs/app-binding-site.png" alt="The haem binding site of haemoglobin's α chain, its residues named and its contacts drawn">

**Binding sites.** *Site* beside a ligand flies into it: the residues within
5 Å drawn and named, with its hydrogen bonds, salt bridges, metal contacts and
clashes as coloured dashes. Above, the haem of haemoglobin's α chain, held
between its proximal His87 and distal His58.

**Selections** use a small language after PyMOL's and VMD's: `chain`, `resi`,
`resn`, `name`, `elem`, `type`, `within N of`, `byres`, comparisons on any value
(`b < 30`, `z > 10`, `c_pe < -3`), joined by `and`, `or` and `not`. Each
selection is drawn its own way over the rest, as sticks, spheres or hidden.

| | |
| :---: | :---: |
| <img src="docs/app-labels.png" alt="Ubiquitin coloured by amino acid, with Lys48, Lys63 and its C terminus labelled"> | <img src="docs/app-ramachandran.png" alt="Ubiquitin's Ramachandran plot in the side panel"> |
| **Labels** pinned to atoms follow them through a trajectory and into images and films: ubiquitin's Lys48, Lys63 and C terminus, coloured by amino acid | **Plots**: a Ramachandran plot, and for AlphaFold models the predicted aligned error, both linked to the view |

Labels come from the picked atom's card, or from the console: `label resi 10-20`,
`label atoms name CA`, `label ligand "the drug"` and `label clear`.

**Biological assemblies** are built from the deposited chains, and **Open ›
Overlay** lays any structure over the one shown, paired by sequence alignment,
with the RMSD and a deviation colouring. ⌘Z and ⇧⌘Z undo and redo how a
structure is shown.

## Opening files

Drop a file on the window, use **Open…** (⌘O), type its path, or name it on the
command line. The format is told from the name, `data.*`, `dump.*`, `CONTCAR`
and `XDATCAR` included, and from the content when the name says nothing.
Gzipped files open as they are.

| Format | Extensions | Notes |
| --- | --- | --- |
| PDB, PQR | `.pdb` `.ent` `.pqr` | Every MODEL, or VMD's END-separated frames, becomes a frame |
| mmCIF, BinaryCIF | `.cif` `.mmcif` `.bcif` | Further models become frames; assemblies built from the file |
| CIF (crystallographic) | `.cif` | Fractional sites expanded by the symmetry operations to fill the cell |
| XYZ, extended XYZ | `.xyz` `.extxyz` | Multi-frame; `Lattice=` and `Properties=` honoured |
| GROMACS | `.gro` | Multi-frame, triclinic boxes, any precision |
| MOL2 | `.mol2` | Bonds from the file; same-atom molecules become frames |
| MDL molfile / SD | `.sdf` `.mol` | V2000 and V3000; conformers become frames |
| LAMMPS data | `.data` `.lmp` `data.*` | Any atom style; elements from Masses; bonds from the file |
| LAMMPS dump | `.lammpstrj` `.dump` `dump.*` | Text dumps; scaled, unwrapped and triclinic all handled; extra columns kept as per-atom properties |
| DCD, XTC | `.dcd` `.xtc` | Trajectories, laid over a structure with as many atoms |
| VASP | `POSCAR` `CONTCAR` `XDATCAR` `.vasp` | Direct or Cartesian; XDATCAR steps become frames |
| Gaussian cube | `.cube` | The molecule, its grid box and an isosurface of the grid |
| CCP4, MRC, OpenDX | `.ccp4` `.map` `.mrc` `.dx` | Volumes, laid over the structure shown or on their own |
| XCrySDen | `.xsf` `.axsf` | Molecules, crystals and animations |

Where a file does not list bonds they are worked out from distance: closer than
the two covalent radii plus 0.4 Å, as VMD does. Metal-to-metal pairs never bond,
so a lattice of iron does not turn into scaffolding. A simulation's atoms often
have no element at all, such as a LAMMPS type with no mass to go on or a
Lennard-Jones bead in reduced units. Those are coloured by type and sized from
their own spacing, the way OVITO draws them.

A trajectory gets a bar along the foot of the viewport: play, step, or drag
through the frames. A simulation cell is drawn as a fine box, and can be
switched off with **Cell**.

## Simulations and volumes

<img src="docs/app-melt.png" alt="The Lennard-Jones melt unwrapped through time, with its RMSD plotted in the side panel">

**Trajectory analysis** in the side panel measures, frame by frame, the RMSD
from the first frame (each fitted onto it by Kabsch), the radius of gyration,
the solvent-accessible surface (Shrake–Rupley) and the polar contacts between
residues; the RMSF of each atom; and helices and strands through time. Each is
worked out only when picked, in the background, of the atoms chosen: Cα, the
backbone, atom types or elements, or `analyse rmsd of <selection>` in the
console. A line marks the frame shown, and a click on the plot goes to that
frame. Above, the Lennard-Jones melt unwrapped through time, its RMSD rising as
the crystal loosens.

| | |
| :---: | :---: |
| <img src="docs/app-property.png" alt="The melt's atoms coloured by a per-atom displacement column"> | <img src="docs/app-periodic.png" alt="The Lennard-Jones box replicated 2×2×2"> |
| **Per-atom properties**: every extra column of a LAMMPS dump or extended XYZ, frame by frame, to colour by or select on. Here each atom's displacement since the first frame | **Periodic boundaries**: *Wrap*, *Whole* (mend molecules the walls cut), *Unwrap* in time, and *Replicate*, under the cell's size. Here the box 2×2×2 |

A property colours along the palette from its 1st to its 99th percentile across
the whole trajectory, so a stray atom does not wash out the rest; `color
property c_pe -5 -2` sets the range by hand. The cell tools are `wrap`,
`whole`, `unwrap`, `replicate 3 3 1` and `original` in the console.

<img src="docs/app-isosurface.png" alt="The 3d z² orbital of hydrogen as two blue lobes and a red torus" align="right" width="420">

**Isosurfaces.** A Gaussian cube file's grid, or a CCP4, MRC or OpenDX map,
is drawn as a surface where its values cross a level: mean plus 2.5 σ to
begin with, − and + to step it, or `iso 0.05` in the console. Signed values,
such as an orbital's, get a second, red lobe where they cross the negative
level. Solid or as a mesh, in the view, traced pictures and films. On the
right, hydrogen's 3d<sub>z²</sub> orbital, computed onto a cube grid.

<br clear="right">

## Pictures and films

**Image** saves the molecule as a PNG at one to three times the size of the
view, on the page or transparent. With **Ray traced** on, the view is lit
with rays once it holds still, with soft shadows, ambient occlusion and edges
smoothed by up to nine rays a pixel, and saved pictures are traced too. Every
style keeps its own look. Tracing runs on the graphics card through a compute
shader where it can (Metal, Vulkan, DirectX 12), and on the CPU otherwise; a
protein at 1200 × 900 takes about a second.

<img src="docs/app-film.png" alt="The Export a film dialog, with a turning preview of the melt">

**Film › Export film…** (⇧⌘E) makes an animation: the trajectory played
through, while the camera holds still, turns, rocks, or replays a path you
recorded with **Record a camera path**. Pick the frames, smoothing between
them, the frame rate, the size (up to 4K), antialiasing and the background,
then the format: GIF or numbered PNGs made in the app, or, with ffmpeg
installed, MP4, and WebM or ProRes 4444 with transparency. Overlays add a
frame counter, the simulation time, a scale bar, the colour legend and a
caption. The preview plays the film as it will be made, and the console's
`movie out.mp4 turn=360 fps=30 counter scalebar` makes one without the dialog.

## Comparing

*Part of the SimpleFold plugin.* In a fold tab, an identifier can be folded or
compared. **Compare**
is the interesting one. It pulls the deposited structure, folds that structure's
*own* sequence with SimpleFold, superimposes the two, and reports the Cα RMSD
between them. That is how close the prediction came to the experiment, for this
protein, on your machine.

Because the prediction is folded from the reference's own sequence, residue *i*
of one is residue *i* of the other and no sequence alignment is needed. The
superposition is Kabsch, by way of Horn's quaternion form. Multi-chain entries
compare the longest chain.

For calibration, SimpleFold-100M at 200 steps against the deposited structures:
crambin (1CRN, 46 aa) came back at **0.36 Å**, and ubiquitin (1UBQ, 76 aa) at
**0.5–0.9 Å** across runs. Ubiquitin's worst residues are the C-terminal tail,
which is flexible in the crystal anyway.

Two colourings come with it. **Pair** gives each structure a flat colour, so you
see the two folds against each other. **Deviation** paints both by how far apart
the counterparts ended up: blue where they agree, red where they do not. That is
where you find the loop the model got wrong.

Either structure can be hidden, from the pills in the viewer or by clicking its
row in the legend, so you can look at one on its own without losing the
alignment or having to fetch it again.

## Folding

Folding happens in a tab of its own, which you open from the **+** in the tab
bar (*Fold a protein*) or with *File › New Fold Tab*. The viewer never folds.
A fold tab takes a sequence, or a PDB ID or UniProt accession whose chains it
folds, and shows the prediction in its own viewport, so you can keep working in
other tabs while it runs.

Folding is a plugin, compiled in by the `simplefold` Cargo feature, which is on
by default. `cargo build --no-default-features` builds the viewer on its own.

SimpleFold is a Python package, so the first fold sets up a private environment
under `~/.profold`: a shallow clone of `apple/ml-simplefold`, a Python 3.11
virtualenv, ESM, and on Apple silicon the MLX backend. It downloads several
gigabytes, once, and progress streams into the left panel. It needs `git` and
[`uv`](https://docs.astral.sh/uv/) on the system:

| | Install uv with |
| --- | --- |
| macOS | `brew install uv`, or `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| Linux | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| Windows | `winget install astral-sh.uv` (and `winget install Git.Git` for git) |

Where it runs depends on the machine. On Apple silicon it runs on the GPU
through MLX. On Linux and Windows with an NVIDIA card it runs on CUDA, since
uv picks the matching PyTorch build for the driver it finds. Anywhere else (an
Intel Mac, an AMD or Intel GPU, no GPU at all) it runs on the CPU, which works
but is slow. Pick the 100M model there.

Model size (100M through 3B), step count and pLDDT are on the fold tab's panel.
100M at 200 steps is the quick setting. Larger models are slower and better.

pLDDT is off by default because it is not free. SimpleFold scores confidence
with the 1.6B checkpoint as a latent encoder, so turning it on downloads ~6.5 GB
the first time, whichever folding model you picked. With it off, SimpleFold
writes a flat B-factor of 100 for every atom, which the app is careful not to
mistake for a real score.

Beyond the folding weights, the first fold also pulls ESM2-3B (~5.3 GB, into
`~/.cache/torch`) for sequence embeddings and a 3.7 GB working cache. All of it
is fetched once. A warm 50-residue fold on the 100M model takes about 25
seconds on Apple silicon. SimpleFold keeps that working cache inside its output
directory, so the app hands it one stable directory rather than a fresh one per
run. Otherwise every fold would re-download the whole 3.7 GB.

Every one of those files is fetched by profold before SimpleFold starts, and
SimpleFold finds them in place. SimpleFold's own downloads write straight to the
final file, so an interrupted one leaves a truncated checkpoint that the next
run believes. profold's are all or nothing: each file is written beside its
destination as `<name>.part`, checked against the size the server declared, and
only then renamed into place. An interrupted download keeps its `.part` and
resumes from there next time, so closing the tab mid-way costs nothing. A file
that was already there is checked against the server once and recorded in
`~/.profold/downloads.json`; one found short is finished rather than fetched
again. The tab says which file it is downloading and how far along it is. The
environment itself counts as installed only once its install has run to the end.

The environment does not go in `~/Library/Application Support`, even though that
is the right place for it. SimpleFold fetches its weights with
`os.system(f"curl -L {url} -o {path}")`, and an unquoted path containing a space
silently downloads the checkpoint to the wrong file. The checkpoint directory is
also passed relative to a working directory the app sets, so a home directory
with a space in it still works.

## Building

```bash
cargo run --release
```

It also takes arguments, each opened in a tab of its own:

```
profold 1UBQ                    # open with that already loaded
profold traj.lammpstrj caffeine # a file, and a molecule found by name
profold --fold                  # an empty fold tab (plugin)
profold --fold 1CRN             # fold an entry's chains (plugin)
profold --compare 1UBQ          # fetch it and fold its sequence alongside (plugin)
```

Without the folding plugin:

```bash
cargo run --release --no-default-features
```

On **macOS** it needs the Xcode Metal toolchain, which GPUI compiles shaders
with:

```bash
xcodebuild -downloadComponent MetalToolchain
```

On **Linux** it needs the development packages GPUI builds against. For
Debian and Ubuntu, [`scripts/linux-deps.sh`](scripts/linux-deps.sh) installs
them (`sudo bash scripts/linux-deps.sh`). On Fedora:

```bash
sudo dnf install gcc clang pkg-config libxkbcommon-x11-devel wayland-devel libxcb-devel fontconfig-devel freetype-devel vulkan-loader alsa-lib-devel openssl-devel
```

On **Windows** it needs the MSVC build tools (Visual Studio's *Desktop
development with C++*), and on Windows on ARM also LLVM (`winget install
LLVM.LLVM`).

`scripts/package.sh <target>` builds the archive a release ships for one
platform, such as `scripts/package.sh aarch64-apple-darwin`.

## Releasing

Releases are built by GitHub Actions
([`.github/workflows/release.yml`](.github/workflows/release.yml)) for six
targets: `aarch64-apple-darwin`, `x86_64-apple-darwin`,
`x86_64-unknown-linux-gnu`, `aarch64-unknown-linux-gnu`,
`x86_64-pc-windows-msvc` and `aarch64-pc-windows-msvc`. To make one, set the
version in `Cargo.toml`, commit, and push a matching tag:

```bash
git tag v0.2.0 && git push origin main v0.2.0
```

The workflow checks that the tag matches `Cargo.toml`, builds and smoke-tests
every target, and publishes the release with the archives, their `.sha256`
files, `install.sh` and `install.ps1`. Installed copies see it the next time
they check. Only these builds update themselves: a build from source shows its
version and says when a release is out, but never overwrites itself.

The gallery pictures in this README are drawn by the renderer itself, ray
traced ones included; the `app-*` ones are screenshots of the window. To
redraw the gallery into `docs/`:

```bash
cargo test --release readme_gallery -- --ignored
```

## How it is put together

| File | |
| --- | --- |
| [`protein.rs`](src/protein.rs) | Structure model: atoms, residues, bonds, cell, frames |
| [`formats/`](src/formats/) | One reader per file format |
| [`elements.rs`](src/elements.rs) | The periodic table: radii, masses, colours |
| [`bonds.rs`](src/bonds.rs) | Bonds from distance, on a grid |
| [`dssp.rs`](src/dssp.rs) | Secondary structure from backbone geometry |
| [`align.rs`](src/align.rs) | Superposition and RMSD |
| [`assembly.rs`](src/assembly.rs) | Biological assemblies from the deposited chains and their symmetry |
| [`select.rs`](src/select.rs) | The selection language |
| [`measure.rs`](src/measure.rs), [`contacts.rs`](src/contacts.rs) | Distances, angles and dihedrals; a ligand's contacts |
| [`analysis.rs`](src/analysis.rs) | Trajectory analysis: RMSD, Rg, RMSF, SASA, contacts, secondary structure |
| [`pbc.rs`](src/pbc.rs) | Periodic boundaries: wrap, make whole, unwrap, replicate |
| [`surface.rs`](src/surface.rs), [`volume.rs`](src/volume.rs) | The molecular surface, and volumes with their isosurfaces |
| [`job.rs`](src/job.rs) | Progress reporting shared by the slow paths |
| [`render.rs`](src/render.rs) | Scene building, camera, and the software rasterizer |
| [`gpu.rs`](src/gpu.rs), [`gpu.wgsl`](src/gpu.wgsl) | The same drawing on the graphics card, with wgpu |
| [`trace.rs`](src/trace.rs), [`gpu_trace.rs`](src/gpu_trace.rs), [`trace.wgsl`](src/trace.wgsl) | Ray tracing, on the CPU and by a compute shader |
| [`movie.rs`](src/movie.rs) | Films: frames, camera moves, overlays and encoders |
| [`settings.rs`](src/settings.rs) | Preferences kept between sessions |
| [`update.rs`](src/update.rs) | Checking GitHub releases, and replacing the app with a new one |
| [`fetch.rs`](src/fetch.rs) | Input classification, RCSB / AlphaFold / UniProt / PubChem |
| [`search.rs`](src/search.rs) | Search by name or gene across the PDB, UniProt and PubChem, and the ranking |
| [`plugins/simplefold/`](src/plugins/simplefold/) | The folding plugin: jobs, settings, and the SimpleFold subprocess |
| [`ui/`](src/ui/) | The GPUI views: the tabbed window, and the viewer each tab is |

GPUI has no 3D pipeline, so profold brings its own. A structure becomes a scene
of sphere and capsule impostors, plus triangles for the cartoon. A camera looks
at it orthographically, and each is resolved against a z-buffer. Two
screen-space passes follow, depth cueing and a depth-discontinuity outline,
which give it the look of a printed molecular illustration. Frames are drawn
with wgpu (Metal, Vulkan, DX12 or GL) on a background thread and uploaded as an
image, cached against the camera so an idle window costs nothing. Translucent
atoms are sorted and blended on top. Without a usable GPU the same frame is
drawn by a software rasterizer, split into bands across every core, and a test
keeps the two images in agreement. The backbone tube is built from capsules
rather than a chain of spheres. Consecutive capsules share an end cap exactly,
so the surface has no seams.

Cartoon needs to know where the helices and sheets are, and most files will not
tell you. A SimpleFold prediction ships no annotation at all. So secondary
structure is assigned from geometry. Amide hydrogens are placed, every
candidate N–H···O=C pair is scored with the Kabsch–Sander energy, and helices
and sheets are read off the resulting n-turn and bridge patterns. It agrees with
the deposited annotations on the structures it was checked against: 1CRN's two
helices and its hairpin, ubiquitin's β-grasp, and porin's barrel at 55% sheet.

The text field is adapted from GPUI's own `input.rs` example, with horizontal
scrolling and change events added.

## License

profold is released under the [MIT License](LICENSE). Use of the app and the
services it talks to is covered by the [Terms of Use](TERMS.md), which include
how it handles your data (no telemetry; files stay on your machine). Credits
for bundled data and dependencies, and the **research-only** license on
SimpleFold's model weights, are in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
