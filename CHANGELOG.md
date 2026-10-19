# Changelog

## Unreleased

### Films
- **Export film** (Image menu, ⇧⌘E) replaces the two fixed GIFs with one
  dialog. A trajectory plays through while the camera holds still, turns,
  rocks, or follows a path you recorded by hand, all at once if you like.
- **Recording.** *Record a camera path* captures what you do in the viewport,
  30 times a second: turning, zooming, panning, and playing or scrubbing the
  trajectory. Stop it, and the film plays it back at full quality.
- **Frames.** Set the first and last frame and keep every Nth. Smoothing
  blends 2×, 4× or 8× between frames (the short way across a periodic cell),
  and ping-pong plays forward then back. Set a frame rate, plus a length in
  seconds or *auto*.
- **Image.** Window size, 720p, 1080p, 4K, square or custom; 1×, 2× or 4×
  antialiasing; a quality setting for video; the background as shown, a
  colour of your own, or transparent.
- **Formats.** GIF and numbered PNGs with transparency, made in the app.
  With ffmpeg installed: MP4 (H.264), and WebM (VP9) and ProRes 4444, both
  with transparency.
- **Overlays.** A frame counter, the simulation step and time (LAMMPS
  timesteps are now read; give the timestep in fs, ps or τ), a scale bar,
  a caption, and the colour legend (hydrophobicity, charge, B-factor,
  sequence, secondary structure, pLDDT, chains, or elements and atom types)
  in a corner of your choice.
- **Preview** plays a small version in the dialog first. While a film is
  being made it shows a progress bar with Cancel, and the settings are
  remembered.
- **Console.** `movie` opens the dialog; `movie out.mp4 frames=1:500:2 fps=30
  turn=360 smooth=3 bg=transparent counter scalebar legend=bottom-left
  caption="…"` makes a
  film without it. `record` and `record stop` capture a path for
  `movie … recorded`.

## 0.3.0 — 2026-10-18

### Looking at a structure
- **Picking.** Hovering the viewport names the atom or residue under the
  pointer; a click picks it out, its residue painted white on the ribbon.
- **Hover card.** The card under the pointer gives a residue's one-letter
  code and full name, its kind, fold, pLDDT or B-factor and position; water,
  ions, common ligands and modified bases go by their names.
- **Sequence.** The whole sequence, chain by chain, in numbered blocks of ten,
  linked both ways with the viewport.
- **Measurements.** Distances, angles and dihedrals between picked atoms,
  drawn as dashes over the molecule with their values.
- **Focus and binding sites.** Fly to a pick; look into a ligand's binding
  site, its residues within 5 Å drawn and named, with its hydrogen bonds,
  salt bridges, metal contacts and clashes as coloured dashes.
- **Selections.** A small selection language after PyMOL's and VMD's
  (`chain`, `resi`, `resn`, `name`, `elem`, `within N of`, `byres`, `and`,
  `or`, `not`…); each selection drawn as sticks or spheres, or hidden.
- **Plots.** A Ramachandran plot, and AlphaFold DB's predicted aligned error
  for AlphaFold models, both linked to the viewport, each folded under a row
  of the facts table.

### Drawing
- **Dark mode**, a ☾/☀ switch in the tab bar, remembered between launches.
- **Colours of your own** for each chain and each kind of ligand, picked
  from its dot in the legend; ligands listed under the chains (HEM ×4).
- **Elements** in their standard colours by default, with a Standard /
  Palette switch; the palette always colours chains and ramps.
- **Surface** style: a Gaussian molecular surface, coloured by any scheme.
- **Nucleic acids** drawn as a tube along the backbone with a rung to each
  base; modified bases found by linkage.
- **Colourings:** B-factor, secondary structure, hydrophobicity (Kyte and
  Doolittle), and charge, which over a surface shows the electrostatic
  potential. Keys for these and for pLDDT in the side panel.
- **Camera:** clipping slab with cut faces drawn solid, perspective, spin
  and rock, saved viewpoints.
- **Water and hydrogens** at a click. Crystal waters are kept now, hidden
  beside a protein by default.

### Structures
- **Biological assemblies** from mmCIF and PDB files.
- **Overlay** any structure on another, paired by sequence alignment, with
  the RMSD and a deviation colouring.

### Files
- **New formats:** BinaryCIF, DCD and XTC. A trajectory opened over a
  structure with as many atoms moves that structure's atoms.
- **Save** the structure, a chain or a binding site as PDB or mmCIF.
- **Sessions** (`.profold`): the structure and how it is shown, in one file.
- **Pictures:** transparent PNGs, and GIFs of a turn or a trajectory.
- **Recent** files and entries in the Samples menu.

### Working
- **A full menu bar**: File, Edit, View, Window and Help, with a shortcut
  for every view and style.
- **Undo and redo** (⌘Z, ⇧⌘Z) of how a structure is shown.
- **Console** (⌘K) of commands, and `profold --run script.txt` to run a
  script of them.
- **Keyboard shortcuts** on one sheet (⌘/).

### Fixes
- Force-field names for histidine and the like (HSD, HIE, CYX…) are read
  as amino acids, not ligands.
- Rainbow colouring runs over the polymer alone.
- Large assemblies build at once, and large DNA cartoons no longer freeze.

## 0.2.1

See the [release](https://github.com/odexion/profold/releases/tag/v0.2.1).
