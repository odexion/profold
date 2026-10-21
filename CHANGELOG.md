# Changelog

## Unreleased

### Fixes
- **The drawing styles on Linux and Windows** are Alt+1 to Alt+4, in place
  of 0.5.0's Ctrl+Shift+1 to Ctrl+Shift+4: with Shift held, what a digit
  key sends depends on the keyboard layout, and Alt and a digit does not.
  The Mac keeps ⌃⌘1 to ⌃⌘4.

## 0.5.0 — 2026-10-21

### Rendering
- **Ray traced.** *Image › Ray traced* lights the molecule with rays: soft
  shadows from the light, ambient occlusion that darkens pockets and clefts,
  and edges smoothed by tracing each pixel up to nine times over. The view is
  traced once it holds still (a quarter of a second), with its progress in
  the corner, and any move goes back to the quick frame; saved PNGs are
  traced too, transparent ones included. Every style keeps its own look, ink
  and fog and all, and slabs, perspective and see-through atoms are traced as
  they are drawn. A protein at 1200 × 900 takes about a second. The console
  has `trace on|off`, and `save x.png` follows it.
- **Ray tracing on the graphics card.** Traced pictures, the traced view
  and traced films are traced by a compute shader where there is a card
  that can (Metal, Vulkan, DirectX 12), the CPU otherwise: a retina window's
  view in 0.6 s rather than 4.5, a picture at three times its size in 3.6 s
  rather than 25, a 1080p film frame in 0.4 s rather than 2.2. The CPU
  tracer is faster too, by a tree built for least cost (1.1–1.6×). The
  pictures are the same, to the last bits of a few grazing rays.
  Tracing stays off until asked for, and *Image › Ray traced · slower* says
  what it costs; the film dialog says about how long a traced film will
  take on this machine.
- **Ray-traced films.** The film dialog's *Ray traced* switch, or
  `movie … trace=on`, traces every frame the same way.

### Films
- **A truer preview.** The film dialog's preview plays every frame in
  order, drawn as the film will draw it, at the film's frame rate, and
  widens with the window (296 to 640 pt) so the counter and time over it
  can be read.
- **Caption and legend.** The caption sits top or bottom, left, centre or
  right, in four sizes, and keeps clear of the counter and scale bar. The
  legend takes its own title, size and backdrop, and the writing a colour
  of its own. The console's `movie` takes `caption-at`, `caption-size`,
  `legend-title`, `legend-size`, `legend-box` and `ink`.
- **Colours from the atoms' picker.** The background and the writing are
  chosen in the same picker as the atoms, with its presets and hex and RGB
  fields.

### Labels
- **Labels in 3D.** Pin a residue's name and number, an atom's name or text
  of your own to atoms: *Label* on the picked atom's card, or `label resi
  10-20`, `label atoms name CA`, `label ligand "the drug"` and `label clear`
  in the console. They follow their atoms through a trajectory, and are
  drawn into saved images and films.

### Trajectories
- **Analysis.** *Trajectory analysis* in the side panel measures, frame
  by frame, the RMSD from the first frame (each fitted onto it by Kabsch),
  the radius of gyration, the solvent-accessible surface (Shrake–Rupley)
  and the polar contacts between residues; the RMSF of each atom; and
  helices and strands through time. Each is worked out only when picked,
  in the background, of the atoms chosen: Cα, the backbone, or any atom
  types or elements, and `analyse rmsd of <selection>` for anything else.
  The plots mark their least and greatest values and their span, a line
  marks the frame shown, and a click goes to that frame. *Fit frames onto
  the first* takes out the tumbling. Selections gain `type 1+2`.

### Volumes
- **Isosurfaces.** A Gaussian cube file's grid is read at last, and CCP4,
  MRC and OpenDX maps open too: laid over the structure shown, or on their
  own in their box. The surface is drawn where the values cross a level
  (mean plus 2.5 σ to begin with, − and + to step it), blue, and for an
  orbital's signed values a red lobe where they cross its negative. Solid
  or as a mesh, in the viewport, traced pictures and films. `iso 0.05`,
  `iso mesh`, `iso off` and `iso auto` in the console.

### Simulations
- **Per-atom properties.** A LAMMPS dump keeps every column beyond type and
  position (`c_pe`, `vx`, `q`, stresses…), and an extended XYZ every numeric
  column of `Properties=`, frame by frame. *Property* colours atoms by any of
  them along the palette, from the 1st to the 99th percentile across the
  whole trajectory, with the columns to pick from and the range in the side
  panel and in a film's legend. `color property c_pe -5 -2` sets a range by
  hand.
- **Periodic boundaries.** Under the cell's size: *Wrap* puts every atom
  back in the cell, *Make whole* mends molecules its walls cut (along a
  polymer's chains, and otherwise along bonds), *Unwrap in time* follows a
  trajectory through the walls rather than letting atoms jump back, and
  *Replicate 2×2×2* repeats the cell, frames and properties included.
  *Original* undoes them. In the console: `wrap`, `whole`, `unwrap`,
  `replicate 3 3 1` and `original`.
- **Large dumps in half the memory.** A LAMMPS dump is read 64 MB at a
  time, the next piece while this one is parsed: 0.80 GB rather than
  1.64 GB for 43k atoms × 801 frames, and no slower. RMSD and RMSF no
  longer hold every fitted frame, and SASA and contacts find neighbours in
  a flat cell list. Each gives the same answer as before, which tests check.
- **Cell tools where they mean something.** Wrap, Whole, Unwrap and
  Replicate show for simulation boxes, trajectories, crystals and small
  molecules, not a protein's crystal cell, on one line with short names.
  A lit button is in effect; a click on it turns it off.
- **Selecting by value.** Selections compare: `c_pe < -3`, `vx >= 0.5 and
  z > 10`, `b < 30`, with `x`, `y`, `z` and `b` always there.

### Fixes
- **Panned pictures.** A PNG saved at 2× or 3×, or from the console, kept the
  view's pan at the screen's size, so a view panned off centre came out
  shifted. The pan now grows with the picture, as it already did in films.
- **Ctrl+1 to Ctrl+4 on Linux and Windows** changed the drawing style
  rather than going to tabs 1 to 4: the styles' ⌃⌘ chord is Ctrl alone off
  the Mac. The styles are now Ctrl+Shift+1 to Ctrl+Shift+4 there (still
  ⌃⌘1 to ⌃⌘4 on the Mac).
- **Menus in a narrow window** hang from their own button and stay inside
  the window; the palette used to run off its right edge.

### Viewing
- **⌘-click to focus** (Ctrl-click on Linux and Windows) picks an atom or
  residue and flies the view to it, as the card's Focus button does.

### Menus
- **The structure's name** is set in Charis SIL Regular; the wordmark
  stays bold.
- **The recording light** sits in the top bar, beside the angles, rather
  than over the style pills it used to cover.
- **The bar along the top** is now Open · Samples · Camera · Image · Film ·
  Save. Overlaying a second structure moved into Open; recording a camera
  path and exporting a film have a Film menu of their own; Image is for
  pictures alone. Camera gains *Focus on the pick*.
- **The menu bar** gains Structure (the representations and styles) and
  Camera (centre, focus, reset, spin, rock, perspective, slab), leaving View
  for the theme, console, shortcuts and full screen. File now overlays,
  records a camera path, and saves the structure or the session.

## 0.4.0 — 2026-10-19

### Viewing
- **Roll off the Mac.** Rolling the view in the screen's plane is Ctrl-drag
  on Linux and Windows, as the hints already said; ⌘-drag on the Mac.
- **Honest angle readout.** X and Y are now measured from square-on (looking
  down z, with x across and y up), so the three-quarter opening view reads
  14° and 32° rather than 0°. Typing 0 and 0 gives the square-on view.

### Searching
- **The field says what it holds.** A chip inside the search field names
  what was typed or pasted (PDB ID, UniProt, PubChem, File, a sequence and
  its length, or Unrecognised), beside a ↵ button that acts on it.
- **SAM synthetase (9SV9)** joins the example proteins.

### Look
- **A new logo:** a ribbon wound into a helix that leaves as a strand,
  indigo to orange from N- to C-terminus, with a two-turn version for small
  icons.

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
