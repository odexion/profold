# Changelog

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
