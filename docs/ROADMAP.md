# Roadmap

What OVITO and VMD offer that profold does not yet, ranked by how much each
would be used against what it costs. Rendering quality (ray-traced pictures
with shadows and ambient occlusion) was the first of these and is done.

| # | Feature | Why here | Size |
| --- | --- | --- | --- |
| 1 | **Per-atom properties.** Keep every column of a LAMMPS dump or an extended XYZ (`c_pe`, `vx`, `q`, stress…) and colour, size and select by any of them, with a ramp and a range. | The gap every simulation user hits first, and what most of the rest builds on. The readers already find the columns; they skip them. | S–M |
| 2 | **Periodic boundaries.** Wrap atoms into the cell, unwrap molecules split across it, unwrap a trajectory through time, and replicate the cell (n×m×k images). | Expected of any viewer of simulation boxes; the cell is drawn but nothing uses it. | S–M |
| 3 | **Volumetric data and isosurfaces.** Read the grid of Gaussian cube (already parsed for its box), CCP4/MRC maps and OpenDX; draw isosurfaces at a level picked on a histogram, two-signed for orbitals. | Orbitals, densities, cryo-EM and electrostatic maps: the largest gap for chemists and structural biologists. `surface.rs` already meshes a field. | M |
| 4 | **Trajectory analysis with plots.** RMSD and RMSF against a reference, fit every frame onto the first, radius of gyration, SASA, hydrogen bonds and salt bridges over time, secondary structure along time (VMD's Timeline). One plot panel, clicking a point goes to its frame. | Most of what VMD is opened for after looking. Kabsch, DSSP, contacts and the plot drawing exist. | M |
| 5 | **Labels in 3D.** Residue and atom names, measurements and free text pinned to atoms, drawn over the frame and into images and films. | Small, and every figure wants it. | S |
| 6 | **Expression selection and a modifier list.** `Position.Z > 10 && c_pe < -3` style selection over properties, and a stack of steps (select, slice, delete, colour code, replicate) re-run on every frame, as in OVITO's pipeline. | Turns 1 and 2 into a workflow; the selection language is half of it. | M–L |
| 7 | **Materials analysis.** RDF, MSD, coordination number, cluster analysis, displacement from a reference frame, atomic strain, spatial binning, histograms and time averages, each exportable as a table. | OVITO's daily tools; needs 1 and the plot panel from 4. | M |
| 8 | **More topology and trajectory formats.** PSF, AMBER prmtop and NetCDF, GROMACS TPR/TOP and TRR, H5MD, binary LAMMPS dumps; writing trajectories out. | Opens MD outputs as people have them; each reader is small, there are many. | M |
| 9 | **More representations.** Lines (fast for huge systems), hydrogen bonds drawn as dashes, bonds recomputed each frame, coordination polyhedra for crystals. | Polyhedra and per-frame bonds matter for materials; lines for scale. | S–M |
| 10 | **Structure identification.** Common neighbour analysis (and adaptive), polyhedral template matching, Voronoi, Wigner–Seitz defects. | Core for metals and crystals, but a narrower audience. Dislocation analysis (DXA) and grain segmentation after. | L |
| 11 | **Vector glyphs.** Arrows for velocities, forces, dipoles, displacements. | Cheap once 1 exists. | S |
| 12 | **Scripting.** A real language behind the console (Python or Rhai): variables, loops, per-frame callbacks, and a headless batch mode. | Power users and reproducible figures; the console commands are the API to wrap. | M–L |
| 13 | **Meshes from atoms.** Alpha-shape surfaces around material (OVITO's construct surface mesh), and pore and cavity surfaces. | Specialised. | M |

Out of scope for a viewer: building systems (psfgen, solvate, ionize,
mutation), running or steering simulations (IMD, QwikMD), stereo and VR, and
multiple sequence alignment.
