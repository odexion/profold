# Roadmap

What OVITO and VMD offer that profold does not yet, ranked by how much each
would be used against what it costs. Rendering quality (ray-traced pictures
and films with shadows and ambient occlusion) came first, then 1–5 below,
each in a first, deliberately small form: what is done is struck through,
with what was left for later.

| # | Feature | Why here | Size |
| --- | --- | --- | --- |
| 1 | ~~**Per-atom properties.**~~ | Done: colour by and select on any column. Left: size by property. | done |
| 2 | ~~**Periodic boundaries.**~~ | Done: wrap, make whole, unwrap in time, replicate. Left: a replicate size picker in the panel (the console takes any). | done |
| 3 | ~~**Volumetric data and isosurfaces.**~~ | Done: cube, CCP4/MRC, OpenDX; solid or mesh, two lobes. Left: a histogram to pick the level on, see-through surfaces. | done |
| 4 | ~~**Trajectory analysis with plots.**~~ | Done: RMSD, Rg, SASA, polar contacts, RMSF, secondary structure, fit frames. Left: true hydrogen bonds and salt bridges, a chosen reference and selection, exporting the series. | done |
| 5 | ~~**Labels in 3D.**~~ | Done: residue, atom and free-text labels in the view, images and films. Left: measurements in images and films, saving labels with a session. | done |
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
