# profold

A small desktop app for [Apple's SimpleFold](https://github.com/apple/ml-simplefold),
written in Rust with [GPUI](https://www.gpui.rs).

Paste a protein sequence and it folds it. Paste an identifier and it fetches the
structure. Either way you get a molecule you can turn around.

## Using it

One field takes everything, and the app works out what you gave it:

| Input | What happens |
| --- | --- |
| `1CRN` — a PDB ID | Downloads the experimental structure from RCSB |
| `P69905` — a UniProt accession | Downloads the predicted model from AlphaFold DB |
| `MVLSPADKTNVK…` — residue letters, FASTA and all | Folds it with SimpleFold |

An accession with no AlphaFold model falls back to fetching its sequence from
UniProt and folding that instead.

## Comparing

Next to **Fetch** there is **Compare**, which is the interesting one. It pulls
the deposited structure, folds that structure's *own* sequence with SimpleFold,
superimposes the two, and reports the Cα RMSD between them — how close the
prediction came to the experiment, for this protein, on your machine.

Because the prediction is folded from the reference's own sequence, residue *i*
of one is residue *i* of the other and no sequence alignment is needed. The
superposition is Kabsch, by way of Horn's quaternion form. Multi-chain entries
compare the longest chain.

For calibration, SimpleFold-100M at 200 steps against the deposited structures:
crambin (1CRN, 46 aa) came back at **0.36 Å**, and ubiquitin (1UBQ, 76 aa) at
**0.5–0.9 Å** across runs, its worst residues being the C-terminal tail that is
flexible in the crystal anyway.

Two colourings come with it. **Pair** gives each structure a flat colour, so you
see the two folds against each other. **Deviation** paints both by how far apart
the counterparts ended up — blue where they agree, red where they do not — which
is where you find the loop the model got wrong.

Either structure can be hidden, from the pills in the viewer or by clicking its
row in the legend, so you can look at one on its own without losing the
alignment or having to fetch it again.

In the viewport: drag to orbit, shift-drag to pan, scroll to zoom, double-click
to reframe. The molecule follows the pointer — dragging right brings the near
face right. Predicted models open coloured by pLDDT.

Two styles. **Soft** is the default: smooth shading on paper, distance washing
out. **Illustrative** is the drawn look — ink outlines on both sides of the
silhouette and across every overlap, saturated colour, light quantised into a
few bands, and no depth fade, because a drawn figure does not have one. The line
weight scales with the frame, so it holds on a retina display.

Three representations. **Cartoon** is the default: ribbons for helices, arrows
for strands, cord for everything else. **Tube** is a plain backbone trace, and
**Spacefill** puts every atom at its van der Waals radius. Colour by chain,
by position along the chain, by pLDDT, or by element.

Anything slow — a download, an install, a fold — takes over the viewport with
its own progress. Where the underlying tool reports a percentage the bar is
real; where it does not, it says so by travelling instead.

## Folding

SimpleFold is a Python package, so the first fold sets up a private environment
under `~/.profold`: a shallow clone of `apple/ml-simplefold`, a Python 3.10
virtualenv, and on Apple silicon the MLX backend and ESM. It needs
[`uv`](https://docs.astral.sh/uv/) on the system (`brew install uv`) and it
downloads several gigabytes, once. Progress streams into the left panel.

Model size (100M through 3B), step count and pLDDT are on the sequence panel.
100M at 200 steps is the quick setting; larger models are slower and better.

pLDDT is off by default because it is not free: SimpleFold scores confidence
with the 1.6B checkpoint as a latent encoder, so turning it on downloads ~6.5 GB
the first time, whichever folding model you picked. With it off, SimpleFold
writes a flat B-factor of 100 for every atom, which the app is careful not to
mistake for a real score.

Beyond the folding weights, the first fold also pulls ESM2-3B (~5.3 GB, into
`~/.cache/torch`) for sequence embeddings and a 3.7 GB working cache. All of it
is fetched once; a warm 50-residue fold on the 100M model takes about 25
seconds on Apple silicon. SimpleFold keeps that working cache inside its output
directory, so the app hands it one stable directory rather than a fresh one per
run — otherwise every fold would re-download the whole 3.7 GB.

Not `~/Library/Application Support`, despite that being the right place for it:
SimpleFold fetches its weights with `os.system(f"curl -L {url} -o {path}")`, and
an unquoted path containing a space silently downloads the checkpoint to the
wrong file. The checkpoint directory is also passed relative to a working
directory the app sets, so a home directory with a space in it still works.

## Building

```
cargo run --release
```

It also takes an argument, which saves typing while iterating:

```
profold 1UBQ              # open with that already loaded
profold --compare 1UBQ    # fetch it and fold its sequence alongside
```

Requires the Xcode Metal toolchain, which GPUI compiles shaders with:

```
xcodebuild -downloadComponent MetalToolchain
```

## How it is put together

| File | |
| --- | --- |
| [`protein.rs`](src/protein.rs) | Structure model, PDB and mmCIF parsers |
| [`dssp.rs`](src/dssp.rs) | Secondary structure from backbone geometry |
| [`align.rs`](src/align.rs) | Superposition and RMSD |
| [`job.rs`](src/job.rs) | Progress reporting shared by the slow paths |
| [`render.rs`](src/render.rs) | The software rasterizer |
| [`fetch.rs`](src/fetch.rs) | Input classification, RCSB / AlphaFold / UniProt |
| [`fold.rs`](src/fold.rs) | The SimpleFold subprocess and its environment |
| [`ui/`](src/ui/) | The GPUI views |

GPUI has no 3D pipeline, so the viewer is a software rasterizer: an orthographic
camera, sphere and capsule impostors resolved against a z-buffer, then two
screen-space passes — depth cueing and a depth-discontinuity outline — that give
it the look of a printed molecular illustration. The result is uploaded as an
image each frame, and cached against the camera so an idle window costs nothing.
The backbone tube is built from capsules rather than a chain of spheres:
consecutive capsules share an end cap exactly, so the surface has no seams.

Cartoon needs to know where the helices and sheets are, and most files will not
tell you — a SimpleFold prediction ships no annotation at all. So secondary
structure is assigned from geometry: amide hydrogens are placed, every candidate
N–H···O=C pair is scored with the Kabsch–Sander energy, and helices and sheets
are read off the resulting n-turn and bridge patterns. It agrees with the
deposited annotations on the structures it was checked against — 1CRN's two
helices and its hairpin, ubiquitin's β-grasp, porin's barrel at 55% sheet.

The text field is adapted from GPUI's own `input.rs` example, with horizontal
scrolling and change events added.
