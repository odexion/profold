# Third-Party Notices

profold's own code is under the [MIT License](LICENSE). It builds on, and works
with, the third-party software and data listed here, each under its own terms.

## Software compiled into profold

| Component | License | Notes |
| --- | --- | --- |
| [GPUI](https://www.gpui.rs) | Apache-2.0 | UI framework. [`src/ui/input.rs`](src/ui/input.rs) is adapted from GPUI's `input.rs` example and keeps the Apache-2.0 terms for the parts taken from it. |
| [wgpu](https://wgpu.rs) | MIT OR Apache-2.0 | GPU rendering |
| [rayon](https://github.com/rayon-rs/rayon), [image](https://github.com/image-rs/image), [ureq](https://github.com/algesten/ureq), [flate2](https://github.com/rust-lang/flate2-rs), [serde_json](https://github.com/serde-rs/json) and others | MIT OR Apache-2.0 | |

profold has nearly 800 crates in its dependency tree, across all platforms. Nearly all of them are
MIT, Apache-2.0, or available under either. The rest are under other permissive
licenses: BSD-2/3-Clause, ISC, Zlib, Unicode-3.0, CC0-1.0, CDLA-Permissive-2.0,
and MPL-2.0 for a few platform and build-time crates. Their license texts ship
with their sources. To list every crate with its license, for example when
you distribute a binary, run:

```bash
cargo install cargo-license && cargo license
```

## Software installed by the SimpleFold plugin

None of this is part of profold or its repository. It is downloaded at run time
into `~/.profold` and `~/.cache/torch`, only when the plugin is used.

| Component | License |
| --- | --- |
| [Apple ml-simplefold](https://github.com/apple/ml-simplefold) code | MIT, Copyright © 2025 Apple Inc. |
| SimpleFold model weights | Apple Machine Learning Research Model license (`LICENSE_MODEL`): **non-commercial research use only** |
| [ESM / ESM2](https://github.com/facebookresearch/esm) code and weights | MIT, Meta Platforms, Inc. |
| [MLX](https://github.com/ml-explore/mlx) | MIT, Apple Inc. |
| Other Python dependencies | Their own licenses, installed with `uv` from PyPI |

## Sample data

The structures bundled in [`samples/`](samples/README.md):

| File | Source | Terms |
| --- | --- | --- |
| `1crn.pdb.gz`, `1ubq.pdb.gz`, `1ema.pdb.gz`, `4hhb.pdb.gz` | [RCSB PDB](https://www.rcsb.org), wwPDB archive | [CC0 1.0](https://www.wwpdb.org/about/usage-policies): free to use, with citation of the entries' primary publications appreciated |
| `caffeine.sdf` | [PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/2519) CID 2519 | NCBI/PubChem data, freely available |
| `c60.xyz`, `rocksalt.cif`, `silicon.vasp`, `lj-melt.lammpstrj.gz` | Generated for profold | MIT, with the rest of the repository |

Structure references for the bundled PDB entries:

- **1CRN**: Teeter, M.M. (1984) *PNAS* 81, 6014–6018.
- **1UBQ**: Vijay-Kumar, S., Bugg, C.E., Cook, W.J. (1987) *J. Mol. Biol.* 194, 531–544.
- **1EMA**: Ormö, M. et al. (1996) *Science* 273, 1392–1395.
- **4HHB**: Fermi, G., Perutz, M.F., Shaanan, B., Fourme, R. (1984) *J. Mol. Biol.* 175, 159–174.

## Data fetched at run time

| Source | License |
| --- | --- |
| RCSB PDB / wwPDB | CC0 1.0 |
| AlphaFold Protein Structure Database | CC BY 4.0: cite Jumper et al. (2021) *Nature* 596, 583–589, and Varadi et al. (2024) *Nucleic Acids Res.* 52, D368–D375 |
| UniProt | CC BY 4.0 |

If you publish an image of a structure you fetched through profold, please
credit the source database and the structure's authors.
