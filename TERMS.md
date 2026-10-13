# Terms of Use

These terms cover your use of profold, the molecular viewer in this
repository. The source code itself is licensed under the [MIT License](LICENSE).
Where these terms and the license differ, the license governs the code.

## 1. What profold is

profold is free, open-source software for looking at molecular structures and,
with the optional SimpleFold plugin, for predicting protein structures on your
own computer. There is no profold account, server or paid service.

## 2. No warranty

profold is provided **as is**, without warranty of any kind, as the
[MIT License](LICENSE) sets out. The authors are not liable for any loss or
damage arising from its use.

## 3. Scientific and medical use

profold is a research and education tool.

- **Predicted structures are models, not measurements.** A SimpleFold
  prediction, an AlphaFold DB model, and any RMSD or pLDDT figure profold
  reports can be wrong. Check results against experimental data before relying
  on them.
- **profold is not a medical device.** Do not use it, or anything it produces,
  to diagnose, treat or make clinical decisions about any person.
- Bonds, secondary structure and element assignments that profold works out
  from geometry are approximations and may not match the deposited or true
  chemistry.

## 4. Third-party services

To fetch structures or set up the folding plugin, profold contacts services
that profold's authors do not run. Your use of each is governed by that
provider's own terms:

| Service | Used for | Terms |
| --- | --- | --- |
| RCSB Protein Data Bank (`files.rcsb.org`) | Downloading PDB entries | [RCSB PDB policies](https://www.rcsb.org/pages/policies) |
| AlphaFold Protein Structure Database (`alphafold.ebi.ac.uk`) | Downloading predicted models | [AlphaFold DB terms](https://alphafold.ebi.ac.uk/faq), [EMBL-EBI terms of use](https://www.ebi.ac.uk/about/terms-of-use) |
| UniProt (`rest.uniprot.org`) | Downloading sequences | [UniProt license and disclaimer](https://www.uniprot.org/help/license) |
| GitHub (`github.com`) | Cloning SimpleFold and ESM (plugin) | [GitHub Terms of Service](https://docs.github.com/site-policy/github-terms/github-terms-of-service) |
| Python Package Index, through `uv` | Installing Python dependencies (plugin) | [PyPI terms](https://policies.python.org/pypi.org/Terms-of-Service/) |
| Apple and Meta model hosting | Downloading SimpleFold and ESM2 weights (plugin) | The licenses shipped with each model |

Please be considerate of these free public services: don't use profold to
make automated bulk downloads.

## 5. SimpleFold model weights are for research only

The SimpleFold plugin installs Apple's
[`ml-simplefold`](https://github.com/apple/ml-simplefold). Its **code** is
MIT-licensed, but its **model weights** come under Apple's separate
`LICENSE_MODEL`, which limits use to **non-commercial scientific research**.
profold does not include or redistribute those weights. They are downloaded to
your computer only when you first fold a sequence. By using the plugin you are
responsible for complying with that license. If your use is commercial, build
profold without the plugin:

```bash
cargo build --release --no-default-features
```

## 6. Your data and privacy

- profold has **no telemetry, analytics, accounts or crash reporting**.
- Files you open are read on your computer and never uploaded.
- Network requests happen only when you ask for something that needs them:
  - fetching a PDB ID or UniProt accession sends that identifier to the
    services above;
  - setting up the folding plugin downloads its code and weights.

  Like any web request, these reveal your IP address to the provider.
- Sequences you fold are processed **locally** by SimpleFold. They are not sent
  anywhere by profold.
- profold stores data only on your computer, under `~/.profold`:
  - `settings.json` holds your viewer preferences;
  - the plugin keeps its environment, model checkpoints and predictions there;
  - ESM2 weights go in `~/.cache/torch`.

  Delete these to remove everything profold has saved.

## 7. Content you open or produce

Structures you open belong to whoever holds rights in them. Images you save
with profold are yours to use, subject to the rights in the underlying
structure data. For example, AlphaFold DB and UniProt data are licensed under
CC BY 4.0 and call for attribution. See
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## 8. Changes

These terms may be updated in this repository. The version in the commit you
are using is the one that applies to it.

## 9. Contact

Questions and reports go to the
[issue tracker](https://github.com/odexion/profold/issues).
