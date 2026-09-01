# External validation of genomic resistance prediction in Helicobacter pylori across independent cohorts

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22234594.svg)](https://doi.org/10.5281/zenodo.22234594)

This public repository contains the frozen configuration, executable analysis workflows, tests, compact derived data, and source data for every manuscript figure in the associated *Microbial Genomics* submission.

## Reproduce the analyses

1. Install the declared Python/R/Conda environments.
2. Obtain the public sequence inputs listed in `metadata/dataset_manifest.csv` and `metadata/sequence_file_manifest.csv` from NCBI BioProjects PRJNA529500 and PRJNA745492 and Zenodo record 10369064.
3. Keep downloaded sequence inputs under the ignored raw-data boundary; this release intentionally does not redistribute large third-party sequence files.
4. Run the frozen workflow (`Snakefile` and analysis scripts) without changing marker catalogues, phenotypes, or external-validation gates.
5. Use `results/source_data/` to reproduce the figures and `results/external_validation/isolate_sequence_phenotype_crosswalk.csv` to audit isolate-level provenance.

## Data and licence boundary

Code is released under the MIT License. Compact author-generated derived data are intended for CC BY 4.0 release. Public source sequences and third-party phenotype resources remain subject to their original repository and publication terms. No new sequence data were generated.

## Repository and archival status

Source repository: https://github.com/benteng-ma/hpylori-amr-transportability

Frozen release: https://github.com/benteng-ma/hpylori-amr-transportability/releases/tag/v1.0.0

Version-specific archive: https://doi.org/10.5281/zenodo.22234594

## Release status

Release `v1.0.0` points to audited commit `531530f2a15ed25585326f4dade85f85975667b9`. A fresh public Windows clone passed 55 public tests with two documented skips for intentionally excluded third-party source files, and all 273 payload-manifest entries matched by path, size, and SHA-256. The archived release is immutable; DOI metadata added to `main` after Zenodo minting does not move the release tag.
