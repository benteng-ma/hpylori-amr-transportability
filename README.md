# Genomic resistance-rule evaluation in Helicobacter pylori

This repository contains the fixed analysis configuration, executable workflows, regression tests and compact derived results for an external evaluation of clarithromycin and levofloxacin resistance catalogues.

The updated six-figure layout, derived provenance audit and current figure/table numbering are documented in [release_updates/v1.0.1](release_updates/v1.0.1/README.md). The original primary predictions are unchanged. The `manuscript/` directory retains the original analytical text; the update directory contains the revised scientific text.

## Reproduce the analyses

1. Install the declared Python/R/Conda environments.
2. Obtain the public sequence inputs listed in `metadata/dataset_manifest.csv` and `metadata/sequence_file_manifest.csv`: NCBI BioProjects PRJNA529500 and PRJNA745492, and Zenodo record 10369064.
3. Keep downloaded inputs under the ignored raw-data directories. This repository does not redistribute large third-party sequences.
4. Run the frozen workflow (`Snakefile` and analysis scripts) without changing marker catalogues, phenotypes or evaluation criteria.
5. Follow the update-directory instructions for the corrected Figure 1 and Figure 5 layouts.

## Original archive

The current citation order is Benteng Ma, Bing Chen, Aojian Deng, Ting Cai, Xiao-ming Liu and Fen Wang. Ma and Chen share first authorship; Liu and Wang are co-corresponding authors. See CITATION.cff for machine-readable author identifiers.

The immutable analytical release is [v1.0.0](https://github.com/benteng-ma/hpylori-amr-transportability/releases/tag/v1.0.0), archived at [10.5281/zenodo.22234594](https://doi.org/10.5281/zenodo.22234594). That DOI identifies the original files, not subsequent presentation or provenance additions. Do not move the old tag.

## Data and licence boundary

Code is MIT licensed. Compact author-generated derived data are CC BY 4.0. Source sequences and third-party phenotype resources retain their original terms. No new sequence data were generated.
