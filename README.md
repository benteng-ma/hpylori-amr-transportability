# Genomic resistance-rule evaluation in Helicobacter pylori

This repository contains the fixed analysis configuration, executable workflows, regression tests and compact derived results for an external evaluation of clarithromycin and levofloxacin resistance catalogues.

The current revision materials are in [release_updates/v1.0.2](release_updates/v1.0.2/README.md). Version 1.0.2 corrects the population used for secondary assembly-quality comparisons by restricting those comparisons to final-genome-QC-passing records. It provides the revised Tables S2 and S6, Supplementary Figure S1, source data, executable scripts, recorded environment and matching supplementary workbook.

## Reproduce the revision

From the repository root, install the recorded Python dependencies and run the following commands.

~~~sh
python release_updates/v1.0.2/reproduce_qc_callability.py --input results/extended_analysis/qc_callability_samples.csv --output outputs/v1.0.2
python -m pytest -q
~~~

The root scripts/run_extended_analyses.py uses the same QC restriction, and results/extended_analysis/qc_callability_associations.csv contains its corrected 24-row output. All 1,052 sample-drug source records remain in the input table. The revised comparison uses 1,034 QC-passing records. The estimable external comparisons are Ningxia clarithromycin (6 callable, 51 uncallable) and paired-read levofloxacin (42, 4). The primary catalogue definitions, phenotype labels and fixed-rule predictions retain their recorded values.

## Full workflow and public inputs

The declared Python/R/Conda environments and Snakefile reproduce sequence processing and the original analyses. The public inputs comprise NCBI BioProjects PRJNA529500 and PRJNA745492 and Zenodo record 10369064. Their records are listed in metadata/dataset_manifest.csv and metadata/sequence_file_manifest.csv. Source sequences are obtained from those providers and are tracked by their manifests and checksums.

The six main figures and their source data are documented in [release_updates/v1.0.1](release_updates/v1.0.1/README.md). The manuscript directory and earlier update directories retain historical analytical text and numbering. The current R1 quality-comparison table is S6, and its updated Figure S1 and table data are in the v1.0.2 directory.

## Citation and version history

The citation order is Benteng Ma, Bing Chen, Aojian Deng, Ting Cai, Xiao-ming Liu and Fen Wang. Ma and Chen share first authorship. Liu and Wang are co-corresponding authors. Author identifiers and the release date are in CITATION.cff.

- [v1.0.2](https://github.com/benteng-ma/hpylori-amr-transportability/releases/tag/v1.0.2), archived at [10.5281/zenodo.24361066](https://doi.org/10.5281/zenodo.24361066), provides the QC-conditioned R1 correction and matching supplementary data.
- [v1.0.1](https://github.com/benteng-ma/hpylori-amr-transportability/releases/tag/v1.0.1), archived at [10.5281/zenodo.22647135](https://doi.org/10.5281/zenodo.22647135), records the six-figure layout and derived provenance audit.
- [v1.0.0](https://github.com/benteng-ma/hpylori-amr-transportability/releases/tag/v1.0.0), archived at [10.5281/zenodo.22234594](https://doi.org/10.5281/zenodo.22234594), records the original analysis.

The v1.0.2 release tag and its Zenodo archive preserve commit `250e4021afdaf011ff9413d25cd38ae9ada7adcf`. The main branch adds the assigned DOI to citation metadata after archiving.

## Data and licence boundary

Code is MIT licensed. Compact author-generated derived data are CC BY 4.0. Source sequences and third-party phenotype resources retain their original terms. No new sequence data were generated.
