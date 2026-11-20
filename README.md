# Genomic resistance-rule evaluation in Helicobacter pylori

This repository contains the fixed analysis configuration, executable workflows, regression tests and compact derived results for an external evaluation of clarithromycin and levofloxacin resistance catalogues.

The current reader-facing workbook is in [release_updates/v1.0.5](release_updates/v1.0.5/README.md). The QC analysis materials remain in [release_updates/v1.0.2](release_updates/v1.0.2/README.md). Version 1.0.2 corrects the population used for secondary assembly-quality comparisons by restricting those comparisons to final-genome-QC-passing records. It provides the revised Tables S2 and S6, Supplementary Figure S1, source data, executable scripts, recorded environment and matching supplementary workbook.

Version 1.0.3 corrects one source-screening annotation in metadata/dataset_manifest.csv: the study-defined levofloxacin disk-zone categories are resistant <13 mm, intermediate 13–17 mm and susceptible >17 mm. These zone categories are not directly comparable with MIC-based breakpoints. The code, data, numerical results, S2/S6 tables, Figure S1 and supplementary workbook are identical to v1.0.2. For source-zone interpretation, the corrected cutoffs above supersede earlier overlap descriptions in historical planning reports and manuscript/update directories.

Version 1.0.4 streamlines the reader workbook by removing redundant inventory rows, administrative fields and an unpopulated collection-year column. It also clarifies the auxiliary-benchmark and external-cohort labels. The scientific records, statistical calculations and numerical results, Tables S2/S6 and Figure S1 are unchanged. An administrative status row and its generating script now record that the common-method laboratory AST comparison was not performed. The v1.0.4 workbook is retained in its original release directory.

Version 1.0.5 clarifies the workbook README description of the 526 linked isolate records and documents the existing input checksums, CheckM2 database provenance and paired-read calling rules. For current data use the v1.0.5 workbook; reproduction scripts, QC source files, Tables S2/S6 and Figure S1 retain their v1.0.2 paths. Scientific data, classification logic and numerical results are unchanged.

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

- [v1.0.5](https://github.com/benteng-ma/hpylori-amr-transportability/releases/tag/v1.0.5) provides the matching reader workbook and reproducibility clarifications.
- [v1.0.4](https://github.com/benteng-ma/hpylori-amr-transportability/releases/tag/v1.0.4), archived at [10.5281/zenodo.24384297](https://doi.org/10.5281/zenodo.24384297), provides the streamlined reader workbook without changing the analyses.
- [v1.0.3](https://github.com/benteng-ma/hpylori-amr-transportability/releases/tag/v1.0.3), archived at [10.5281/zenodo.24361451](https://doi.org/10.5281/zenodo.24361451), corrects the dataset-manifest levofloxacin zone annotation; the scientific payload is unchanged from v1.0.2.
- [v1.0.2](https://github.com/benteng-ma/hpylori-amr-transportability/releases/tag/v1.0.2), archived at [10.5281/zenodo.24361066](https://doi.org/10.5281/zenodo.24361066), provides the QC-conditioned R1 correction and matching supplementary data.
- [v1.0.1](https://github.com/benteng-ma/hpylori-amr-transportability/releases/tag/v1.0.1), archived at [10.5281/zenodo.22647135](https://doi.org/10.5281/zenodo.22647135), records the six-figure layout and derived provenance audit.
- [v1.0.0](https://github.com/benteng-ma/hpylori-amr-transportability/releases/tag/v1.0.0), archived at [10.5281/zenodo.22234594](https://doi.org/10.5281/zenodo.22234594), records the original analysis.

The v1.0.2 release tag and its Zenodo archive preserve commit `250e4021afdaf011ff9413d25cd38ae9ada7adcf`.

The v1.0.3 release tag and its Zenodo archive preserve commit `8b90c137e4904fa8a678221a7f00eea9fcb7f370`.

The v1.0.4 release tag and its Zenodo archive preserve commit `60498854bc9259c85113db6935ff20ebc9fa818d`.

## Data and licence boundary

Code is MIT licensed. Compact author-generated derived data are CC BY 4.0. Source sequences and third-party phenotype resources retain their original terms. No new sequence data were generated.
