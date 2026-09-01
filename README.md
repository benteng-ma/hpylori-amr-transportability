# External validation of genomic resistance prediction in Helicobacter pylori across independent cohorts

This versioned release candidate contains the frozen configuration, executable analysis workflows, tests, compact derived data, and source data for every manuscript figure in the associated *Microbial Genomics* submission.

## Reproduce the analyses

1. Install the declared Python/R/Conda environments.
2. Obtain the public sequence inputs listed in `metadata/dataset_manifest.csv` and `metadata/sequence_file_manifest.csv` from NCBI BioProjects PRJNA529500 and PRJNA745492 and Zenodo record 10369064.
3. Keep downloaded sequence inputs under the ignored raw-data boundary; this release intentionally does not redistribute large third-party sequence files.
4. Run the frozen workflow (`Snakefile` and analysis scripts) without changing marker catalogues, phenotypes, or external-validation gates.
5. Use `results/source_data/` to reproduce the figures and `results/external_validation/isolate_sequence_phenotype_crosswalk.csv` to audit isolate-level provenance.

## Data and licence boundary

Code is released under the MIT License. Compact author-generated derived data are intended for CC BY 4.0 release. Public source sequences and third-party phenotype resources remain subject to their original repository and publication terms. No new sequence data were generated.

## Pre-release status

This is an author-review candidate. It may first be published as a private GitHub repository to establish the final repository URL. Insert the final repository URL and version DOI in `CITATION.cff` and the manuscript, complete the release checklist, and verify a clean clone before making the repository public or submitting the article.
