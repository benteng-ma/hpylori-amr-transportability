# Presentation and provenance update

This directory supplements the original analytical release. It does not add a cohort, retrain a predictor, change a phenotype, or revise a fixed-rule prediction.

## Contents

- Figures 1-6 correspond to the current six-figure manuscript layout. Figure 1 preserves the genome-QC guides but removes the unrelated 90% target-callability guide; Figure 5 separates the two external cohorts.
- `S12_provenance_audit.csv` is the 57-isolate derived audit formerly supplied as Supplementary Data S21. One source-recall discrepancy is attributable to N87Y rule scope; 25 remain unresolved. Independent callers agree on all 57 deposited assemblies. This does not verify the original specimen linkage.
- `numbering_map.json` maps legacy labels to the current Figures S1-S12 and Tables S1-S20. Legacy Table S12 is now an unnumbered technical manifest. Source-file identifiers in the original analysis remain unchanged.
- `manuscript_scientific_text.md` provides the corresponding scientific text and current cross-references. The root `manuscript/` directory remains the original analytical source, not the current submission layout.

## Reproduce the corrected figures

From the repository root, using the declared Python environment:

```sh
python release_updates/v1.0.1/render_figure1.py
python release_updates/v1.0.1/render_figure5.py --repository-root .
python -m pytest -q
```

The renderers use unchanged compact source tables and do not fit models. Outputs go into dedicated rendered_figure1 and rendered_figure5 directories. The complete original sequence-processing workflow remains available in the root repository and requires the public inputs listed in its manifests. A compact-output regression test is not a new independent validation or a full fresh genome-processing run.

Code is MIT licensed. Author-generated compact derived data are CC BY 4.0; third-party inputs retain their original terms.

The Snippy runtime safety test is explicitly Linux-only. Windows regression checks skip that test rather than resolving `/tmp` as a drive-relative Windows path. The Linux runtime helper and its path-safety checks are unchanged.

The descriptive `canonical_root` field is relative to the checkout, not an author-specific drive. Scripts continue to resolve their own repository paths. This portability correction changes no scientific setting or result.
