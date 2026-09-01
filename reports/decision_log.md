# Detailed decision log

See the repository-level `DECISION_LOG.md` for the signed Phase 0 decisions. The machine-readable gate is `reports/phase0_decision.json` and the frozen downstream plan is `reports/analysis_plan.md`.

## 2026-09-01 — public-repository sanitization

- Excluded local-path-bearing per-sample fastp HTML/JSON reports and the resolved Conda prefix file from the public repository while retaining compact QC summaries and exact environment lock files.
- Added fail-closed checks for Windows and Unix local paths, excluded cloud-sync sidecars, and shipped a protective `.gitignore`. The repository must be created privately first and reviewed before public release.
- This changes only the public packaging boundary. No analytical data, result, figure, table, authorship, or conclusion changed.

## 2026-09-01 — natural-language article title

- Replaced the coined `Audit-first` title with `External validation of genomic resistance prediction in Helicobacter pylori across independent cohorts` on all journal-facing and release-facing surfaces.
- The new wording has no hyphen, dash, colon, or branded framework label and makes no new scientific claim. All frozen analyses, results, authorship, figures, tables, and evidence boundaries remain unchanged.

## 2026-09-01 — supplementary presentation and open-repository boundary

- Supplementary Figures S1-S8 use panel letters, axes, legends, thresholds, and statistical annotations but no narrative panel or figure titles; captions carry the narrative explanation.
- A future GitHub/Zenodo release is a reproducibility layer, not a substitute for reviewer-facing supplementary methods, figures, and scientific tables. Only the file-path inventory is a candidate for later compression after stable repository URLs and a version DOI exist.
