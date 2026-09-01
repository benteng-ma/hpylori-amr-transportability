# Analysis plan

Status: `FROZEN_2026-08-30_CONDITIONAL_GO_CLR_LVX_BINARY`

This plan was frozen before full-cohort sequence processing or performance inspection. Any change requires an entry in `DECISION_LOG.md` and `reports/deviations.md`; result-driven threshold changes are prohibited.

## Scope and estimands

Primary drugs are clarithromycin and levofloxacin. The baseline rules are the frozen HpGP 23S and gyrA panels in `config/panels.yaml`. The primary estimand is out-of-study diagnostic performance on phenotype-linked isolates, with false-susceptible rate as the safety-priority error. Marker discovery is outside the primary analysis.

HpGP (414 exactly linked isolates) is the catalogue/global internal cohort. Ningxia 2022 (60 exact MIC-linked isolates) is the primary external cohort. Zenodo 10369064 (52 exact categorical/read-linked isolates) is a secondary phenotype-method stress test. ML296 is overlap-audit only. US Mayo, Linqu, Israel, Swiss, Eastern China, and China landscape datasets cannot contribute performance observations until an exact public phenotype crosswalk exists.

## Locked phenotype rules

- Preserve every study's original label and breakpoint provenance.
- Primary binary analyses exclude I/ambiguous observations; I-as-S and I-as-R are separate sensitivity analyses.
- Ningxia clarithromycin uses the publication P-R label in the primary analysis; the six MIC=0.5 mg/L observations form a prespecified borderline sensitivity set.
- Zenodo original categorical labels are used only in the secondary stress test. No threshold is reconstructed from model performance. The overlapping 13/17 mm levofloxacin endpoints are reported as a source limitation.
- MIC censoring operators are retained. MIC-gradient analyses use censored/ordinal methods rather than replacing `>` or `<` by an exact value.
- Continuous MIC-gradient analysis is cohort-limited (principally Ningxia) and cannot be presented as a complete multi-cohort endpoint because HpGP public phenotypes are binary and Zenodo public phenotypes are categorical.

## Sequence and marker calling

- Reference coordinate system: 26695 / NC_000915.1. Common 23S labels A2142/A2143 map to U27270.1 feature positions 2143/2144.
- Raw reads are the preferred source for 23S because assemblies may collapse or split the two copies. A locus is callable only with depth >=10, base quality >=20, mapping quality >=20, and major-allele fraction >=0.80. Allele fractions 0.20-0.80 are reported as mixed/heteroresistant and never forced to wild type.
- Assembly calls require a target-spanning alignment. No hit or a broken repeat is `uncallable`, not susceptible.
- gyrA residues 87, 88, and 91 are translated relative to the 26695 CDS. The published baseline is applied without marker additions.

## Independence, lineage, and splits

Before modeling, exact isolate/BioSample/run/assembly overlaps are removed and all genomes receive Mash/skani screening. A near-clone group is connected by ANI >=99.9%, both aligned fractions >=0.90, and Mash distance <=0.001; the group, not isolate, is the resampling unit.

Validation priority is leave-study-out, then leave-country-out within HpGP where feasible, then leave-lineage-out, clone-grouped split, and finally random-isolate split as an optimism comparator only. A held-out lineage result is estimable only when it contains at least 30 isolates and at least 10 S and 10 R outcomes for the drug; otherwise it is `INSUFFICIENT_DATA`.

## Metrics and uncertainty

Report sensitivity, specificity, false-susceptible rate, false-resistant rate, balanced accuracy, MCC, PPV, and NPV with exact/binomial intervals. Use cluster bootstrap by patient when available, otherwise near-clone group, with 2,000 replicates and seed 20260830. Report cohort-specific results before any pooled estimate. Probability-model calibration is reported only for prespecified probabilistic models; deterministic panels receive empirical resistance risk by marker group and MIC-gradient analysis, not pseudo-calibration.

## Transportability labels

The labels in `config/validation.yaml` are assigned using frozen rules. A panel is `ROBUSTLY_TRANSPORTABLE` only if every eligible external cohort has estimable sensitivity and specificity >=0.90, false-susceptible rate <=0.10, no cohort differs by >0.10 in either sensitivity or specificity, and no eligible lineage violates those gates. Failure attributable to phenotype method, study, or lineage is labeled accordingly; sparse gates yield `INSUFFICIENT_DATA` rather than success.

## Negative controls and robustness

- Compare random-isolate, near-clone-grouped, leave-country-out, leave-lineage-out, and HpGP-to-external validation under the same frozen-feature logistic specification.
- Repeat primary deterministic-panel estimates after retaining the lexicographically first isolate per near-clone component; patient-level thinning is `NOT_ESTIMABLE` because public patient identifiers are unavailable.
- Keep Zenodo raw-read predictions authoritative and quantify their concordance/callability against support-only assemblies without counting either sample twice.
- Report original labels, I-excluded, I-as-S, I-as-R, prespecified borderline-MIC exclusion, cohort/AST-method strata, and high-quality-genome attrition.
- Run 100 training-label permutations for mutation-only models; only HpGP training labels are permuted and every test label remains untouched.
- Construct 100 size- and prevalence-matched random core-SNP panels. Candidate SNP features and frequency matching are defined from HpGP sequences without phenotype use; panels are trained in HpGP and evaluated under the same split/external framework. External sequences do not alter SNP feature selection.
- Random seeds are 20260830. Negative controls cannot promote markers, alter the frozen panel, or replace the primary estimates.

## Locked until conditions are met

No final Results prose, exploratory marker promotion, performance-driven breakpoint choice, or strong non-Chinese geographic claim is allowed until complete sequence acquisition, full duplicate screening, lineage assignment, and the prespecified QC attrition table are committed.
