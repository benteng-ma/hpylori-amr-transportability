# Data audit

Audit cutoff: 2026-08-30.

## Outcome

Three datasets satisfy the minimum public sequence-to-phenotype linkage requirement for clarithromycin and levofloxacin:

| Dataset | Public phenotype | Public sequence | Exact links | Phenotype grade | Phase 1 role |
|---|---:|---:|---:|---|---|
| HpGP global | 419 binary | 1,017 current assembly records | 414 | A, exact public binary | frozen catalogue source and global internal cohort |
| Ningxia 2022 | 60 MIC | 60 assemblies | 60 | A, exact MIC | primary external cohort |
| Zenodo 10369064 | 52 categorical | 52 paired-read files | 52 | B, custom disk thresholds | secondary external stress test |

The five unmatched HpGP phenotype IDs are ISR-003 and ISR-007 through ISR-010. They remain excluded unless an exact accession crosswalk is recovered; no fuzzy linkage is permitted.

## Phenotype audit

For Ningxia, Tables S1 and S9 yielded 60 MIC rows and 60 published phenotype-resistance rows. Recalculation agreed exactly for levofloxacin, amoxicillin, furazolidone, tetracycline, and metronidazole after confirming that the furazolidone resistance rule is inclusive at 4 mg/L. Six clarithromycin isolates at exactly 0.5 mg/L are S in the publication but I under the audited three-class interpretation. Both fields are retained; neither is overwritten.

For Zenodo 10369064, the detailed methods and Supplementary Figure 1 establish Kirby-Bauer disk diffusion, despite an abstract-level description elsewhere. The levofloxacin figure assigns 13 mm to both R and I and 17 mm to both I and S. Original categorical labels are therefore preserved; this cohort cannot define or optimize MIC breakpoints and is secondary only.

US Mayo 2025, Linqu 2025, Swiss 2019, Eastern China 2025, and China landscape 2026 contain useful sequence or aggregate phenotype information but no public isolate-level sequence-to-AST crosswalk. Israel PRJEB37854 returned zero ENA/NCBI sequence records on the audit date. None is eligible for performance estimation in the current public-only benchmark.

## Sequence and lineage smoke tests

- Nine assemblies (three HpGP, three Ningxia, three Linqu) passed basic size/contiguity audit.
- All nine yielded full-length gyrA target calls. Seven yielded callable 23S target loci; SHZY01 and SHZY02 were uncallable at 23S because the repeated rRNA region is fragmented in the public assemblies. They are missing, not wild type.
- Thirty-six pairwise Mash/skani comparisons found no representative near-clone pair using the frozen smoke threshold: ANI >=99.9%, both aligned fractions >=0.90, and Mash distance <=0.001.
- Zenodo isolate S27 contained 137,016 reads; 59.11% of primary reads mapped to 26695. All seven audited reference marker loci had nonzero depth, but depth was only 4-14, reinforcing a prespecified minimum-depth rule for Phase 1.

Machine-readable evidence is in `metadata/phase0/`, `metadata/sequence_file_manifest.csv`, and `results/smoke/`.
