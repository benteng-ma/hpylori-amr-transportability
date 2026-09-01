## Supplementary tables



The tables below provide reviewer-readable results within this Word/PDF supplement. Complete row-level and replicate-level data remain in Supplementary Data S1-S20.xlsx; the workbook is the machine-readable companion to, not a replacement for, the tables printed here. Worksheet S12 is a provenance manifest for figure source-data files and is retained only in the machine-readable workbook and repository; because it contains file inventory metadata rather than scientific results, it is not reproduced as a paginated table. NE, not estimable; CLR, clarithromycin; LVX, levofloxacin; FSR, false-susceptible rate; BA, balanced accuracy; QC, quality control.

### Table S1. Isolate-sequence-phenotype-prediction crosswalk: cohort summary

| Cohort | Geography | Materialized | Final-QC pass | CLR phenotype | LVX phenotype | Collection years |
| --- | --- | --- | --- | --- | --- | --- |
| HpGP | Multiple (39) | 414 | 414 | 414 | 414 | NE |
| Ningxia | China | 60 | 57 | 60 | 60 | NE |
| Read cohort | China | 52 | 46 | 52 | 52 | NE |

*Note:* The complete 526-row isolate-level crosswalk, accessions, checksums, QC fields, phenotypes, marker calls and predictions are provided in workbook sheet S01_Crosswalk.

### Table S2. Complete genome quality and target callability

| Cohort | Drug | Phenotype-linked | Callable | Uncallable | Callability |
| --- | --- | --- | --- | --- | --- |
| Ningxia | Clarithromycin | 60 | 6 | 54 | 10.0% |
| Ningxia | Levofloxacin | 60 | 57 | 3 | 95.0% |
| HpGP | Clarithromycin | 414 | 414 | 0 | 100.0% |
| HpGP | Levofloxacin | 414 | 414 | 0 | 100.0% |
| Read cohort | Clarithromycin | 52 | 46 | 6 | 88.5% |
| Read cohort | Levofloxacin | 52 | 42 | 10 | 80.8% |

### Table S3a. Phenotype-blind core-SNP lineage composition

| Cohort | Fixed SNP cluster | n | Within-cohort proportion |
| --- | --- | --- | --- |
| Ningxia | SNP_CLUSTER_01 | 4 | 7.0% |
| Ningxia | SNP_CLUSTER_03 | 52 | 91.2% |
| Ningxia | SNP_CLUSTER_04 | 1 | 1.8% |
| HpGP | SNP_CLUSTER_01 | 122 | 29.5% |
| HpGP | SNP_CLUSTER_02 | 2 | 0.5% |
| HpGP | SNP_CLUSTER_03 | 60 | 14.5% |
| HpGP | SNP_CLUSTER_04 | 109 | 26.3% |
| HpGP | SNP_CLUSTER_05 | 10 | 2.4% |
| HpGP | SNP_CLUSTER_06 | 65 | 15.7% |
| HpGP | SNP_CLUSTER_07 | 19 | 4.6% |
| HpGP | SNP_CLUSTER_08 | 27 | 6.5% |
| Read cohort | SNP_CLUSTER_01 | 2 | 4.3% |
| Read cohort | SNP_CLUSTER_03 | 44 | 95.7% |

### Table S3b. Near-clone pairs meeting the frozen conjunction

| Isolate pair | Mash distance | ANI (%) | Aligned A | Aligned B |
| --- | --- | --- | --- | --- |
| FRA-012 / FRA-014 | 0.000144 | 99.99 | 100.0% | 99.3% |
| GRE-041 / GRE-046 | 0.000096 | 99.99 | 100.0% | 100.0% |
| POL-106 / POL-108 | 0.000048 | 99.99 | 100.0% | 100.0% |
| TUR-002 / TUR-005 | 0.000459 | 99.92 | 99.9% | 99.9% |

*Note:* All four pairs were within HpGP; no near-clone component crossed an internal/external cohort boundary.

### Table S4. Frozen-panel sample-prediction accounting

| Cohort | Drug | All records | Primary | Pred R | Pred S | Correct | Errors |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Ningxia | Clarithromycin | 60 | 6 | 2 | 4 | 6 | 0 |
| Ningxia | Levofloxacin | 60 | 57 | 19 | 38 | 34 | 23 |
| HpGP | Clarithromycin | 414 | 414 | 108 | 306 | 414 | 0 |
| HpGP | Levofloxacin | 414 | 414 | 131 | 283 | 412 | 2 |
| Read cohort | Clarithromycin | 52 | 44 | 19 | 25 | 41 | 3 |
| Read cohort | Levofloxacin | 52 | 39 | 15 | 24 | 36 | 3 |

*Note:* The full 1,052-row sample-level table is provided in workbook sheet S04_Predictions.

### Table S5. Frozen catalogue performance by cohort

| Cohort | Drug | n (R/S) | Sensitivity (95% CI) | Specificity (95% CI) | FSR (95% CI) | BA |
| --- | --- | --- | --- | --- | --- | --- |
| Ningxia | Clarithromycin | 6 (2/4) | 100.0% (15.8%-100.0%) | 100.0% (39.8%-100.0%) | 0.0% (0.0%-84.2%) | 100.0% |
| Ningxia | Levofloxacin | 57 (28/29) | 42.9% (24.5%-62.8%) | 75.9% (56.5%-89.7%) | 57.1% (37.2%-75.5%) | 59.4% |
| HpGP | Clarithromycin | 414 (108/306) | 100.0% (96.6%-100.0%) | 100.0% (98.8%-100.0%) | 0.0% (0.0%-3.4%) | 100.0% |
| HpGP | Levofloxacin | 414 (133/281) | 98.5% (94.7%-99.8%) | 100.0% (98.7%-100.0%) | 1.5% (0.2%-5.3%) | 99.2% |
| Read cohort | Clarithromycin | 44 (20/24) | 90.0% (68.3%-98.8%) | 95.8% (78.9%-99.9%) | 10.0% (1.2%-31.7%) | 92.9% |
| Read cohort | Levofloxacin | 39 (16/23) | 87.5% (61.7%-98.4%) | 95.7% (78.1%-99.9%) | 12.5% (1.6%-38.3%) | 91.6% |

### Table S6a. Validation-design benchmark: Clarithromycin

| Model | Validation design | Folds | Sensitivity | Specificity | FSR | BA |
| --- | --- | --- | --- | --- | --- | --- |
| Lineage only | Clone Grouped Split | 5 | 0.0% | 100.0% | 100.0% | 50.0% |
| Lineage only | Hpgp To External | 1 | 0.0% | 100.0% | 100.0% | 50.0% |
| Lineage only | Leave Country Out | 1 | 0.0% | 100.0% | 100.0% | 50.0% |
| Lineage only | Leave Lineage Out | 3 | 0.0% | 100.0% | 100.0% | 50.0% |
| Lineage only | Random Isolate Split | 100 | 0.0% | 100.0% | 100.0% | 50.0% |
| Mutation only | Clone Grouped Split | 5 | 100.0% | 100.0% | 0.0% | 100.0% |
| Mutation only | Hpgp To External | 1 | 90.0% | 95.8% | 10.0% | 92.9% |
| Mutation only | Leave Country Out | 1 | 100.0% | 100.0% | 0.0% | 100.0% |
| Mutation only | Leave Lineage Out | 3 | 91.4% | 100.0% | 8.6% | 95.7% |
| Mutation only | Random Isolate Split | 100 | 100.0% | 100.0% | 0.0% | 100.0% |
| Mutation + lineage | Clone Grouped Split | 5 | 94.2% | 100.0% | 5.8% | 97.1% |
| Mutation + lineage | Hpgp To External | 1 | 90.0% | 95.8% | 10.0% | 92.9% |
| Mutation + lineage | Leave Country Out | 1 | 100.0% | 100.0% | 0.0% | 100.0% |
| Mutation + lineage | Leave Lineage Out | 3 | 100.0% | 100.0% | 0.0% | 100.0% |
| Mutation + lineage | Random Isolate Split | 100 | 99.4% | 100.0% | 0.6% | 99.7% |

### Table S6b. Validation-design benchmark: Levofloxacin

| Model | Validation design | Folds | Sensitivity | Specificity | FSR | BA |
| --- | --- | --- | --- | --- | --- | --- |
| Lineage only | Clone Grouped Split | 5 | 2.4% | 97.0% | 97.6% | 49.7% |
| Lineage only | Hpgp To External | 2 | 0.0% | 100.0% | 100.0% | 50.0% |
| Lineage only | Leave Lineage Out | 4 | 0.0% | 100.0% | 100.0% | 50.0% |
| Lineage only | Random Isolate Split | 100 | 3.4% | 97.5% | 96.6% | 50.5% |
| Mutation only | Clone Grouped Split | 5 | 90.8% | 100.0% | 9.2% | 95.4% |
| Mutation only | Hpgp To External | 2 | 65.2% | 85.8% | 34.8% | 75.5% |
| Mutation only | Leave Lineage Out | 4 | 98.0% | 100.0% | 2.0% | 99.0% |
| Mutation only | Random Isolate Split | 100 | 96.0% | 100.0% | 4.0% | 98.0% |
| Mutation + lineage | Clone Grouped Split | 5 | 92.3% | 100.0% | 7.7% | 96.1% |
| Mutation + lineage | Hpgp To External | 2 | 65.2% | 85.8% | 34.8% | 75.5% |
| Mutation + lineage | Leave Lineage Out | 4 | 98.0% | 100.0% | 2.0% | 99.0% |
| Mutation + lineage | Random Isolate Split | 100 | 95.3% | 100.0% | 4.7% | 97.7% |

### Table S7a. Phenotype and breakpoint sensitivity: Clarithromycin

| Scenario | Cohort | n | Sensitivity | Specificity | FSR | BA |
| --- | --- | --- | --- | --- | --- | --- |
| Primary Original I Excluded | Ningxia | 6 | 100.0% | 100.0% | 0.0% | 100.0% |
| Primary Original I Excluded | HpGP | 414 | 100.0% | 100.0% | 0.0% | 100.0% |
| Primary Original I Excluded | Read cohort | 44 | 90.0% | 95.8% | 10.0% | 92.9% |
| Exclude Borderline Mic | Ningxia | 5 | 100.0% | 100.0% | 0.0% | 100.0% |
| Exclude Borderline Mic | HpGP | 414 | 100.0% | 100.0% | 0.0% | 100.0% |
| Exclude Borderline Mic | Read cohort | 44 | 90.0% | 95.8% | 10.0% | 92.9% |
| I As S | Ningxia | 6 | 100.0% | 100.0% | 0.0% | 100.0% |
| I As S | HpGP | 414 | 100.0% | 100.0% | 0.0% | 100.0% |
| I As S | Read cohort | 46 | 90.0% | 92.3% | 10.0% | 91.2% |
| I As R | Ningxia | 6 | 100.0% | 100.0% | 0.0% | 100.0% |
| I As R | HpGP | 414 | 100.0% | 100.0% | 0.0% | 100.0% |
| I As R | Read cohort | 46 | 86.4% | 95.8% | 13.6% | 91.1% |
| Recomputed I Excluded | Ningxia | 5 | 100.0% | 100.0% | 0.0% | 100.0% |
| Recomputed I As S | Ningxia | 6 | 100.0% | 100.0% | 0.0% | 100.0% |
| Recomputed I As R | Ningxia | 6 | 66.7% | 100.0% | 33.3% | 83.3% |

### Table S7b. Phenotype and breakpoint sensitivity: Levofloxacin

| Scenario | Cohort | n | Sensitivity | Specificity | FSR | BA |
| --- | --- | --- | --- | --- | --- | --- |
| Primary Original I Excluded | Ningxia | 57 | 42.9% | 75.9% | 57.1% | 59.4% |
| Primary Original I Excluded | HpGP | 414 | 98.5% | 100.0% | 1.5% | 99.2% |
| Primary Original I Excluded | Read cohort | 39 | 87.5% | 95.7% | 12.5% | 91.6% |
| Exclude Borderline Mic | Ningxia | 53 | 42.9% | 76.0% | 57.1% | 59.4% |
| Exclude Borderline Mic | HpGP | 414 | 98.5% | 100.0% | 1.5% | 99.2% |
| Exclude Borderline Mic | Read cohort | 39 | 87.5% | 95.7% | 12.5% | 91.6% |
| I As S | Ningxia | 57 | 42.9% | 75.9% | 57.1% | 59.4% |
| I As S | HpGP | 414 | 98.5% | 100.0% | 1.5% | 99.2% |
| I As S | Read cohort | 42 | 87.5% | 92.3% | 12.5% | 89.9% |
| I As R | Ningxia | 57 | 42.9% | 75.9% | 57.1% | 59.4% |
| I As R | HpGP | 414 | 98.5% | 100.0% | 1.5% | 99.2% |
| I As R | Read cohort | 42 | 78.9% | 95.7% | 21.1% | 87.3% |
| Recomputed I Excluded | Ningxia | 57 | 42.9% | 75.9% | 57.1% | 59.4% |
| Recomputed I As S | Ningxia | 57 | 42.9% | 75.9% | 57.1% | 59.4% |
| Recomputed I As R | Ningxia | 57 | 42.9% | 75.9% | 57.1% | 59.4% |

### Table S8a. Negative-control completion inventory

| ID | Control | Status | Interpretation/implementation |
| --- | --- | --- | --- |
| 1 | Random Vs Leave Study Out | Computed Elsewhere | random HpGP splits versus HpGP-to-external tests |
| 2 | Random Vs Leave Lineage Out | Computed Elsewhere | same frozen-feature models |
| 3 | Exclude Near Clones | Computed | lexicographic first isolate per connected component |
| 4 | One Isolate Per Patient | Not Estimable | public patient identifiers unavailable |
| 5 | Raw Read Supported Only | Computed Elsewhere | Zenodo raw-read predictions authoritative |
| 6 | Assembly Supported Only | Computed Elsewhere | HpGP and Ningxia assembly predictions; Zenodo support assembly separate |
| 7 | Alternative Breakpoint Interpretation | Computed Elsewhere | original and recomputed Ningxia labels with I sensitivities |
| 8 | Exclude Borderline Mic | Computed Elsewhere | prespecified borderline flag |
| 9 | Ast Method Strata | Computed Elsewhere | method-specific descriptive estimates |
| 10 | Unified Laboratory Ast | Computed Elsewhere | cohort-specific estimates preserve laboratory/study separation |
| 11 | High Quality Genomes Only | Computed Elsewhere | frozen complete genome-QC pass |
| 12 | Country Strata | Computed Elsewhere | country-specific descriptive estimates |
| 13 | Lineage Strata | Computed Elsewhere | predeclared estimability gate plus descriptive sensitivity table |
| 14 | Lineage Only Model | Computed Elsewhere | prespecified negative-control model |
| 15 | Random Core Snp Panels 100 | Computed | HpGP sequence frequency matching without phenotype selection |
| 16 | Training Label Permutation 100 | Computed | HpGP training labels only |

### Table S8b. Near-clone-thinned performance

| Cohort | Drug | n | Sensitivity | Specificity | FSR | BA |
| --- | --- | --- | --- | --- | --- | --- |
| Ningxia | Clarithromycin | 6 | 100.0% | 100.0% | 0.0% | 100.0% |
| Ningxia | Levofloxacin | 57 | 42.9% | 75.9% | 57.1% | 59.4% |
| HpGP | Clarithromycin | 410 | 100.0% | 100.0% | 0.0% | 100.0% |
| HpGP | Levofloxacin | 410 | 98.5% | 100.0% | 1.5% | 99.2% |
| Read cohort | Clarithromycin | 44 | 90.0% | 95.8% | 10.0% | 92.9% |
| Read cohort | Levofloxacin | 39 | 87.5% | 95.7% | 12.5% | 91.6% |

### Table S8c. Raw-read versus support-assembly concordance

| Drug | Raw callable | Assembly callable | Both callable | Concordant | Discordant | Agreement |
| --- | --- | --- | --- | --- | --- | --- |
| Clarithromycin | 52 | 44 | 44 | 44 | 0 | 100.0% |
| Levofloxacin | 45 | 45 | 41 | 41 | 0 | 100.0% |

### Table S8d. Random-panel and label-permutation balanced-accuracy distributions

| Control | Drug | Design | Replicates | Median (2.5th-97.5th percentile) |
| --- | --- | --- | --- | --- |
| Label permutation | Clarithromycin | Hpgp To External Label Permutation | 100 | 50.0% (50.0%-50.0%) |
| Label permutation | Clarithromycin | Random Isolate Label Permutation | 100 | 50.0% (50.0%-50.0%) |
| Label permutation | Levofloxacin | Hpgp To External Label Permutation | 200 | 50.0% (50.0%-51.8%) |
| Label permutation | Levofloxacin | Random Isolate Label Permutation | 100 | 50.0% (50.0%-61.1%) |
| Random SNP panel | Clarithromycin | Hpgp To External Random Snp Panel | 100 | 50.0% (50.0%-50.0%) |
| Random SNP panel | Clarithromycin | Random Isolate Random Snp Panel | 100 | 50.0% (50.0%-50.0%) |
| Random SNP panel | Levofloxacin | Hpgp To External Random Snp Panel | 200 | 50.0% (46.5%-53.1%) |
| Random SNP panel | Levofloxacin | Random Isolate Random Snp Panel | 100 | 50.0% (46.9%-51.9%) |

### Table S9. Ningxia MIC lower-bound summaries

| Drug | Marker prediction | n | Right-censored | Median (mg/L) | IQR (mg/L) | Mann-Whitney P |
| --- | --- | --- | --- | --- | --- | --- |
| Clarithromycin | R | 2 | 0 | 8 | 8-8 | NE |
| Clarithromycin | S | 4 | 0 | 0.09 | 0.06-0.22 | NE |
| Levofloxacin | R | 19 | 7 | 8 | 0.5-32 | NE |
| Levofloxacin | S | 38 | 6 | 0.75 | 0.12-8 | NE |
| Clarithromycin | R_VS_S_TEST | 6 | 0 | NE | NE-NE | 0.095 |
| Levofloxacin | R_VS_S_TEST | 57 | 13 | NE | NE-NE | 0.060 |

### Table S10. Ningxia source-genotype and independent-caller audit

| Audit | Drug | Estimate/call source | Scope | n | Result |
| --- | --- | --- | --- | --- | --- |
| Phenotype Performance | CLR | Source Reported G R | All 60 Source Rows | 60 | Sens 87.1%; spec 100.0%; FSR 12.9% |
| Phenotype Performance | CLR | Source Reported G R | Final Genome Qc Pass | 57 | Sens 90.0%; spec 100.0%; FSR 10.0% |
| Phenotype Performance | CLR | Independent Recall | Final Genome Qc Pass And Target Callable | 6 | Sens 100.0%; spec 100.0%; FSR 0.0% |
| Call Reproducibility | CLR | Source Reported G R Vs Independent Recall | Final Genome Qc Pass And Both Callable | 6 | 6/6 agree (100.0%) |
| Phenotype Performance | LVX | Source Reported G R | All 60 Source Rows | 60 | Sens 89.7%; spec 93.5%; FSR 10.3% |
| Phenotype Performance | LVX | Source Reported G R | Final Genome Qc Pass | 57 | Sens 89.3%; spec 93.1%; FSR 10.7% |
| Phenotype Performance | LVX | Independent Recall | Final Genome Qc Pass And Target Callable | 57 | Sens 42.9%; spec 75.9%; FSR 57.1% |
| Call Reproducibility | LVX | Source Reported G R Vs Independent Recall | Final Genome Qc Pass And Both Callable | 57 | 31/57 agree (54.4%) |
| Independent Caller Reproducibility | LVX | Target Blast Vs Snippy Whole Genome Alignment | Final Genome Qc Pass And Both Callable | 57 | 57/57 agree (100.0%) |

### Table S11. Frozen transportability classification

| Drug | Classification | Frozen-rule reason | Scope |
| --- | --- | --- | --- |
| Clarithromycin | Insufficient Data | fewer than two estimable external cohorts with at least 10 S and 10 R | two Chinese external cohorts; no global claim |
| Levofloxacin | High False Susceptible Risk | at least one eligible external cohort exceeded the frozen 0.10 false-susceptible gate | two Chinese external cohorts; no global claim |

### Table S13a. Phenotype-marker prevalence shift

| Cohort | Drug | n | Phenotypic R | Marker R | Marker minus phenotype |
| --- | --- | --- | --- | --- | --- |
| Ningxia | Levofloxacin | 57 | 49.1% | 33.3% | -15.8 pp |
| Ningxia | Clarithromycin | 6 | 33.3% | 33.3% | +0.0 pp |
| HpGP | Clarithromycin | 414 | 26.1% | 26.1% | +0.0 pp |
| HpGP | Levofloxacin | 414 | 32.1% | 31.6% | -0.5 pp |
| Read cohort | Clarithromycin | 44 | 45.5% | 43.2% | -2.3 pp |
| Read cohort | Levofloxacin | 39 | 41.0% | 38.5% | -2.6 pp |

### Table S13b. Marker prevalence and error by fixed SNP cluster

| Cohort | Drug | Cluster | n | Phenotypic R | Marker R | Error rate |
| --- | --- | --- | --- | --- | --- | --- |
| Ningxia | CLR | C03 | 6 | 33.3% | 33.3% | 0.0% |
| Ningxia | LVX | C01 | 4 | 100.0% | 0.0% | 100.0% |
| Ningxia | LVX | C03 | 52 | 44.2% | 36.5% | 34.6% |
| Ningxia | LVX | C04 | 1 | 100.0% | 0.0% | 100.0% |
| HpGP | CLR | C01 | 122 | 25.4% | 25.4% | 0.0% |
| HpGP | CLR | C02 | 2 | 0.0% | 0.0% | 0.0% |
| HpGP | CLR | C03 | 60 | 38.3% | 38.3% | 0.0% |
| HpGP | CLR | C04 | 109 | 38.5% | 38.5% | 0.0% |
| HpGP | CLR | C05 | 10 | 20.0% | 20.0% | 0.0% |
| HpGP | CLR | C06 | 65 | 6.2% | 6.2% | 0.0% |
| HpGP | CLR | C07 | 19 | 5.3% | 5.3% | 0.0% |
| HpGP | CLR | C08 | 27 | 18.5% | 18.5% | 0.0% |
| HpGP | LVX | C01 | 122 | 29.5% | 29.5% | 0.0% |
| HpGP | LVX | C02 | 2 | 100.0% | 0.0% | 100.0% |
| HpGP | LVX | C03 | 60 | 36.7% | 36.7% | 0.0% |
| HpGP | LVX | C04 | 109 | 25.7% | 25.7% | 0.0% |
| HpGP | LVX | C05 | 10 | 20.0% | 20.0% | 0.0% |
| HpGP | LVX | C06 | 65 | 33.8% | 33.8% | 0.0% |
| HpGP | LVX | C07 | 19 | 36.8% | 36.8% | 0.0% |
| HpGP | LVX | C08 | 27 | 51.9% | 51.9% | 0.0% |
| Read cohort | CLR | C01 | 2 | 50.0% | 50.0% | 0.0% |
| Read cohort | CLR | C03 | 42 | 45.2% | 42.9% | 7.1% |
| Read cohort | LVX | C01 | 2 | 50.0% | 50.0% | 0.0% |
| Read cohort | LVX | C03 | 37 | 40.5% | 37.8% | 8.1% |

*Note:* C01-C08 correspond to the fixed labels SNP_CLUSTER_01-SNP_CLUSTER_08.

### Table S14. Assembly-quality associations with target callability

| Cohort | Drug | Metric | Callable/uncallable | Median callable | Median uncallable | P | Cliff's delta |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Ningxia | CLR | Contigs | 6/54 | 24.5 | 25.5 | 0.815 | -0.062 |
| Ningxia | CLR | Log10 N50 Bp | 6/54 | 5.012 | 5.025 | 0.639 | -0.123 |
| Ningxia | CLR | Completeness Percent | 6/54 | 99.99 | 99.98 | 0.199 | 0.315 |
| Ningxia | CLR | Contamination Percent | 6/54 | 0.01 | 0.045 | 0.499 | -0.17 |
| Ningxia | LVX | Contigs | 57/3 | 25 | 94 | 0.075 | -0.62 |
| Ningxia | LVX | Log10 N50 Bp | 57/3 | 5.025 | 4.313 | 0.061 | 0.649 |
| Ningxia | LVX | Completeness Percent | 57/3 | 99.99 | 91.49 | 0.007 | 0.912 |
| Ningxia | LVX | Contamination Percent | 57/3 | 0.01 | 10.11 | 0.003 | -1 |
| HpGP | CLR | Contigs | 414/0 | 1 | NE | NE | NE |
| HpGP | CLR | Log10 N50 Bp | 414/0 | 6.214 | NE | NE | NE |
| HpGP | CLR | Completeness Percent | 414/0 | 99.99 | NE | NE | NE |
| HpGP | CLR | Contamination Percent | 414/0 | 0.06 | NE | NE | NE |
| HpGP | LVX | Contigs | 414/0 | 1 | NE | NE | NE |
| HpGP | LVX | Log10 N50 Bp | 414/0 | 6.214 | NE | NE | NE |
| HpGP | LVX | Completeness Percent | 414/0 | 99.99 | NE | NE | NE |
| HpGP | LVX | Contamination Percent | 414/0 | 0.06 | NE | NE | NE |
| Read cohort | CLR | Contigs | 46/6 | 80.5 | 423 | 8.15e-05 | -1 |
| Read cohort | CLR | Log10 N50 Bp | 46/6 | 4.624 | 4.106 | 0.135 | 0.384 |
| Read cohort | CLR | Completeness Percent | 46/6 | 99.99 | 95.7 | 0.106 | 0.399 |
| Read cohort | CLR | Contamination Percent | 46/6 | 0.025 | 3.965 | 8.01e-05 | -0.993 |
| Read cohort | LVX | Contigs | 42/10 | 78 | 313 | 0.009 | -0.538 |
| Read cohort | LVX | Log10 N50 Bp | 42/10 | 4.624 | 4.584 | 0.508 | 0.138 |
| Read cohort | LVX | Completeness Percent | 42/10 | 99.99 | 99.97 | 0.318 | 0.2 |
| Read cohort | LVX | Contamination Percent | 42/10 | 0.02 | 2.375 | 0.002 | -0.629 |

### Table S15a. Development-manifold and lineage-composition shift

| Cohort | n | Median 5-NN distance | Median HpGP percentile | P vs HpGP | Cliff's delta | Lineage JS distance |
| --- | --- | --- | --- | --- | --- | --- |
| HpGP | 414 | 0.406 | 0.501 | NE | NE | 0 |
| Ningxia | 57 | 0.135 | 0.058 | 8.82e-19 | -0.723 | 0.72 |
| Read cohort | 46 | 0.116 | 0.022 | 1.06e-19 | -0.816 | 0.773 |

### Table S15b. Development-manifold distance by correct versus incorrect call

| Cohort | Drug | Correct/error n | Median correct | Median error | P | Cliff's delta |
| --- | --- | --- | --- | --- | --- | --- |
| Ningxia | LVX | 34/23 | 0.135 | 0.134 | 0.968 | 0.008 |
| Ningxia | CLR | 6/0 | 0.134 | NE | NE | NE |
| Read cohort | CLR | 41/3 | 0.117 | 0.139 | 0.565 | 0.22 |
| Read cohort | LVX | 36/3 | 0.118 | 0.106 | 0.271 | -0.407 |

### Table S16a. Ningxia levofloxacin false-susceptible decomposition

| Mechanism category | n | Proportion |
| --- | --- | --- |
| Qrdr Wild Type In Deposited Assembly | 15 | 93.8% |
| Off Panel N87Y | 1 | 6.2% |

### Table S16b. Twelve smallest adjusted P values in the full-gyrA missense residual scan

| Variant | Frozen marker | R with/without | S with/without | FN/TP with | Fisher P | BH q |
| --- | --- | --- | --- | --- | --- | --- |
| V655I | False | 21/7 | 29/0 | 10/11 | 0.004 | 0.313 |
| R784K | False | 21/7 | 29/0 | 10/11 | 0.004 | 0.313 |
| G755S | False | 22/6 | 29/0 | 11/11 | 0.010 | 0.484 |
| L817S | False | 2/26 | 10/19 | 1/1 | 0.021 | 0.505 |
| S539N | False | 23/5 | 29/0 | 12/11 | 0.023 | 0.505 |
| V805A | False | 23/5 | 29/0 | 11/12 | 0.023 | 0.505 |
| N797D | False | 21/7 | 28/1 | 10/11 | 0.025 | 0.505 |
| M803V | False | 22/6 | 28/1 | 11/11 | 0.052 | 0.800 |
| N87I | True | 4/24 | 0/29 | 0/4 | 0.052 | 0.800 |
| S492A | False | 21/7 | 27/2 | 12/9 | 0.079 | 0.800 |
| K694R | False | 5/23 | 1/28 | 4/1 | 0.102 | 0.800 |
| P484Q | False | 23/5 | 28/1 | 12/11 | 0.102 | 0.800 |

*Note:* All 140 tested substitution rows are retained in workbook sheet S16_Residual_Mechanisms; none passed false-discovery correction.

### Table S16c. Exploratory add-only N87Y scenario

| Panel | n | TP/TN/FP/FN | Sensitivity | Specificity | FSR | BA |
| --- | --- | --- | --- | --- | --- | --- |
| Frozen Catalogue | 57 | 12/22/7/16 | 42.9% | 75.9% | 57.1% | 59.4% |
| Exploratory Add N87Y | 57 | 13/22/7/15 | 46.4% | 75.9% | 53.6% | 61.1% |

### Table S17a. Near-clone-group bootstrap distributions

| Cohort | Drug | Replicates | Sensitivity | Specificity | FSR | BA |
| --- | --- | --- | --- | --- | --- | --- |
| Ningxia | CLR | 1809 | 100.0% (100.0%-100.0%) | 100.0% (100.0%-100.0%) | 0.0% (0.0%-0.0%) | 100.0% (100.0%-100.0%) |
| Ningxia | LVX | 2000 | 42.9% (24.1%-61.3%) | 76.5% (59.4%-90.0%) | 57.1% (38.7%-75.9%) | 59.4% (47.4%-71.6%) |
| Read cohort | CLR | 2000 | 90.5% (75.0%-100.0%) | 96.0% (86.4%-100.0%) | 9.5% (0.0%-25.0%) | 93.2% (84.1%-100.0%) |
| Read cohort | LVX | 2000 | 88.2% (68.8%-100.0%) | 95.8% (86.4%-100.0%) | 11.8% (0.0%-31.2%) | 92.1% (81.2%-100.0%) |

*Note:* Cells report median (2.5th-97.5th percentile).

### Table S17b. Ningxia levofloxacin leave-one-out influence

| Analysis | FSR | BA |
| --- | --- | --- |
| Baseline | 57.1% | 59.4% |
| Leave-one-out range | 55.6%-59.3% | 58.3%-60.7% |

### Table S17c. Predictive values under selected assumed prevalences

| Cohort | Drug | Assumed resistance prevalence | PPV | NPV |
| --- | --- | --- | --- | --- |
| Ningxia | CLR | 10.0% | 100.0% | 100.0% |
| Ningxia | CLR | 25.0% | 100.0% | 100.0% |
| Ningxia | CLR | 50.0% | 100.0% | 100.0% |
| Ningxia | LVX | 10.0% | 16.5% | 92.3% |
| Ningxia | LVX | 25.0% | 37.2% | 79.9% |
| Ningxia | LVX | 50.0% | 64.0% | 57.0% |
| Read cohort | CLR | 10.0% | 70.6% | 98.9% |
| Read cohort | CLR | 25.0% | 87.8% | 96.6% |
| Read cohort | CLR | 50.0% | 95.6% | 90.6% |
| Read cohort | LVX | 10.0% | 69.1% | 98.6% |
| Read cohort | LVX | 25.0% | 87.0% | 95.8% |
| Read cohort | LVX | 50.0% | 95.3% | 88.4% |

*Note:* These are mathematical prevalence scenarios, not empirical recalibration.

### Table S17d. Transport-gate robustness-grid classifications

| Drug | Classification | Grid cells |
| --- | --- | --- |
| Clarithromycin | High False Susceptible Risk | 1 |
| Clarithromycin | Insufficient Data | 32 |
| Clarithromycin | Passes Exploratory Grid | 7 |
| Levofloxacin | High False Susceptible Risk | 32 |
| Levofloxacin | Insufficient Data | 8 |

*Note:* The prespecified rule remained the primary classification; this grid is a post-freeze robustness display.

### Table S18a. Four-domain transport-shift summary

| Domain | Cohort | Drug | Measure | Estimate | Detail |
| --- | --- | --- | --- | --- | --- |
| Analytic Availability | Ningxia | CLR | Target Callability | 10.0% | 6/60 |
| Analytic Availability | Ningxia | LVX | Target Callability | 95.0% | 57/60 |
| Analytic Availability | Read cohort | CLR | Target Callability | 88.5% | 46/52 |
| Analytic Availability | Read cohort | LVX | Target Callability | 80.8% | 42/52 |
| Population Shift | Ningxia | - | Pc1-Pc10 Cohort Discriminator Oof Auc | 89.9% | permutation P=0.000999 |
| Population Shift | Read cohort | - | Pc1-Pc10 Cohort Discriminator Oof Auc | 91.8% | permutation P=0.000999 |
| Conditional Transport | Ningxia | LVX | Snp Cluster 03 Sensitivity | 52.2% | 12/23 |
| Conditional Transport | Ningxia | LVX | Snp Cluster 03 Specificity | 75.9% | 22/29 |
| Conditional Transport | Ningxia | LVX | Snp Cluster 03 Marker Negative Resistance | 33.3% | 11/33 |
| Conditional Transport | HpGP | LVX | Snp Cluster 03 Sensitivity | 100.0% | 22/22 |
| Conditional Transport | HpGP | LVX | Snp Cluster 03 Specificity | 100.0% | 38/38 |
| Conditional Transport | HpGP | LVX | Snp Cluster 03 Marker Negative Resistance | 0.0% | 0/38 |
| Conditional Transport | Read cohort | LVX | Snp Cluster 03 Sensitivity | 86.7% | 13/15 |
| Conditional Transport | Read cohort | LVX | Snp Cluster 03 Specificity | 95.5% | 21/22 |
| Conditional Transport | Read cohort | LVX | Snp Cluster 03 Marker Negative Resistance | 8.7% | 2/23 |
| Clinical Severity | Ningxia | LVX | False Susceptible Median Lower Bound Mic Mg L | 16 | n=16; right-censored=6 |

### Table S18b. Dominant-lineage conditional comparisons

| Metric | Dataset A | Dataset B | A success/fail | B success/fail | Odds ratio | Fisher P |
| --- | --- | --- | --- | --- | --- | --- |
| Sensitivity | Ningxia | HpGP | 12/11 | 22/0 | 0 | 2.03e-04 |
| Specificity | Ningxia | HpGP | 22/7 | 38/0 | 0 | 0.002 |
| Marker Negative Resistance | Ningxia | HpGP | 11/22 | 0/38 | NE | 7.56e-05 |
| Sensitivity | Ningxia | Read cohort | 12/11 | 13/2 | 0.168 | 0.039 |
| Specificity | Ningxia | Read cohort | 22/7 | 21/1 | 0.15 | 0.117 |
| Marker Negative Resistance | Ningxia | Read cohort | 11/22 | 2/21 | 5.25 | 0.052 |

### Table S18c. Phenotype-blind cohort-discriminator audit

| External cohort | HpGP/external n | Features | OOF AUC | Group-bootstrap 95% CI | Permutation P |
| --- | --- | --- | --- | --- | --- |
| Ningxia | 414/57 | PC1-PC10 | 89.9% | 87.1%-92.4% | 9.99e-04 |
| Read cohort | 414/46 | PC1-PC10 | 91.8% | 89.3%-94.1% | 9.99e-04 |

### Table S18d. Ningxia levofloxacin resistant-MIC severity

| Outcome | n | Median lower-bound MIC (mg/L) | IQR (mg/L) | Right-censored |
| --- | --- | --- | --- | --- |
| TP | 12 | 32 | 14-32 | 7 |
| FN | 16 | 16 | 3.5-32 | 6 |

### Table S19a. Frozen 23S configuration and recovery-class summary

| Section | Metric | Value | Denominator | Note |
| --- | --- | --- | --- | --- |
| Configuration | Frozen Blast Task | megablast | NE | primary caller |
| Configuration | Frozen Minimum Identity | 0.9 | NE | fraction |
| Configuration | Frozen Minimum Query Coverage | 0.05 | NE | fraction |
| Result | Frozen Callable Final Qc | 6 | 57 | unchanged primary callability |
| Recovery Class | Callable Marker Spanning | 7 | 60 | all Ningxia assemblies |
| Recovery Class | Partial 23S No Marker Span | 53 | 60 | all Ningxia assemblies |

### Table S19b. 23S callability sensitivity-grid audit

| BLAST task | Grid cells | Identity range | Coverage range | Final-QC callable counts |
| --- | --- | --- | --- | --- |
| blastn | 30 | 70.0%-98.0% | 0.0%-90.0% | 6 |
| megablast | 30 | 70.0%-98.0% | 0.0%-90.0% | 6 |

*Note:* Every grid setting required both resistance-marker bases to be spanned; no threshold relaxation rescued an additional final-QC assembly.

### Table S20a. End-to-end correct-result yield

| Cohort | Drug | Binary phenotype | QC pass | Primary evaluable | Correct | Unresolved | Correct yield |
| --- | --- | --- | --- | --- | --- | --- | --- |
| HpGP | CLR | 414 | 414 | 414 | 414 | 0 | 100.0% |
| HpGP | LVX | 414 | 414 | 414 | 412 | 0 | 99.5% |
| Ningxia | CLR | 60 | 57 | 6 | 6 | 54 | 10.0% |
| Ningxia | LVX | 60 | 57 | 57 | 34 | 3 | 56.7% |
| Read cohort | CLR | 50 | 44 | 44 | 41 | 6 | 82.0% |
| Read cohort | LVX | 48 | 43 | 39 | 36 | 9 | 75.0% |

*Note:* Correct yield uses every phenotype-linked binary isolate as the denominator; it is not accuracy conditional on receiving a result.

### Table S20b. Callability-aware logical performance bounds

| Cohort | Drug | Metric | Resolved/total | Evaluable estimate | Logical lower-upper bound |
| --- | --- | --- | --- | --- | --- |
| HpGP | CLR | Sensitivity | 108/108 | 100.0% | 100.0%-100.0% |
| HpGP | CLR | Specificity | 306/306 | 100.0% | 100.0%-100.0% |
| HpGP | LVX | Sensitivity | 133/133 | 98.5% | 98.5%-98.5% |
| HpGP | LVX | Specificity | 281/281 | 100.0% | 100.0%-100.0% |
| Ningxia | CLR | Sensitivity | 2/31 | 100.0% | 6.5%-100.0% |
| Ningxia | CLR | Specificity | 4/29 | 100.0% | 13.8%-100.0% |
| Ningxia | LVX | Sensitivity | 28/29 | 42.9% | 41.4%-44.8% |
| Ningxia | LVX | Specificity | 29/31 | 75.9% | 71.0%-77.4% |
| Read cohort | CLR | Sensitivity | 20/23 | 90.0% | 78.3%-91.3% |
| Read cohort | CLR | Specificity | 24/27 | 95.8% | 85.2%-96.3% |
| Read cohort | LVX | Sensitivity | 16/20 | 87.5% | 70.0%-90.0% |
| Read cohort | LVX | Specificity | 23/28 | 95.7% | 78.6%-96.4% |

*Note:* Bounds assign unresolved isolates to the least or most favourable compatible prediction; they are not confidence intervals.

### Table S20c. Phenotype-blind development-manifold abstention for levofloxacin

| Cohort | HpGP distance cutoff | Accepted/primary | Coverage | Sensitivity | Specificity | FSR | Passes gate |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Ningxia | 5.0% | 26/57 | 45.6% | 36.4% | 80.0% | 63.6% | no |
| Ningxia | 10.0% | 37/57 | 64.9% | 47.1% | 70.0% | 52.9% | no |
| Ningxia | 20.0% | 50/57 | 87.7% | 50.0% | 75.0% | 50.0% | no |
| Ningxia | 50.0% | 50/57 | 87.7% | 50.0% | 75.0% | 50.0% | no |
| Ningxia | 90.0% | 56/57 | 98.2% | 44.4% | 75.9% | 55.6% | no |
| Ningxia | 100.0% | 57/57 | 100.0% | 42.9% | 75.9% | 57.1% | no |
| Read cohort | 5.0% | 29/39 | 74.4% | 83.3% | 94.1% | 16.7% | no |
| Read cohort | 10.0% | 34/39 | 87.2% | 84.6% | 95.2% | 15.4% | no |
| Read cohort | 20.0% | 36/39 | 92.3% | 85.7% | 95.5% | 14.3% | no |
| Read cohort | 50.0% | 37/39 | 94.9% | 86.7% | 95.5% | 13.3% | no |
| Read cohort | 90.0% | 39/39 | 100.0% | 87.5% | 95.7% | 12.5% | no |
| Read cohort | 100.0% | 39/39 | 100.0% | 87.5% | 95.7% | 12.5% | no |

*Note:* Cutoffs were fixed HpGP empirical distance percentiles and did not use external phenotypes; a subset also required at least 10 resistant and 10 susceptible isolates to pass the frozen gate.

### Table S20d. External-cohort levofloxacin performance differences

| Metric | Ningxia | Read cohort | Absolute difference | 95% interval | Fisher P |
| --- | --- | --- | --- | --- | --- |
| Sensitivity | 42.9% | 87.5% | -44.6 pp | -63.3 to -15.0 pp | 0.005 |
| Specificity | 75.9% | 95.7% | -19.8 pp | -38.1 to +0.7 pp | 0.064 |
| False Susceptible Rate | 57.1% | 12.5% | +44.6 pp | +15.0 to +63.3 pp | 0.005 |
| Balanced Accuracy | 59.4% | 91.6% | -32.2 pp | -47.1 to -17.4 pp | NE |

*Note:* Differences are Ningxia minus read cohort. Balanced-accuracy uncertainty uses the near-clone-group bootstrap and has no single Fisher exact P value.
