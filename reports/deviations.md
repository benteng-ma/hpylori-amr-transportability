# Deviations from the launch assumptions

The Phase 1 plan was frozen on 2026-08-30. Post-freeze execution corrections are recorded below; none changed cohorts, phenotypes, frozen panels, validation gates, or estimands.

Pre-freeze data-audit corrections:

- US Mayo 2025 was initially listed as the anticipated primary multi-drug external cohort. It was made ineligible after the public supplement was shown to contain only aggregate AST and no isolate-to-assembly crosswalk.
- Linqu 2025 was changed from a phenotype stress-test candidate to a sequence-only lineage resource because its public Zenodo archive contains assemblies but no public isolate-level phenotype crosswalk.
- Israel PRJEB37854 was changed from an external candidate to ineligible because both ENA filereport and NCBI SRA queries returned zero records on 2026-08-30.
- Ningxia 2022 became the primary external cohort after exact 60/60 MIC-to-assembly linkage and independent extraction of the publication's 60-row P-R table.
- Zenodo 10369064 became a secondary external phenotype-method stress test after exact 52/52 read-to-label mapping, with its custom and internally overlapping levofloxacin disk thresholds explicitly retained as a limitation.

These changes were made before performance inspection and therefore do not constitute result-driven cohort selection.

## 2026-08-30 — Phase 2 acquisition implementation

- Added a resumable, checksum-verified downloader for the already frozen 474 public assemblies and 52 Zenodo paired-read archives.
- Large Zenodo files may use a native-curl transport wrapper with identical frozen task construction, resume behavior, expected-byte checks, repository MD5 checks, full ZIP validation, and recorded SHA-256; this changes transfer throughput only.
- The repository-supplied 52-sample merged SNP table was acquired as a supporting population-structure resource and verified against the Zenodo byte count and MD5 plus a locally recorded SHA-256. It cannot replace the raw-read marker calls or create performance observations.
- This is an implementation specification, not a change to cohorts, phenotypes, markers, QC thresholds, estimands, or validation gates.
- No performance output had been inspected when this implementation was added.
- Basic sequence QC materializes the frozen public assemblies and applies the already frozen genome-size and contig-count gates. Completeness and contamination remain explicitly `NOT_COMPUTED` until an independent marker-based genome-quality tool is added; they are not silently inferred from NCBI submission status.
- The Phase 0 targeted caller is reused without changing marker definitions; the Phase 2 wrapper only adds parallel execution and explicit dataset/isolate provenance.
- Read-derived Zenodo assemblies enter the same size, contig, completeness, contamination, relatedness, and core-SNP QC tables, but their panel predictions remain raw-read authoritative; assembly marker calls are retained only as support and never duplicated in performance denominators.
- Full relatedness uses the same frozen Mash/skani near-clone conjunction as Phase 0. A permissive Mash distance of 0.005 is used only to decide which pairs receive the more expensive skani audit; it is not the near-clone threshold and cannot create an edge unless the frozen 0.001/99.9%/0.90 conjunction is satisfied.
- Full Mash/skani relatedness is restricted to assemblies passing the complete frozen genome-QC path; failed or unestimated genomes cannot define near-clone components used in validation.

## 2026-08-30 — raw-read implementation specification

- Before formal performance inspection, fixed fastp processing to paired-end adapter detection, tail-window mean quality 20, minimum length 50, and overlap correction. These settings do not use phenotype labels.
- Added SKESA assembly of the 52 public read sets so they can enter the prespecified genome-QC, Mash/skani, and core-SNP workflows. This does not replace raw-read marker evidence.
- SKESA is isolated in its own environment so adding the assembler cannot downgrade the frozen principal analysis toolchain.
- The read-processing worker count was increased from two to four after three samples completed successfully. Per-sample tools, thread counts, thresholds, and outputs were unchanged; this is a throughput-only execution correction.
- Processing of the final high-depth read archives exposed severe mounted-drive temporary-I/O contention during BAM sorting. Four incomplete samples were interrupted before producing calls and restarted with execution-only FASTQ/BAM/SKESA intermediates on the native WSL temporary filesystem. All formal fastp reports, mapping logs, pileups, predictions, assemblies, and completion records remain on E drive; no tool option, threshold, read, or result was changed.
- The same four high-depth alignments exposed a second scale-only defect: BWA progress messages were captured in a fixed-size stderr pipe that was read only after `samtools sort`, allowing the pipe to fill and pause both processes. The live pipes were drained without terminating or altering the alignments, and subsequent invocations write BWA stderr directly to each formal per-sample log. This changes logging transport only; reads, alignments, callers, thresholds, and phenotypes are unchanged. The recovered pipe fragments are retained with the execution logs.
- Three high-depth SKESA support assemblies then overlapped in memory and the WSL 15-GB limit terminated S56 before an assembly or completion record was produced. S56 was designated for complete rerun under identical parameters. Subsequent SKESA graph-construction steps are serialized with an execution lock, and a temporary 24-GB WSL swap file was enabled only to prevent operating-system termination. This is a resource-scheduling correction: raw-read marker calls remain authoritative, and no sequence, assembly parameter, QC threshold, phenotype, or estimand changed.
- Raw-read 23S evidence is reported in two complementary forms: unique whole-genome mappings for copy-specific audit and a single-reference pooled 23S alignment for allele fractions. The frozen base/mapping quality 20, depth 10, major fraction 0.80, and mixed fraction 0.20–0.80 gates apply.
- A pooled resistant 23S allele fraction of at least 0.20 yields a resistant panel call; a susceptible call requires both pooled target positions to pass depth and major-allele gates with resistant fraction below 0.20. Copy-specific low coverage remains visible and is never converted to wild type.
- A raw-read gyrA call requires all three bases of residues 87, 88, and 91 to pass the frozen depth and major-allele gates. Any uncallable codon makes the panel prediction uncallable rather than susceptible.
- Assembly-based 23S susceptibility requires at least one target-spanning copy/hit and no frozen resistance marker; complete absence or repeat breakage remains uncallable. This preserves the predeclared assembly-collapse limitation instead of pretending that two copies were resolved.
- Assembly-based gyrA susceptibility requires callable residues 87, 88, and 91 and no frozen marker. Diagnostic denominators exclude uncallable predictions but callability/attrition is reported separately.
- Phenotype sensitivities are fixed as original S/R with I excluded, exclusion of prespecified borderline MICs, I-as-S, and I-as-R. Ningxia MIC summaries retain censoring operators and are explicitly cohort-limited lower-bound rank analyses, never a full multi-cohort continuous endpoint.
- Amoxicillin, furazolidone, metronidazole, and tetracycline are assigned `INSUFFICIENT_DATA` for rule transportability because no exactly reproducible panel was frozen for them; no post hoc marker discovery is used to manufacture a secondary result.
- Classification priority is frozen as: insufficient estimability, high false-susceptible risk, robust transportability, lineage dependence, phenotype-method sensitivity, study dependence, and otherwise non-transportability. This order prevents a safety failure from being hidden by a more favorable descriptive label.

## 2026-08-30 — leakage benchmark implementation

- Before performance inspection, fixed HpGP as the development cohort; external cohorts never contribute to feature selection or hyperparameter tuning.
- Fixed L2 logistic regression with C=1, liblinear, threshold 0.5, and no tuning for mutation-only, lineage-only, and mutation-plus-lineage models. Mutation features are limited to the frozen panel.
- Fixed 100 stratified 80/20 random-isolate repeats, five-fold near-clone-grouped validation, phenotype-eligible leave-country-out and leave-core-SNP-cluster-out folds, and direct HpGP-to-each-external-cohort tests.
- A study-only model is not fitted because HpGP is a single development study and the two external cohorts cannot be used to estimate a study effect without contaminating the external test.
- Prespecified logistic models emit probabilities as well as threshold-0.5 classifications. AUROC, AUPRC, Brier score, and logistic recalibration intercept/slope use those probabilities; deterministic frozen panels retain explicitly labeled binary-score summaries and are not assigned pseudo-calibration.

## 2026-08-30 — phenotype-blind lineage implementation

- Before full sequence processing and any performance inspection, fixed Snippy 4.6.0 against 26695 for the equivalent core-SNP population-structure path, with mapping/base quality 20 and minimum depth 10.
- In Snippy's assembly mode, contigs are deterministically shredded to 20x pseudo-reads and remapped to 26695. Mapping quality 20, depth 10, and VCF quality 20 gates are active; pseudo-read base qualities are fixed high by Snippy, so the base-quality-20 option is mechanically satisfied and is not presented as original-read quality evidence.
- Frozen eight de novo genetic clusters from a 20-component sparse core-SNP embedding and K-means (50 starts; seed 20260830). These labels are used only as phenotype-blind leakage groups and are not presented as formal hp population assignments unless independently concordant with published metadata.
- Mash/skani remains the independent whole-genome distance path; core-SNP clustering is the second population-structure method required by the protocol.
- Before running the full core alignment, hardened the leakage boundary: variable-site selection, SVD fitting, and K-means fitting use HpGP development genomes only. Ningxia and Zenodo genomes are transformed in that fixed feature space and assigned to the nearest fixed centroid, so external genomes do not determine the population clusters used to evaluate them.
- The first single-isolate Snippy smoke test exposed two execution issues before any phenotype performance was inspected: helper-script paths containing spaces were not quoted internally, and a complex non-SNP FreeBayes event failed reference normalization. The formal runtime is therefore materialized under a no-space WSL temporary path from the frozen environment specification, and FreeBayes is explicitly restricted to SNP observations (`-i -X -u --haplotype-length 0`). This changes neither the reference, frozen quality gates, samples, nor phenotype-independent feature selection.
- Per-isolate Snippy working directories are staged under no-space WSL temporary storage to avoid severe mounted-drive small-I/O overhead. The VCF, aligned FASTA, raw/filtered VCF, BAM/index, BED, tabular calls, parameters, and log are copied back to the formal E-drive sample directory before the temporary directory is removed; `snippy-core` consumes only those committed-back artifacts.
- Per-isolate phenotype-blind Snippy calls may be precomputed while CheckM2 is running, but no provisional core alignment is built in that mode. The final `snippy-core` input list is regenerated only from assemblies passing the completed CheckM2-integrated QC table.

## 2026-08-30 — negative-control implementation

- Before performance inspection, operationalized 100 HpGP-training-label permutations, 100 phenotype-blind HpGP-frequency-matched random core-SNP panels, lexicographic one-per-near-clone thinning, and raw-read-versus-support-assembly concordance.
- Random SNP panels match the frozen panel feature count and each frozen marker's HpGP prevalence. Candidate sites, allele frequencies, and selection use HpGP sequences without phenotype labels; external sequences are projected only after selection.
- Patient thinning remains explicitly `NOT_ESTIMABLE` because public patient identifiers are unavailable; near-clone connected components are the frozen fallback, not a claim of patient identity.

## 2026-08-30 — genome completeness and contamination implementation

- Before performance inspection, selected CheckM2 1.1.0 in an isolated environment to operationalize the frozen completeness >=90% and contamination <=10% gates.
- The current Bioconda 1.1.0 build declares Python >3.12, so the isolated environment leaves Python resolution to the signed channel metadata rather than imposing an incompatible legacy pin.
- CheckM2 is constrained to TensorFlow 2.17 CPU builds; an automatically proposed CUDA transaction was cancelled before installation because this host has no need for GPU libraries and the prediction algorithm is unchanged.
- Because the project path contains spaces, CheckM2 is invoked through its Python module rather than the generated shell entry-point; the latter does not quote the environment path correctly. The downloaded database directory is deterministically resolved to its version-validated `uniref100.KO.1.dmnd` file. These are execution-path corrections only.
- CheckM2 1.1.0's official test run further showed that its internal Prodigal and DIAMOND commands interpolate paths into unquoted shell strings. Formal inputs, outputs, and a symlink to the checksum-verified database are therefore staged under a no-space system temporary directory, after which the complete CheckM2 output is copied back to the project. The first official test attempt failed before producing estimates and was retained in `results/logs/checkm2_testrun.log`.
- The CheckM2 downloader could not reach DOI resolution from WSL on this host, so the exact database archive referenced by CheckM2 1.1.0 (Zenodo 14897628, v3) is instead transferred with resumable native curl, checked against the repository byte count and MD5, assigned a local SHA-256, safely extracted, and then subjected to CheckM2's own database-version validation.
- CheckM2 uses its documented low-memory DIAMOND mode to cap memory use; this changes DIAMOND block size and runtime, not the completeness or contamination models or frozen gates.
- Basic size/contig QC remains separately visible. NCBI species labels are not used as a substitute for completeness or contamination estimates.
