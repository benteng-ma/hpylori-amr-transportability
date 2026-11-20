# Matching reader workbook and reproducibility notes

The current reader-facing supplementary workbook is `MGEN_Supplementary_Data.xlsx`. Version 1.0.5 clarifies its README cell B4: the 526 linked isolate records comprise 414 HpGP assemblies, 60 Ningxia assemblies and 52 paired-read datasets. All data sheets, numerical values, formulas and executable worksheet content are unchanged from v1.0.4.

## Reproduction and checksums

Use the existing scripts from the repository root:

~~~sh
python release_updates/v1.0.2/reproduce_qc_tables_from_workbook.py --input release_updates/v1.0.5/MGEN_Supplementary_Data.xlsx --output outputs/v1.0.5_workbook
python release_updates/v1.0.2/reproduce_qc_callability.py --input results/extended_analysis/qc_callability_samples.csv --output outputs/v1.0.5_csv
~~~

The [v1.0.2 directory](../v1.0.2/README.md) retains the reproduction scripts, recorded environment, revised Tables S2/S6, Figure S1 PDF/PNG and figure source data. Historical workbook versions remain in their original directories.

The archived CSV input is `results/extended_analysis/qc_callability_samples.csv`. Its LF-form SHA-256 is `49eb8aa893a5ffbfe1d4173c8f6696881902bce370913a303d94e2d16fb111b6`. The original analysis input used CRLF line endings and had SHA-256 `b0464ebb4df0e49db5db0e71322a78e0b0a853a581cc3d073a3beae57dbc4758`. These files contain the same 1,052 parsed records and field values. The historical `release_updates/v1.0.2/analysis_environment.json` records these CSV checksums.

The same sample input is embedded in sheet `S06_QC_Associations`: select rows whose `source_file` is `results/extended_analysis/qc_callability_samples.csv`. The workbook itself has a separate checksum, recorded under `release_updates/v1.0.5/MGEN_Supplementary_Data.xlsx` in the release-root `SHA256SUMS.csv`. When the workbook reproduction script is run, its new output `analysis_environment.json` records the SHA-256 of the actual XLSX supplied through `--input`; this is a workbook checksum, not the historical CSV checksum.

## Database provenance

CheckM2 1.1.0 used the CM2 Reference Database Uniref100/KO version 3 (20 February 2025; [10.5281/zenodo.14897628](https://doi.org/10.5281/zenodo.14897628)). The source archive, size and checksums are retained in `metadata/phase2/tool_database_manifest.csv`. This database DOI identifies the reference database rather than this project's code archive.

## Primary paired-read clarithromycin rule

Both pooled 23S marker positions require depth at least 10, with minimum base and mapping qualities of 20. The resistance check comes first: an A2142 combined G+C fraction or an A2143 G fraction at least 0.20 yields a resistant call. If both resistance fractions are below 0.20, susceptibility additionally requires a major-base fraction at least 0.80 at both sites; otherwise the result is unresolved. At adequate depth, the site-level `MIXED` annotation means a major-base fraction below 0.80 and is separate from the primary prediction. Copy-specific pileups are retained for audit and do not override the pooled classification. Final-genome-QC failures remain uncallable in the benchmark.

These descriptions document the existing implementation in `scripts/process_zenodo_reads.py` and its benchmark integration. Classification logic, fixed predictions and statistical results have not changed.
