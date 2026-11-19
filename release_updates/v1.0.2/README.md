# QC-conditioned R1 correction

Version 1.0.2 accompanies the revised quality-control and target-callability accounting in Tables S2 and S6 and Supplementary Figure S1.

## Correction and data

The assembly-quality association analysis now excludes final-genome-QC failures before comparing callable and uncallable targets. It retains the original two-sided Mann-Whitney U test with automatic exact/asymptotic selection and continuity correction, Cliff's delta direction, descriptive AUC and the minimum of three available records per group.

The 1,052 sample-drug input records are in results/extended_analysis/qc_callability_samples.csv. The input values are preserved, and 1,034 QC-passing records enter the revised comparisons. Target callability for paired reads comes from the primary raw-read calls, while assembly-quality metrics describe the support assemblies.

| Cohort and drug | Linked | QC pass | QC fail | Callable after QC | Uncallable after QC |
| --- | ---: | ---: | ---: | ---: | ---: |
| Ningxia CLR | 60 | 57 | 3 | 6 | 51 |
| Ningxia LVX | 60 | 57 | 3 | 57 | 0 |
| HpGP CLR | 414 | 414 | 0 | 414 | 0 |
| HpGP LVX | 414 | 414 | 0 | 414 | 0 |
| Paired reads CLR | 52 | 46 | 6 | 46 | 0 |
| Paired reads LVX | 52 | 46 | 6 | 42 | 4 |

Conditional callability uses QC-passing records as its denominator. End-to-end callability uses all linked records. Groups without QC-passing uncallable records have no estimable within-QC association.

## Files

- Table_S2_QC_accounting.csv contains both QC stages and conditional/end-to-end callability.
- Table_S6_QC_conditioned.csv contains all 24 metric summaries, with empty statistics for non-estimable comparisons.
- Figure_S1_source_data.csv contains the eight estimable effect-size rows.
- Supplementary_Figure_S1_QC_conditioned.pdf and .png are the two-panel data-driven figure included in the supplement.
- MGEN_Supplementary_Data.xlsx matches the R1 supplementary workbook, including the sample-level inputs, revised summaries and visible S06_Reanalysis_Code.
- reproduce_qc_callability.py regenerates the tables, figure data and figure from the CSV or workbook.
- reproduce_qc_tables_from_workbook.py reproduces the tables and figure data directly from the workbook. The same source is in its code worksheet.
- analysis_environment.json records the environment, original run checksum and archived input checksum. The archived input has LF line endings and the original run input has CRLF, with identical records and field values.

## Reproduction

Run the following commands from the repository root.

~~~sh
python release_updates/v1.0.2/reproduce_qc_callability.py --input results/extended_analysis/qc_callability_samples.csv --output outputs/v1.0.2
python release_updates/v1.0.2/reproduce_qc_tables_from_workbook.py --input release_updates/v1.0.2/MGEN_Supplementary_Data.xlsx --output outputs/v1.0.2_workbook
~~~

The recorded environment uses Python 3.13.11, NumPy 2.3.5, pandas 2.3.3, SciPy 1.16.3, scikit-learn 1.8.0 and Matplotlib 3.10.8. Workbook input additionally requires openpyxl.

The release also corrects the descriptive Z52_LEV breakpoint note. The published levofloxacin zone definitions are R <13 mm, I 13-17 mm and S >17 mm. The metadata preserves the supplied source categories and reports the source's conflicting assay descriptions.

Code is MIT licensed and author-generated compact derived data are CC BY 4.0. Historical release directories preserve their original recorded material.
