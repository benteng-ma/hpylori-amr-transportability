The current reader-facing supplementary workbook is MGEN_Supplementary_Data.xlsx.

Version 1.0.4 removes redundant workbook inventory rows and administrative fields, omits an unpopulated collection-year field, and clarifies the auxiliary-benchmark and external-cohort scope labels. One administrative status now records that no common-method laboratory AST comparison was performed. Scientific observations, predictions, numerical results, Tables S2/S6 and Figure S1 data are preserved. The historical workbooks remain available in their original release directories.

Use the existing reproduction script from the repository root with this workbook:

~~~sh
python release_updates/v1.0.2/reproduce_qc_tables_from_workbook.py --input release_updates/v1.0.4/MGEN_Supplementary_Data.xlsx --output outputs/v1.0.4
~~~

The QC-conditioned CSV inputs, scripts, source-data tables and Figure S1 PDF/PNG remain in release_updates/v1.0.2.
