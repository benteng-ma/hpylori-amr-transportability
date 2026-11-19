#!/usr/bin/env python3
"""Reproduce QC-pass Tables S2/S6 and Figure S1 data from the companion workbook."""
import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import openpyxl
import scipy
import sklearn
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

COHORTS = ["CHINA_NINGXIA_2022", "HPGP_GLOBAL", "ZENODO_10369064"]
DRUGS = ["clarithromycin", "levofloxacin"]
METRICS = ["contigs", "log10_n50_bp", "completeness_percent", "contamination_percent"]
SOURCE = "results/extended_analysis/qc_callability_samples.csv"


def read_samples(path):
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = book["S06_QC_Associations"].iter_rows(values_only=True)
    header = next(rows)
    data = [dict(zip(header, row)) for row in rows if row[0] == SOURCE]
    assert len(data) == 1052
    assert len({(r["dataset_id"], r["isolate_id"], r["antibiotic"]) for r in data}) == 1052
    for row in data:
        for field in METRICS + ["callable_binary"]:
            row[field] = float(row[field])
        assert row["final_qc_status"] in {"PASS", "FAIL"}
        assert row["callable_binary"] == (row["callable"] == "yes")
        assert row["final_qc_status"] == "PASS" or row["callable_binary"] == 0
        assert all(np.isfinite(row[field]) for field in METRICS)
    return data


def summarize(data):
    accounting, associations = [], []
    for cohort in COHORTS:
        for drug in DRUGS:
            linked = [r for r in data if r["dataset_id"] == cohort and r["antibiotic"] == drug]
            eligible = [r for r in linked if r["final_qc_status"] == "PASS"]
            n_call = sum(int(r["callable_binary"]) for r in eligible)
            accounting.append(dict(
                dataset_id=cohort, antibiotic=drug, n_linked=len(linked),
                n_qc_pass=len(eligible), n_qc_fail=len(linked) - len(eligible),
                n_callable_after_qc=n_call, n_uncallable_after_qc=len(eligible) - n_call,
                conditional_callability=n_call / len(eligible),
                end_to_end_callability=n_call / len(linked)))
            for metric in METRICS:
                yes = [r[metric] for r in eligible if r["callable_binary"] == 1]
                no = [r[metric] for r in eligible if r["callable_binary"] == 0]
                row = dict(
                    dataset_id=cohort, antibiotic=drug, metric=metric,
                    n_callable=len(yes), n_uncallable=len(no),
                    callable_median=float(np.median(yes)) if yes else None,
                    uncallable_median=float(np.median(no)) if no else None,
                    mann_whitney_u=None, mann_whitney_p=None,
                    cliffs_delta_callable_minus_uncallable=None,
                    auc_metric_predicting_callable=None,
                    analysis_scope="final_genome_QC_PASS")
                if len(yes) >= 3 and len(no) >= 3:
                    test = mannwhitneyu(yes, no, alternative="two-sided", method="auto")
                    row["mann_whitney_u"] = float(test.statistic)
                    row["mann_whitney_p"] = float(test.pvalue)
                    row["cliffs_delta_callable_minus_uncallable"] = (
                        2 * float(test.statistic) / (len(yes) * len(no)) - 1)
                    row["auc_metric_predicting_callable"] = float(roc_auc_score(
                        [r["callable_binary"] for r in eligible],
                        [r[metric] for r in eligible]))
                associations.append(row)
    return accounting, associations


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path("MGEN_Supplementary_Data.xlsx"))
    parser.add_argument("--output", type=Path, default=Path("qc_reanalysis_output"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    samples = read_samples(args.input)
    accounting, associations = summarize(samples)
    write_csv(args.output / "Table_S2_QC_accounting.csv", accounting)
    write_csv(args.output / "Table_S6_QC_conditioned.csv", associations)
    write_csv(args.output / "Figure_S1_source_data.csv",
              [r for r in associations if r["cliffs_delta_callable_minus_uncallable"] is not None])
    environment = dict(
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
        scikit_learn=sklearn.__version__, openpyxl=openpyxl.__version__,
        input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        input_rows=len(samples), QC_pass_sample_drug_rows=sum(
            r["final_qc_status"] == "PASS" for r in samples),
        test="two-sided Mann-Whitney U, method=auto, continuity=True, minimum 3 per group",
        scope="final-genome-QC-pass rows only", randomness="none")
    (args.output / "analysis_environment.json").write_text(
        json.dumps(environment, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

