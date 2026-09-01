#!/usr/bin/env python3
"""Prespecified phenotype, callability, and cohort-limited MIC sensitivities."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

from benchmark_frozen_panels import summarize, write_csv


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def scenario_label(row: dict[str, str], scenario: str) -> str:
    original = row["phenotype"]
    if scenario == "PRIMARY_ORIGINAL_I_EXCLUDED":
        return original if original in {"S", "R"} else ""
    if scenario == "EXCLUDE_BORDERLINE_MIC":
        return "" if row["borderline_mic"] == "yes" else original if original in {"S", "R"} else ""
    if scenario == "I_AS_S":
        return "S" if original == "I" else original if original in {"S", "R"} else ""
    if scenario == "I_AS_R":
        return "R" if original == "I" else original if original in {"S", "R"} else ""
    recomputed = row.get("phenotype_recomputed", "")
    if scenario == "RECOMPUTED_I_EXCLUDED":
        return recomputed if recomputed in {"S", "R"} else ""
    if scenario == "RECOMPUTED_I_AS_S":
        return "S" if recomputed == "I" else recomputed if recomputed in {"S", "R"} else ""
    if scenario == "RECOMPUTED_I_AS_R":
        return "R" if recomputed == "I" else recomputed if recomputed in {"S", "R"} else ""
    raise ValueError(scenario)


def callability_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["dataset_id"], row["antibiotic"])].append(row)
    output = []
    for (dataset, antibiotic), values in sorted(grouped.items()):
        callable_n = sum(row["prediction"] in {"S", "R"} for row in values)
        output.append({
            "dataset_id": dataset, "antibiotic": antibiotic, "n_phenotype_linked": len(values),
            "n_callable": callable_n, "n_uncallable": len(values) - callable_n,
            "callability": callable_n / len(values) if values else math.nan,
            "sequence_support": ";".join(sorted({row["sequence_support"] for row in values})),
        })
    return output


def sensitivity_metrics(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    scenarios = [
        "PRIMARY_ORIGINAL_I_EXCLUDED", "EXCLUDE_BORDERLINE_MIC", "I_AS_S", "I_AS_R",
        "RECOMPUTED_I_EXCLUDED", "RECOMPUTED_I_AS_S", "RECOMPUTED_I_AS_R",
    ]
    output: list[dict[str, object]] = []
    for scenario in scenarios:
        grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for original in rows:
            label = scenario_label(original, scenario)
            if not label or original["prediction"] not in {"S", "R"}:
                continue
            row = dict(original)
            row["phenotype"] = label
            grouped[(row["dataset_id"], row["antibiotic"])].append(row)
        for (dataset, antibiotic), values in sorted(grouped.items()):
            if not any(row["phenotype"] == "S" for row in values) or not any(row["phenotype"] == "R" for row in values):
                continue
            output.append(summarize(values, {
                "scenario": scenario, "stratum_type": "COHORT", "stratum": dataset,
                "dataset_id": dataset, "antibiotic": antibiotic,
            }))
    primary = [row for row in rows if row["phenotype"] in {"S", "R"} and row["prediction"] in {"S", "R"}]
    for column, stratum_type in (
        ("ast_method", "AST_METHOD"), ("sequence_support", "SEQUENCE_SUPPORT"),
        ("country", "COUNTRY"), ("lineage_recomputed", "LINEAGE_DESCRIPTIVE"),
    ):
        grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for row in primary:
            if row.get(column, ""):
                grouped[(row[column], row["antibiotic"])].append(row)
        for (stratum, antibiotic), values in sorted(grouped.items()):
            if any(row["phenotype"] == "S" for row in values) and any(row["phenotype"] == "R" for row in values):
                output.append(summarize(values, {
                    "scenario": "PRIMARY_ORIGINAL_I_EXCLUDED", "stratum_type": stratum_type,
                    "stratum": stratum, "dataset_id": "MULTI_COHORT", "antibiotic": antibiotic,
                }))
    high_quality = [row for row in primary if row.get("basic_qc_status") == "PASS"]
    for antibiotic in ("clarithromycin", "levofloxacin"):
        values = [row for row in high_quality if row["antibiotic"] == antibiotic]
        if values and any(row["phenotype"] == "S" for row in values) and any(row["phenotype"] == "R" for row in values):
            output.append(summarize(values, {
                "scenario": "PRIMARY_ORIGINAL_I_EXCLUDED", "stratum_type": "HIGH_QUALITY_GENOMES_ONLY",
                "stratum": "FROZEN_COMPLETE_GENOME_QC_PASS", "dataset_id": "MULTI_COHORT", "antibiotic": antibiotic,
            }))
    return output


def mic_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    eligible = [
        row for row in rows
        if row["dataset_id"] == "CHINA_NINGXIA_2022" and row["mic_numeric"] and row["prediction"] in {"S", "R"}
    ]
    output: list[dict[str, object]] = []
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    censoring: dict[tuple[str, str], int] = defaultdict(int)
    for row in eligible:
        key = (row["antibiotic"], row["prediction"])
        grouped[key].append(float(row["mic_numeric"]))
        censoring[key] += row["mic_operator"] != "="
    for (antibiotic, prediction), values in sorted(grouped.items()):
        array = np.asarray(values)
        output.append({
            "dataset_id": "CHINA_NINGXIA_2022", "antibiotic": antibiotic, "marker_prediction": prediction,
            "n": len(values), "n_right_censored": censoring[(antibiotic, prediction)],
            "mic_median_lower_bound": float(np.median(array)), "mic_q1_lower_bound": float(np.quantile(array, 0.25)),
            "mic_q3_lower_bound": float(np.quantile(array, 0.75)),
            "mann_whitney_u": "", "mann_whitney_p": "",
            "interpretation": "cohort-limited lower-bound rank summary; not a multi-cohort continuous endpoint",
        })
    for antibiotic in sorted({row["antibiotic"] for row in eligible}):
        susceptible = grouped.get((antibiotic, "S"), [])
        resistant = grouped.get((antibiotic, "R"), [])
        if susceptible and resistant:
            test = mannwhitneyu(susceptible, resistant, alternative="two-sided", method="auto")
            output.append({
                "dataset_id": "CHINA_NINGXIA_2022", "antibiotic": antibiotic, "marker_prediction": "R_VS_S_TEST",
                "n": len(susceptible) + len(resistant), "n_right_censored": censoring.get((antibiotic, "S"), 0) + censoring.get((antibiotic, "R"), 0),
                "mic_median_lower_bound": "", "mic_q1_lower_bound": "", "mic_q3_lower_bound": "",
                "mann_whitney_u": float(test.statistic), "mann_whitney_p": float(test.pvalue),
                "interpretation": "exploratory lower-bound rank test with censoring retained in source data",
            })
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.root.resolve()
    rows = read_csv(root / "results/external_validation/sample_level_predictions.csv")
    write_csv(root / "results/qc/panel_callability.csv", callability_rows(rows))
    write_csv(root / "results/sensitivity/phenotype_sensitivity_metrics.csv", sensitivity_metrics(rows))
    write_csv(root / "results/mic/ningxia_mic_lower_bound_summary.csv", mic_rows(rows))
    boundaries = [
        {"antibiotic": drug, "result_class": "INSUFFICIENT_DATA", "reason": "no exactly reproducible frozen published genomic panel across eligible cohorts"}
        for drug in ("amoxicillin", "furazolidone", "metronidazole", "rifampin", "tetracycline")
    ]
    write_csv(root / "results/external_validation/secondary_drug_boundaries.csv", boundaries)


if __name__ == "__main__":
    main()
