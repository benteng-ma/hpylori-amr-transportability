#!/usr/bin/env python3
"""Apply the frozen transportability gates to completed panel metrics."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


EXTERNAL = ("CHINA_NINGXIA_2022", "ZENODO_10369064")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def passed(row: dict[str, str]) -> bool:
    return float(row["sensitivity"]) >= 0.90 and float(row["specificity"]) >= 0.90 and float(row["false_susceptible_rate"]) <= 0.10


def classify(drug: str, cohort_rows: list[dict[str, str]], lineage_rows: list[dict[str, str]]) -> dict[str, object]:
    by_cohort = {row["dataset_id"]: row for row in cohort_rows if row["antibiotic"] == drug and row["dataset_id"] in EXTERNAL}
    if set(by_cohort) != set(EXTERNAL) or any(int(row["n_resistant"]) < 10 or int(row["n_susceptible"]) < 10 for row in by_cohort.values()):
        label = "INSUFFICIENT_DATA"
        reason = "fewer than two estimable external cohorts with at least 10 S and 10 R"
    elif any(float(row["false_susceptible_rate"]) > 0.10 for row in by_cohort.values()):
        label = "HIGH_FALSE_SUSCEPTIBLE_RISK"
        reason = "at least one eligible external cohort exceeded the frozen 0.10 false-susceptible gate"
    else:
        metric_differences = {
            metric: abs(float(by_cohort[EXTERNAL[0]][metric]) - float(by_cohort[EXTERNAL[1]][metric]))
            for metric in ("sensitivity", "specificity")
        }
        eligible_lineages = [row for row in lineage_rows if row["antibiotic"] == drug]
        lineage_failure = any(not passed(row) for row in eligible_lineages)
        if all(passed(row) for row in by_cohort.values()) and max(metric_differences.values()) <= 0.10 and not lineage_failure:
            label = "ROBUSTLY_TRANSPORTABLE"
            reason = "all eligible external cohorts and lineages passed frozen performance and heterogeneity gates"
        elif lineage_failure:
            label = "LINEAGE_DEPENDENT"
            reason = "at least one prespecified estimable core-SNP cluster failed the frozen gate"
        elif passed(by_cohort[EXTERNAL[0]]) != passed(by_cohort[EXTERNAL[1]]):
            label = "PHENOTYPE_METHOD_SENSITIVE"
            reason = "performance gate differed between MIC-linked and custom disk-diffusion cohorts"
        elif max(metric_differences.values()) > 0.10:
            label = "STUDY_DEPENDENT"
            reason = "between-cohort sensitivity or specificity differed by more than 0.10"
        else:
            label = "NON_TRANSPORTABLE"
            reason = "external performance failed without an attributable eligible lineage or method contrast"
    return {
        "antibiotic": drug, "transportability_label": label, "reason": reason,
        "external_cohorts_required": ";".join(EXTERNAL), "geographic_scope": "two Chinese external cohorts; no global claim",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.root.resolve()
    metrics = read_csv(root / "results/external_validation/frozen_panel_metrics.csv")
    cohorts = [row for row in metrics if row["stratum_type"] == "COHORT"]
    lineages = [row for row in metrics if row["stratum_type"] == "LINEAGE"]
    rows = [classify(drug, cohorts, lineages) for drug in ("clarithromycin", "levofloxacin")]
    output = root / "results/external_validation/transportability_classification.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    (root / "results/external_validation/transportability_classification.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
