#!/usr/bin/env python3
"""Evaluate the frozen deterministic panels without refitting or marker selection."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import beta


BOOTSTRAP_REPLICATES = 2000
BOOTSTRAP_SEED = 20260830


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        if not rows:
            raise ValueError(f"cannot infer fields for empty output {path}")
        fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def exact_interval(successes: int, total: int, alpha: float = 0.05) -> tuple[float, float]:
    if total == 0:
        return math.nan, math.nan
    low = 0.0 if successes == 0 else float(beta.ppf(alpha / 2, successes, total - successes + 1))
    high = 1.0 if successes == total else float(beta.ppf(1 - alpha / 2, successes + 1, total - successes))
    return low, high


def confusion(y_true: list[int], y_pred: list[int]) -> tuple[int, int, int, int]:
    tp = sum(truth == 1 and pred == 1 for truth, pred in zip(y_true, y_pred))
    tn = sum(truth == 0 and pred == 0 for truth, pred in zip(y_true, y_pred))
    fp = sum(truth == 0 and pred == 1 for truth, pred in zip(y_true, y_pred))
    fn = sum(truth == 1 and pred == 0 for truth, pred in zip(y_true, y_pred))
    return tp, tn, fp, fn


def safe_ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else math.nan


def point_metrics(y_true: list[int], y_pred: list[int]) -> dict[str, float | int]:
    tp, tn, fp, fn = confusion(y_true, y_pred)
    sensitivity = safe_ratio(tp, tp + fn)
    specificity = safe_ratio(tn, tn + fp)
    balanced = (sensitivity + specificity) / 2 if not math.isnan(sensitivity) and not math.isnan(specificity) else math.nan
    has_two_classes = bool(tp + fn) and bool(tn + fp)
    mcc_denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / mcc_denominator if mcc_denominator else 0.0
    prevalence = safe_ratio(tp + fn, len(y_true))
    if has_two_classes:
        predicted_positive_precision = safe_ratio(tp, tp + fp)
        if math.isnan(predicted_positive_precision):
            predicted_positive_precision = 0.0
        predicted_positive_recall = tp / (tp + fn)
        average_precision = (
            predicted_positive_recall * predicted_positive_precision
            + (1 - predicted_positive_recall) * prevalence
        )
    else:
        average_precision = math.nan
    return {
        "n": len(y_true), "n_resistant": tp + fn, "n_susceptible": tn + fp,
        "tp": tp, "tn": tn, "fp": fp, "fn": fn,
        "prevalence": prevalence, "sensitivity": sensitivity, "specificity": specificity,
        "false_susceptible_rate": safe_ratio(fn, tp + fn), "false_resistant_rate": safe_ratio(fp, tn + fp),
        "balanced_accuracy": balanced, "mcc": mcc if has_two_classes else math.nan,
        "ppv": safe_ratio(tp, tp + fp), "npv": safe_ratio(tn, tn + fn),
        "auroc_binary_score": balanced if has_two_classes else math.nan,
        "auprc_binary_score": average_precision,
        "brier_binary_score": safe_ratio(fp + fn, len(y_true)),
    }


def bootstrap_intervals(rows: list[dict[str, str]], replicates: int = BOOTSTRAP_REPLICATES, seed: int = BOOTSTRAP_SEED) -> dict[str, tuple[float, float]]:
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row["bootstrap_unit"]].append(row)
    names = sorted(groups)
    rng = np.random.default_rng(seed)
    values: dict[str, list[float]] = defaultdict(list)
    for _ in range(replicates):
        sampled = rng.choice(names, size=len(names), replace=True)
        sample_rows = [row for name in sampled for row in groups[str(name)]]
        metrics = point_metrics(
            [1 if row["phenotype"] == "R" else 0 for row in sample_rows],
            [1 if row["prediction"] == "R" else 0 for row in sample_rows],
        )
        for metric in ("balanced_accuracy", "mcc", "auroc_binary_score", "auprc_binary_score", "brier_binary_score"):
            value = float(metrics[metric])
            if not math.isnan(value):
                values[metric].append(value)
    return {metric: (float(np.quantile(observed, 0.025)), float(np.quantile(observed, 0.975))) for metric, observed in values.items() if observed}


def summarize(group_rows: list[dict[str, str]], group_values: dict[str, str]) -> dict[str, object]:
    y_true = [1 if row["phenotype"] == "R" else 0 for row in group_rows]
    y_pred = [1 if row["prediction"] == "R" else 0 for row in group_rows]
    metrics = point_metrics(y_true, y_pred)
    result: dict[str, object] = {**group_values, **metrics}
    proportions = {
        "sensitivity": (int(metrics["tp"]), int(metrics["tp"]) + int(metrics["fn"])),
        "specificity": (int(metrics["tn"]), int(metrics["tn"]) + int(metrics["fp"])),
        "false_susceptible_rate": (int(metrics["fn"]), int(metrics["tp"]) + int(metrics["fn"])),
        "false_resistant_rate": (int(metrics["fp"]), int(metrics["tn"]) + int(metrics["fp"])),
        "ppv": (int(metrics["tp"]), int(metrics["tp"]) + int(metrics["fp"])),
        "npv": (int(metrics["tn"]), int(metrics["tn"]) + int(metrics["fn"])),
    }
    for metric, (successes, total) in proportions.items():
        result[f"{metric}_ci_low"], result[f"{metric}_ci_high"] = exact_interval(successes, total)
    for metric, (low, high) in bootstrap_intervals(group_rows).items():
        result[f"{metric}_bootstrap_low"] = low
        result[f"{metric}_bootstrap_high"] = high
    return result


def assembly_predictions(root: Path, include_zenodo: bool = False) -> list[dict[str, str]]:
    calls = read_csv(root / "results/panels/assembly_marker_calls.csv")
    final_qc_path = root / "results/qc/assembly_qc_with_checkm2.csv"
    qc_path = final_qc_path if final_qc_path.exists() else root / "results/qc/assembly_qc.csv"
    qc = {(row["dataset_id"], row["isolate_id"]): row for row in read_csv(qc_path)}
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in calls:
        grouped[(row["dataset_id"], row["isolate_id"])].append(row)
    predictions: list[dict[str, str]] = []
    for key, sample_calls in grouped.items():
        if key[0] == "ZENODO_10369064" and not include_zenodo:
            continue
        quality = qc.get(key, {})
        basic_qc = quality.get("final_qc_status", quality.get("basic_qc_status", "FAIL"))
        rrna = [row for row in sample_calls if row["gene"] == "23S_rRNA"]
        rrna_callable = any(row["status"] == "PASS" for row in rrna)
        rrna_markers = sorted({row["change"] for row in rrna if row["known_resistance_marker"] == "yes"})
        clr = "R" if basic_qc == "PASS" and rrna_markers else "S" if basic_qc == "PASS" and rrna_callable else "UNCALLABLE"
        predictions.append({
            "dataset_id": key[0], "isolate_id": key[1], "antibiotic": "clarithromycin", "prediction": clr,
            "sequence_support": "ASSEMBLY_ONLY", "marker_summary": ";".join(rrna_markers),
            "callable": "yes" if clr != "UNCALLABLE" else "no", "basic_qc_status": basic_qc,
        })
        gyr = [row for row in sample_calls if row["gene"] == "gyrA"]
        positions = {int(row["position"]): row for row in gyr if row["status"] == "PASS"}
        gyr_callable = all(position in positions for position in (87, 88, 91))
        gyr_markers = sorted({row["change"] for row in gyr if row["known_resistance_marker"] == "yes"})
        lvx = "R" if basic_qc == "PASS" and gyr_callable and gyr_markers else "S" if basic_qc == "PASS" and gyr_callable else "UNCALLABLE"
        predictions.append({
            "dataset_id": key[0], "isolate_id": key[1], "antibiotic": "levofloxacin", "prediction": lvx,
            "sequence_support": "ASSEMBLY_ONLY", "marker_summary": ";".join(gyr_markers),
            "callable": "yes" if lvx != "UNCALLABLE" else "no", "basic_qc_status": basic_qc,
        })
    return predictions


def raw_read_predictions(root: Path) -> list[dict[str, str]]:
    """Attach complete genome QC to authoritative raw-read marker calls."""
    raw_path = root / "results/panels/zenodo_read_marker_predictions.csv"
    if not raw_path.exists():
        return []
    final_qc_path = root / "results/qc/assembly_qc_with_checkm2.csv"
    qc_path = final_qc_path if final_qc_path.exists() else root / "results/qc/assembly_qc.csv"
    qc = {(row["dataset_id"], row["isolate_id"]): row for row in read_csv(qc_path)}
    predictions = []
    for original in read_csv(raw_path):
        row = dict(original)
        quality = qc.get((row["dataset_id"], row["isolate_id"]), {})
        quality_status = quality.get("final_qc_status", quality.get("basic_qc_status", "FAIL"))
        row["basic_qc_status"] = quality_status
        if quality_status != "PASS":
            # Keep the raw evidence in its source table, but prevent a genome
            # that failed the frozen sample-level gate from entering a
            # diagnostic denominator or a leakage model.
            row["prediction"] = "UNCALLABLE"
            row["callable"] = "no"
        predictions.append(row)
    return predictions


def merge_benchmark(root: Path) -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    predictions = assembly_predictions(root)
    predictions.extend(raw_read_predictions(root))
    phenotypes = {
        (row["dataset_id"], row["isolate_id"], row["antibiotic"]): row
        for row in read_csv(root / "metadata/phenotype_manifest.csv")
        if row["antibiotic"] in {"clarithromycin", "levofloxacin"}
    }
    clones_path = root / "results/qc/near_clone_groups.csv"
    clones = {(row["dataset_id"], row["isolate_id"]): row["near_clone_group"] for row in read_csv(clones_path)} if clones_path.exists() else {}
    lineages_path = root / "results/lineage_validation/lineage_assignments.csv"
    lineages = {(row["dataset_id"], row["isolate_id"]): row["lineage_recomputed"] for row in read_csv(lineages_path)} if lineages_path.exists() else {}
    sample_rows: list[dict[str, str]] = []
    for prediction in predictions:
        key = (prediction["dataset_id"], prediction["isolate_id"], prediction["antibiotic"])
        phenotype = phenotypes.get(key)
        if phenotype is None:
            continue
        label = phenotype["susceptibility_original"]
        row = {
            **prediction, "phenotype": label, "phenotype_quality": phenotype["phenotype_quality"],
            "phenotype_recomputed": phenotype["susceptibility_recomputed"],
            "mic_raw": phenotype["mic_raw"], "mic_numeric": phenotype["mic_numeric"], "mic_operator": phenotype["mic_operator"],
            "borderline_mic": phenotype["borderline_mic"],
            "ast_method": phenotype["ast_method"], "medium": phenotype["medium"],
            "breakpoint_standard": phenotype["breakpoint_standard"],
            "breakpoint_version": phenotype["breakpoint_version"],
            "near_clone_group": clones.get((prediction["dataset_id"], prediction["isolate_id"]), ""),
            "lineage_recomputed": lineages.get((prediction["dataset_id"], prediction["isolate_id"]), ""),
        }
        row["country"] = prediction["isolate_id"].split("-", 1)[0] if prediction["dataset_id"] == "HPGP_GLOBAL" else "China"
        row["bootstrap_unit"] = row["near_clone_group"] or f"{row['dataset_id']}::{row['isolate_id']}"
        row["analysis_status"] = "PRIMARY" if label in {"S", "R"} and prediction["prediction"] in {"S", "R"} else "EXCLUDED_AMBIGUOUS_OR_UNCALLABLE"
        row["correct"] = "yes" if row["analysis_status"] == "PRIMARY" and label == prediction["prediction"] else "no" if row["analysis_status"] == "PRIMARY" else ""
        sample_rows.append(row)

    metric_rows: list[dict[str, object]] = []
    eligible = [row for row in sample_rows if row["analysis_status"] == "PRIMARY"]
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in eligible:
        groups[(row["dataset_id"], row["antibiotic"])].append(row)
    for (dataset, antibiotic), rows in sorted(groups.items()):
        metric_rows.append(summarize(rows, {"stratum_type": "COHORT", "stratum": dataset, "dataset_id": dataset, "antibiotic": antibiotic}))
    country_groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in eligible:
        country_groups[(row["country"], row["antibiotic"])].append(row)
    for (country, antibiotic), rows in sorted(country_groups.items()):
        if any(row["phenotype"] == "S" for row in rows) and any(row["phenotype"] == "R" for row in rows):
            metric_rows.append(summarize(rows, {"stratum_type": "COUNTRY", "stratum": country, "dataset_id": "MULTI_COHORT", "antibiotic": antibiotic}))
    lineage_groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in eligible:
        if row["lineage_recomputed"]:
            lineage_groups[(row["lineage_recomputed"], row["antibiotic"])].append(row)
    for (lineage, antibiotic), rows in sorted(lineage_groups.items()):
        n_s = sum(row["phenotype"] == "S" for row in rows)
        n_r = sum(row["phenotype"] == "R" for row in rows)
        if len(rows) >= 30 and n_s >= 10 and n_r >= 10:
            metric_rows.append(summarize(rows, {"stratum_type": "LINEAGE", "stratum": lineage, "dataset_id": "MULTI_COHORT", "antibiotic": antibiotic}))
    return sample_rows, metric_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.root.resolve()
    sample_rows, metric_rows = merge_benchmark(root)
    write_csv(root / "results/external_validation/sample_level_predictions.csv", sample_rows)
    write_csv(root / "results/external_validation/frozen_panel_metrics.csv", metric_rows)


if __name__ == "__main__":
    main()
