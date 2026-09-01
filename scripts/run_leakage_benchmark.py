#!/usr/bin/env python3
"""Compare random and independence-aware validation without feature selection."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from benchmark_frozen_panels import point_metrics


SEED = 20260830
MARKERS = {
    "clarithromycin": ["CLR_A2142_ANY", "CLR_A2143G"],
    "levofloxacin": ["LVX_A88V", "LVX_A88P", "LVX_N87K", "LVX_N87I", "LVX_D91G", "LVX_D91N", "LVX_D91Y"],
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def assembly_features(root: Path) -> dict[tuple[str, str, str], dict[str, int]]:
    calls = read_csv(root / "results/panels/assembly_marker_calls.csv")
    output: dict[tuple[str, str, str], dict[str, int]] = {}
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in calls:
        grouped[(row["dataset_id"], row["isolate_id"])].append(row)
    for (dataset, isolate), values in grouped.items():
        clr_changes = {row["change"] for row in values if row["gene"] == "23S_rRNA" and row["known_resistance_marker"] == "yes"}
        output[(dataset, isolate, "clarithromycin")] = {
            "CLR_A2142_ANY": int(bool(clr_changes & {"A2142G", "A2142C"})),
            "CLR_A2143G": int("A2143G" in clr_changes),
        }
        gyr_changes = {row["change"] for row in values if row["gene"] == "gyrA" and row["known_resistance_marker"] == "yes"}
        output[(dataset, isolate, "levofloxacin")] = {f"LVX_{change}": int(change in gyr_changes) for change in ("A88V", "A88P", "N87K", "N87I", "D91G", "D91N", "D91Y")}
    return output


def raw_features(row: dict[str, str]) -> dict[str, int]:
    if row["antibiotic"] == "clarithromycin":
        fractions = {}
        for token in row["marker_summary"].split(";"):
            if ":" in token:
                label, value = token.split(":", 1)
                fractions[label] = float(value)
        return {"CLR_A2142_ANY": int(fractions.get("A2142", 0) >= 0.20), "CLR_A2143G": int(fractions.get("A2143", 0) >= 0.20)}
    changes = set(row["marker_summary"].split(";"))
    return {f"LVX_{change}": int(change in changes) for change in ("A88V", "A88P", "N87K", "N87I", "D91G", "D91N", "D91Y")}


def feature_table(root: Path) -> pd.DataFrame:
    sample_rows = read_csv(root / "results/external_validation/sample_level_predictions.csv")
    assembly = assembly_features(root)
    rows = []
    for sample in sample_rows:
        if sample["analysis_status"] != "PRIMARY" or not sample["lineage_recomputed"]:
            continue
        key = (sample["dataset_id"], sample["isolate_id"], sample["antibiotic"])
        features = assembly.get(key) if sample["dataset_id"] != "ZENODO_10369064" else raw_features(sample)
        if features is None:
            continue
        row = {
            **sample, **features, "outcome": int(sample["phenotype"] == "R"),
            "country": sample["isolate_id"].split("-", 1)[0] if sample["dataset_id"] == "HPGP_GLOBAL" else sample["dataset_id"],
        }
        rows.append(row)
    return pd.DataFrame(rows)


def model_pipeline(drug: str, model_name: str) -> Pipeline:
    marker_columns = MARKERS[drug]
    if model_name == "MUTATION_ONLY_LOGISTIC":
        transformer = ColumnTransformer([("markers", "passthrough", marker_columns)], remainder="drop")
    elif model_name == "LINEAGE_ONLY_LOGISTIC":
        transformer = ColumnTransformer([("lineage", OneHotEncoder(handle_unknown="ignore"), ["lineage_recomputed"])], remainder="drop")
    elif model_name == "MUTATION_PLUS_LINEAGE_LOGISTIC":
        transformer = ColumnTransformer([
            ("markers", "passthrough", marker_columns),
            ("lineage", OneHotEncoder(handle_unknown="ignore"), ["lineage_recomputed"]),
        ], remainder="drop")
    else:
        raise ValueError(model_name)
    return Pipeline([
        ("features", transformer),
        ("model", LogisticRegression(C=1.0, solver="liblinear", random_state=SEED, max_iter=1000)),
    ])


def probability_metrics(y_true: list[int], probabilities: list[float]) -> dict[str, float]:
    truth = np.asarray(y_true, dtype=int)
    probability = np.clip(np.asarray(probabilities, dtype=float), 1e-8, 1 - 1e-8)
    result = {
        "auroc_probability": float(roc_auc_score(truth, probability)) if len(set(truth)) == 2 else np.nan,
        "auprc_probability": float(average_precision_score(truth, probability)) if len(set(truth)) == 2 else np.nan,
        "brier_probability": float(brier_score_loss(truth, probability)),
        "calibration_intercept": np.nan,
        "calibration_slope": np.nan,
    }
    if len(set(truth)) == 2 and float(np.ptp(probability)) > 0:
        logits = np.log(probability / (1 - probability)).reshape(-1, 1)
        try:
            # A vanishingly small L2 penalty keeps separated small validation
            # folds finite while remaining an effectively unregularized
            # logistic recalibration fit.
            calibration = LogisticRegression(C=1e6, solver="lbfgs", max_iter=500)
            calibration.fit(logits, truth)
            result["calibration_intercept"] = float(calibration.intercept_[0])
            result["calibration_slope"] = float(calibration.coef_[0, 0])
        except (ValueError, FloatingPointError):
            pass
    return result


def evaluate(train: pd.DataFrame, test: pd.DataFrame, drug: str, model_name: str, split_type: str, fold: str) -> dict[str, object]:
    pipeline = model_pipeline(drug, model_name)
    pipeline.fit(train, train["outcome"])
    prediction = pipeline.predict(test).astype(int).tolist()
    probabilities = pipeline.predict_proba(test)[:, 1].astype(float).tolist()
    metrics = point_metrics(test["outcome"].astype(int).tolist(), prediction)
    probability_summary = probability_metrics(test["outcome"].astype(int).tolist(), probabilities)
    return {
        "antibiotic": drug, "model": model_name, "split_type": split_type, "fold": fold,
        "n_train": len(train), "n_test": len(test), **metrics, **probability_summary,
    }


def eligible_group(values: pd.DataFrame) -> bool:
    return len(values) >= 30 and int(values["outcome"].sum()) >= 10 and int((1 - values["outcome"]).sum()) >= 10


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.root.resolve()
    table = feature_table(root)
    output_rows: list[dict[str, object]] = []
    models = ("MUTATION_ONLY_LOGISTIC", "LINEAGE_ONLY_LOGISTIC", "MUTATION_PLUS_LINEAGE_LOGISTIC")
    for drug in ("clarithromycin", "levofloxacin"):
        development = table[(table["dataset_id"] == "HPGP_GLOBAL") & (table["antibiotic"] == drug)].reset_index(drop=True)
        splitter = StratifiedShuffleSplit(n_splits=100, test_size=0.20, random_state=SEED)
        for repeat, (train_index, test_index) in enumerate(splitter.split(development, development["outcome"]), start=1):
            for model_name in models:
                output_rows.append(evaluate(development.iloc[train_index], development.iloc[test_index], drug, model_name, "RANDOM_ISOLATE_SPLIT", str(repeat)))

        groups = development["near_clone_group"].where(development["near_clone_group"].astype(bool), development["isolate_id"])
        grouped = GroupKFold(n_splits=5)
        for fold, (train_index, test_index) in enumerate(grouped.split(development, development["outcome"], groups), start=1):
            for model_name in models:
                output_rows.append(evaluate(development.iloc[train_index], development.iloc[test_index], drug, model_name, "CLONE_GROUPED_SPLIT", str(fold)))

        for group_column, split_type in (("country", "LEAVE_COUNTRY_OUT"), ("lineage_recomputed", "LEAVE_LINEAGE_OUT")):
            for group_name, test in development.groupby(group_column):
                if not eligible_group(test):
                    continue
                train = development[development[group_column] != group_name]
                if train["outcome"].nunique() < 2:
                    continue
                for model_name in models:
                    output_rows.append(evaluate(train, test, drug, model_name, split_type, str(group_name)))

        for external_name in ("CHINA_NINGXIA_2022", "ZENODO_10369064"):
            test = table[(table["dataset_id"] == external_name) & (table["antibiotic"] == drug)]
            if not eligible_group(test):
                continue
            for model_name in models:
                output_rows.append(evaluate(development, test, drug, model_name, "HPGP_TO_EXTERNAL", external_name))

    output = root / "results/lineage_validation/leakage_benchmark_folds.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0].keys()))
        writer.writeheader()
        writer.writerows(output_rows)

    summary_rows = []
    frame = pd.DataFrame(output_rows)
    metric_columns = [
        "sensitivity", "specificity", "false_susceptible_rate", "balanced_accuracy", "mcc",
        "auroc_probability", "auprc_probability", "brier_probability",
        "calibration_intercept", "calibration_slope",
    ]
    for keys, group in frame.groupby(["antibiotic", "model", "split_type"]):
        row = {"antibiotic": keys[0], "model": keys[1], "split_type": keys[2], "folds": len(group)}
        for metric in metric_columns:
            row[f"{metric}_mean"] = group[metric].mean()
            row[f"{metric}_sd"] = group[metric].std(ddof=1) if len(group) > 1 else np.nan
        summary_rows.append(row)
    summary = root / "results/lineage_validation/leakage_benchmark_summary.csv"
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)


if __name__ == "__main__":
    main()
