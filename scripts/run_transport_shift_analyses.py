#!/usr/bin/env python3
"""Post-freeze transport-shift analyses for the H. pylori AMR benchmark.

These analyses explain, but never refit, the frozen resistance catalogues. They
separate three failure layers: analytic availability, population/case-mix
shift, and conditional genotype-phenotype transport within a shared lineage.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, fisher_exact
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import auc, roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


HPGP = "HPGP_GLOBAL"
NINGXIA = "CHINA_NINGXIA_2022"
READ_COHORT = "ZENODO_10369064"
DOMINANT_LINEAGE = "SNP_CLUSTER_03"
SEED = 20260830


def exact_interval(successes: int, total: int, alpha: float = 0.05) -> tuple[float, float]:
    """Clopper-Pearson interval with defined 0/1 boundaries."""
    if total <= 0:
        return np.nan, np.nan
    low = 0.0 if successes == 0 else float(beta.ppf(alpha / 2, successes, total - successes + 1))
    high = 1.0 if successes == total else float(beta.ppf(1 - alpha / 2, successes + 1, total - successes))
    return low, high


def add_binary_columns(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["phenotype_binary"] = (frame["phenotype"] == "R").astype(int)
    frame["marker_binary"] = (frame["prediction"] == "R").astype(int)
    frame["outcome"] = np.select(
        [
            (frame["phenotype"] == "R") & (frame["prediction"] == "R"),
            (frame["phenotype"] == "R") & (frame["prediction"] == "S"),
            (frame["phenotype"] == "S") & (frame["prediction"] == "R"),
        ],
        ["TP", "FN", "FP"],
        default="TN",
    )
    return frame


def primary_levofloxacin(predictions: pd.DataFrame) -> pd.DataFrame:
    keep = (
        (predictions["antibiotic"] == "levofloxacin")
        & (predictions["analysis_status"] == "PRIMARY")
        & predictions["phenotype"].isin(["R", "S"])
        & predictions["prediction"].isin(["R", "S"])
    )
    return add_binary_columns(predictions.loc[keep])


def metric_record(group: pd.DataFrame, metric: str) -> dict[str, float | int]:
    if metric == "sensitivity":
        subset = group[group["phenotype"] == "R"]
        success = int((subset["prediction"] == "R").sum())
    elif metric == "specificity":
        subset = group[group["phenotype"] == "S"]
        success = int((subset["prediction"] == "S").sum())
    elif metric == "false_susceptible_rate":
        subset = group[group["phenotype"] == "R"]
        success = int((subset["prediction"] == "S").sum())
    elif metric == "marker_negative_resistance":
        subset = group[group["prediction"] == "S"]
        success = int((subset["phenotype"] == "R").sum())
    elif metric == "marker_positive_resistance":
        subset = group[group["prediction"] == "R"]
        success = int((subset["phenotype"] == "R").sum())
    else:
        raise ValueError(metric)
    total = int(len(subset))
    estimate = success / total if total else np.nan
    low, high = exact_interval(success, total)
    return {
        "metric": metric,
        "successes": success,
        "total": total,
        "estimate": estimate,
        "ci_low": low,
        "ci_high": high,
    }


def lineage_performance(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    metrics = [
        "sensitivity",
        "specificity",
        "false_susceptible_rate",
        "marker_negative_resistance",
        "marker_positive_resistance",
    ]
    for (dataset, lineage), group in frame.groupby(["dataset_id", "lineage_recomputed"], dropna=False):
        for metric in metrics:
            rows.append({"dataset_id": dataset, "lineage_recomputed": lineage, **metric_record(group, metric)})
    for dataset, group in frame.groupby("dataset_id"):
        for metric in metrics:
            rows.append({"dataset_id": dataset, "lineage_recomputed": "ALL", **metric_record(group, metric)})
    return pd.DataFrame(rows)


def comparison_table(frame: pd.DataFrame, metric: str) -> tuple[int, int]:
    record = metric_record(frame, metric)
    return int(record["successes"]), int(record["total"] - record["successes"])


def dominant_lineage_comparisons(frame: pd.DataFrame) -> pd.DataFrame:
    dominant = frame[frame["lineage_recomputed"] == DOMINANT_LINEAGE]
    rows: list[dict] = []
    for comparator in [HPGP, READ_COHORT]:
        for metric in ["sensitivity", "specificity", "marker_negative_resistance"]:
            a = comparison_table(dominant[dominant["dataset_id"] == NINGXIA], metric)
            b = comparison_table(dominant[dominant["dataset_id"] == comparator], metric)
            odds_ratio, p_value = fisher_exact([a, b], alternative="two-sided")
            rows.append(
                {
                    "lineage_recomputed": DOMINANT_LINEAGE,
                    "metric": metric,
                    "dataset_a": NINGXIA,
                    "dataset_b": comparator,
                    "dataset_a_successes": a[0],
                    "dataset_a_failures": a[1],
                    "dataset_b_successes": b[0],
                    "dataset_b_failures": b[1],
                    "odds_ratio": odds_ratio,
                    "fisher_p": p_value,
                }
            )
    return pd.DataFrame(rows)


def standardized_rate(source: pd.DataFrame, target: pd.DataFrame, metric: str) -> tuple[float, list[dict]]:
    phenotype = "R" if metric == "sensitivity" else "S"
    target_subset = target[target["phenotype"] == phenotype]
    source_subset = source[source["phenotype"] == phenotype]
    target_counts = target_subset["lineage_recomputed"].value_counts()
    total = int(target_counts.sum())
    pieces: list[dict] = []
    estimate = 0.0
    for lineage, count in target_counts.items():
        src = source_subset[source_subset["lineage_recomputed"] == lineage]
        if src.empty:
            return np.nan, []
        correct = (src["prediction"] == phenotype).mean()
        weight = count / total
        estimate += weight * correct
        pieces.append(
            {
                "lineage_recomputed": lineage,
                "target_count": int(count),
                "target_weight": weight,
                "source_count": int(len(src)),
                "source_rate": float(correct),
            }
        )
    return float(estimate), pieces


def lineage_standardization(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = frame[frame["dataset_id"] == HPGP]
    rows: list[dict] = []
    weights: list[dict] = []
    for target_id in [NINGXIA, READ_COHORT]:
        target = frame[frame["dataset_id"] == target_id]
        for metric in ["sensitivity", "specificity"]:
            observed = metric_record(target, metric)
            source_crude = metric_record(source, metric)
            standardized, pieces = standardized_rate(source, target, metric)
            rows.extend(
                [
                    {
                        "target_dataset": target_id,
                        "metric": metric,
                        "estimate_type": "target_observed",
                        "estimate": observed["estimate"],
                        "ci_low": observed["ci_low"],
                        "ci_high": observed["ci_high"],
                        "numerator": observed["successes"],
                        "denominator": observed["total"],
                    },
                    {
                        "target_dataset": target_id,
                        "metric": metric,
                        "estimate_type": "HpGP_crude",
                        "estimate": source_crude["estimate"],
                        "ci_low": source_crude["ci_low"],
                        "ci_high": source_crude["ci_high"],
                        "numerator": source_crude["successes"],
                        "denominator": source_crude["total"],
                    },
                    {
                        "target_dataset": target_id,
                        "metric": metric,
                        "estimate_type": "HpGP_standardized_to_target_lineages",
                        "estimate": standardized,
                        "ci_low": np.nan,
                        "ci_high": np.nan,
                        "numerator": np.nan,
                        "denominator": np.nan,
                    },
                ]
            )
            for piece in pieces:
                weights.append({"target_dataset": target_id, "metric": metric, **piece})
    return pd.DataFrame(rows), pd.DataFrame(weights)


def bootstrap_auc_by_group(frame: pd.DataFrame, repeats: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    grouped = {label: list(group["cv_group"].unique()) for label, group in frame.groupby("cohort_binary")}
    by_group = {key: value.index.to_numpy() for key, value in frame.groupby("cv_group")}
    values: list[float] = []
    for _ in range(repeats):
        indices: list[int] = []
        for label in [0, 1]:
            choices = rng.choice(grouped[label], size=len(grouped[label]), replace=True)
            for choice in choices:
                indices.extend(by_group[choice])
        sample = frame.loc[indices]
        if sample["cohort_binary"].nunique() == 2:
            values.append(roc_auc_score(sample["cohort_binary"], sample["oof_probability"]))
    return tuple(np.quantile(values, [0.025, 0.975]))


def adversarial_validation(
    lineage: pd.DataFrame,
    group_map: pd.DataFrame,
    permutations: int,
    bootstrap_repeats: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pc_cols = [f"PC{i}" for i in range(1, 11)]
    merged = lineage.merge(group_map, on=["dataset_id", "isolate_id"], how="left", validate="one_to_one")
    merged["cv_group"] = merged["near_clone_group"].fillna(
        merged["dataset_id"].astype(str) + "|" + merged["isolate_id"].astype(str)
    )
    prediction_rows: list[pd.DataFrame] = []
    summary_rows: list[dict] = []
    null_rows: list[dict] = []
    coefficient_rows: list[dict] = []
    for pair_index, external in enumerate([NINGXIA, READ_COHORT]):
        data = merged[merged["dataset_id"].isin([HPGP, external])].copy().reset_index(drop=True)
        data["cohort_binary"] = (data["dataset_id"] == external).astype(int)
        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(C=1.0, class_weight="balanced", solver="liblinear", max_iter=2000),
        )
        cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED + pair_index)
        probabilities = np.full(len(data), np.nan)
        fold_ids = np.full(len(data), -1)
        fold_aucs: list[float] = []
        splits = list(cv.split(data[pc_cols], data["cohort_binary"], groups=data["cv_group"]))
        for fold, (train, test) in enumerate(splits):
            model.fit(data.loc[train, pc_cols], data.loc[train, "cohort_binary"])
            probabilities[test] = model.predict_proba(data.loc[test, pc_cols])[:, 1]
            fold_ids[test] = fold
            fold_aucs.append(roc_auc_score(data.loc[test, "cohort_binary"], probabilities[test]))
        data["oof_probability"] = probabilities
        data["fold"] = fold_ids
        data["comparison"] = f"{HPGP}_vs_{external}"
        pooled_auc = roc_auc_score(data["cohort_binary"], probabilities)
        low, high = bootstrap_auc_by_group(data, bootstrap_repeats, SEED + 100 + pair_index)

        rng = np.random.default_rng(SEED + 1000 + pair_index)
        null_scores: list[float] = []
        observed_labels = data["cohort_binary"].to_numpy()
        for repeat in range(permutations):
            permuted = rng.permutation(observed_labels)
            fold_scores: list[float] = []
            for train, test in splits:
                if len(np.unique(permuted[train])) < 2 or len(np.unique(permuted[test])) < 2:
                    continue
                model.fit(data.loc[train, pc_cols], permuted[train])
                prob = model.predict_proba(data.loc[test, pc_cols])[:, 1]
                fold_scores.append(roc_auc_score(permuted[test], prob))
            null_scores.append(float(np.mean(fold_scores)) if fold_scores else np.nan)
        null = np.asarray(null_scores, dtype=float)
        observed_mean_fold = float(np.mean(fold_aucs))
        permutation_p = (1 + int(np.nansum(null >= observed_mean_fold))) / (1 + int(np.isfinite(null).sum()))
        summary_rows.append(
            {
                "comparison": f"{HPGP}_vs_{external}",
                "external_dataset": external,
                "n_hpGP": int((data["dataset_id"] == HPGP).sum()),
                "n_external": int((data["dataset_id"] == external).sum()),
                "features": "PC1-PC10",
                "cv": "5-fold stratified near-clone-grouped",
                "pooled_oof_auc": pooled_auc,
                "group_bootstrap_ci_low": low,
                "group_bootstrap_ci_high": high,
                "mean_fold_auc": observed_mean_fold,
                "min_fold_auc": float(np.min(fold_aucs)),
                "max_fold_auc": float(np.max(fold_aucs)),
                "permutations": permutations,
                "permutation_p": permutation_p,
            }
        )
        for repeat, score in enumerate(null):
            null_rows.append(
                {
                    "comparison": f"{HPGP}_vs_{external}",
                    "repeat": repeat,
                    "null_mean_fold_auc": score,
                    "observed_mean_fold_auc": observed_mean_fold,
                }
            )
        model.fit(data[pc_cols], data["cohort_binary"])
        coefs = model.named_steps["logisticregression"].coef_[0]
        for feature, coefficient in zip(pc_cols, coefs):
            coefficient_rows.append(
                {
                    "comparison": f"{HPGP}_vs_{external}",
                    "feature": feature,
                    "standardized_log_odds_coefficient": coefficient,
                }
            )
        prediction_rows.append(
            data[
                [
                    "comparison",
                    "dataset_id",
                    "isolate_id",
                    "lineage_recomputed",
                    "cv_group",
                    "cohort_binary",
                    "fold",
                    "oof_probability",
                    *pc_cols,
                ]
            ]
        )
    return (
        pd.concat(prediction_rows, ignore_index=True),
        pd.DataFrame(summary_rows),
        pd.DataFrame(null_rows),
        pd.DataFrame(coefficient_rows),
    )


def mic_severity(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    ningxia = frame[frame["dataset_id"] == NINGXIA].copy()
    resistant = ningxia[ningxia["phenotype"] == "R"]
    rows: list[dict] = []
    summaries: list[dict] = []
    for outcome in ["TP", "FN"]:
        group = resistant[resistant["outcome"] == outcome]
        values = pd.to_numeric(group["mic_numeric"], errors="coerce").dropna()
        summaries.append(
            {
                "outcome": outcome,
                "n": int(len(group)),
                "median_lower_bound_mic_mg_L": float(values.median()),
                "q1_lower_bound_mic_mg_L": float(values.quantile(0.25)),
                "q3_lower_bound_mic_mg_L": float(values.quantile(0.75)),
                "right_censored_n": int((group["mic_operator"] == ">").sum()),
            }
        )
        for threshold in [2, 4, 8, 16, 32]:
            count = int((values >= threshold).sum())
            rows.append(
                {
                    "outcome": outcome,
                    "threshold_mg_L": threshold,
                    "n_at_or_above": count,
                    "denominator": int(len(values)),
                    "proportion_at_or_above": count / len(values) if len(values) else np.nan,
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(summaries)


def failure_atlas(
    frame: pd.DataFrame,
    manifold: pd.DataFrame,
    audit: pd.DataFrame,
) -> pd.DataFrame:
    ningxia = frame[frame["dataset_id"] == NINGXIA].copy()
    manifold_cols = [
        "dataset_id",
        "isolate_id",
        "PC1",
        "PC2",
        "development_5nn_distance",
        "development_5nn_percentile",
    ]
    audit_cols = [
        "dataset_id",
        "isolate_id",
        "full_gyrA_missense_variants",
        "published_genotype_call",
        "recalled_frozen_panel_call",
        "snippy_frozen_panel_call",
        "calls_concordant",
        "false_susceptible_mechanism",
    ]
    merged = ningxia.merge(manifold[manifold_cols], on=["dataset_id", "isolate_id"], how="left")
    merged = merged.merge(audit[audit_cols], on=["dataset_id", "isolate_id"], how="left")
    order = pd.Categorical(merged["outcome"], categories=["FN", "TP", "FP", "TN"], ordered=True)
    merged = merged.assign(_outcome_order=order).sort_values(
        ["_outcome_order", "mic_numeric", "isolate_id"], ascending=[True, False, True]
    )
    merged["atlas_order"] = np.arange(1, len(merged) + 1)
    return merged.drop(columns="_outcome_order")


def three_layer_summary(
    callability: pd.DataFrame,
    discriminator: pd.DataFrame,
    performance: pd.DataFrame,
    mic_summary: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict] = []
    for dataset in [NINGXIA, READ_COHORT]:
        for antibiotic in ["clarithromycin", "levofloxacin"]:
            item = callability[(callability["dataset_id"] == dataset) & (callability["antibiotic"] == antibiotic)].iloc[0]
            rows.append(
                {
                    "layer": "analytic_availability",
                    "dataset_id": dataset,
                    "antibiotic": antibiotic,
                    "measure": "target_callability",
                    "estimate": float(item["callability"]),
                    "detail": f"{int(item['n_callable'])}/{int(item['n_phenotype_linked'])}",
                }
            )
    for _, item in discriminator.iterrows():
        rows.append(
            {
                "layer": "population_shift",
                "dataset_id": item["external_dataset"],
                "antibiotic": "not_applicable",
                "measure": "PC1-PC10_cohort_discriminator_oof_AUC",
                "estimate": float(item["pooled_oof_auc"]),
                "detail": f"permutation P={item['permutation_p']:.4g}",
            }
        )
    dominant = performance[
        (performance["lineage_recomputed"] == DOMINANT_LINEAGE)
        & performance["metric"].isin(["sensitivity", "specificity", "marker_negative_resistance"])
    ]
    for _, item in dominant.iterrows():
        rows.append(
            {
                "layer": "conditional_transport",
                "dataset_id": item["dataset_id"],
                "antibiotic": "levofloxacin",
                "measure": f"{DOMINANT_LINEAGE}_{item['metric']}",
                "estimate": float(item["estimate"]),
                "detail": f"{int(item['successes'])}/{int(item['total'])}",
            }
        )
    fn = mic_summary[mic_summary["outcome"] == "FN"].iloc[0]
    rows.append(
        {
            "layer": "clinical_severity",
            "dataset_id": NINGXIA,
            "antibiotic": "levofloxacin",
            "measure": "false_susceptible_median_lower_bound_MIC_mg_L",
            "estimate": float(fn["median_lower_bound_mic_mg_L"]),
            "detail": f"n={int(fn['n'])}; right-censored={int(fn['right_censored_n'])}",
        }
    )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--bootstrap-repeats", type=int, default=2000)
    args = parser.parse_args()
    root = args.root.resolve()
    output = root / "results/transport_shift"
    output.mkdir(parents=True, exist_ok=True)

    predictions = pd.read_csv(root / "results/external_validation/sample_level_predictions.csv")
    lineage = pd.read_csv(root / "results/lineage_validation/lineage_assignments.csv")
    manifold = pd.read_csv(root / "results/extended_analysis/development_manifold_distance.csv")
    audit = pd.read_csv(root / "results/extended_analysis/ningxia_lvx_error_audit_enriched.csv")
    callability = pd.read_csv(root / "results/qc/panel_callability.csv")
    primary = primary_levofloxacin(predictions)

    performance = lineage_performance(primary)
    comparisons = dominant_lineage_comparisons(primary)
    standardized, standardization_weights = lineage_standardization(primary)
    group_map = primary[["dataset_id", "isolate_id", "near_clone_group"]].drop_duplicates()
    adv_predictions, adv_summary, adv_null, adv_coefficients = adversarial_validation(
        lineage,
        group_map,
        permutations=args.permutations,
        bootstrap_repeats=args.bootstrap_repeats,
    )
    severity, severity_summary = mic_severity(primary)
    atlas = failure_atlas(primary, manifold, audit)
    layers = three_layer_summary(callability, adv_summary, performance, severity_summary)

    outputs = {
        "levofloxacin_lineage_performance.csv": performance,
        "dominant_lineage_comparisons.csv": comparisons,
        "lineage_standardized_performance.csv": standardized,
        "lineage_standardization_weights.csv": standardization_weights,
        "cohort_discriminator_oof_predictions.csv": adv_predictions,
        "cohort_discriminator_summary.csv": adv_summary,
        "cohort_discriminator_permutations.csv": adv_null,
        "cohort_discriminator_coefficients.csv": adv_coefficients,
        "ningxia_resistant_mic_severity.csv": severity,
        "ningxia_resistant_mic_summary.csv": severity_summary,
        "ningxia_failure_atlas.csv": atlas,
        "three_layer_transport_summary.csv": layers,
    }
    for name, frame in outputs.items():
        frame.to_csv(output / name, index=False)
    print(f"Wrote {len(outputs)} transport-shift outputs to {output}")


if __name__ == "__main__":
    main()
