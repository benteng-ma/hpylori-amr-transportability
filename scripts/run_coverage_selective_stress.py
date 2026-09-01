#!/usr/bin/env python3
"""Post-freeze coverage-aware and selective-reporting stress analyses.

These analyses do not refit either frozen mutation catalogue. They quantify
the full phenotype-linked diagnostic cascade, logical performance bounds when
some isolates are unresolved, phenotype-blind manifold-distance abstention,
and between-external-cohort levofloxacin heterogeneity.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


DATASET_ORDER = ["HPGP_GLOBAL", "CHINA_NINGXIA_2022", "ZENODO_10369064"]
DRUG_ORDER = ["clarithromycin", "levofloxacin"]
ABSTENTION_CUTOFFS = [0.05, 0.10, 0.20, 0.50, 0.75, 0.90, 0.95, 0.99, 1.00]
Z_975 = 1.959963984540054


def safe_rate(successes: int, total: int) -> float:
    return successes / total if total else float("nan")


def confusion(frame: pd.DataFrame) -> dict[str, int | float]:
    resistant = frame["phenotype"].eq("R")
    susceptible = frame["phenotype"].eq("S")
    predicted_resistant = frame["prediction"].eq("R")
    predicted_susceptible = frame["prediction"].eq("S")
    tp = int((resistant & predicted_resistant).sum())
    fn = int((resistant & predicted_susceptible).sum())
    tn = int((susceptible & predicted_susceptible).sum())
    fp = int((susceptible & predicted_resistant).sum())
    sensitivity = safe_rate(tp, tp + fn)
    specificity = safe_rate(tn, tn + fp)
    return {
        "n": int(len(frame)),
        "n_resistant": int(tp + fn),
        "n_susceptible": int(tn + fp),
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "false_susceptible_rate": 1 - sensitivity if math.isfinite(sensitivity) else float("nan"),
        "balanced_accuracy": (sensitivity + specificity) / 2
        if math.isfinite(sensitivity) and math.isfinite(specificity)
        else float("nan"),
    }


def wilson_interval(successes: int, total: int) -> tuple[float, float]:
    if total <= 0:
        return float("nan"), float("nan")
    estimate = successes / total
    denominator = 1 + Z_975**2 / total
    centre = (estimate + Z_975**2 / (2 * total)) / denominator
    half = (
        Z_975
        * math.sqrt(estimate * (1 - estimate) / total + Z_975**2 / (4 * total**2))
        / denominator
    )
    return max(0.0, centre - half), min(1.0, centre + half)


def newcombe_difference_interval(
    success_a: int, total_a: int, success_b: int, total_b: int
) -> tuple[float, float]:
    """Newcombe hybrid-score interval for independent risk difference A-B."""
    rate_a = safe_rate(success_a, total_a)
    rate_b = safe_rate(success_b, total_b)
    low_a, high_a = wilson_interval(success_a, total_a)
    low_b, high_b = wilson_interval(success_b, total_b)
    estimate = rate_a - rate_b
    lower = estimate - math.sqrt((rate_a - low_a) ** 2 + (high_b - rate_b) ** 2)
    upper = estimate + math.sqrt((high_a - rate_a) ** 2 + (rate_b - low_b) ** 2)
    return lower, upper


def end_to_end_tables(predictions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    binary = predictions[predictions["phenotype"].isin(["R", "S"])].copy()
    yield_rows: list[dict[str, object]] = []
    bound_rows: list[dict[str, object]] = []
    for dataset_id in DATASET_ORDER:
        for antibiotic in DRUG_ORDER:
            group = binary[
                binary["dataset_id"].eq(dataset_id) & binary["antibiotic"].eq(antibiotic)
            ].copy()
            if group.empty:
                continue
            evaluable = group["analysis_status"].eq("PRIMARY")
            correct = group["correct"].eq("yes")
            resistant = group["phenotype"].eq("R")
            susceptible = group["phenotype"].eq("S")
            evaluated = group[evaluable]
            counts = confusion(evaluated)
            linked = int(len(group))
            evaluable_n = int(evaluable.sum())
            correct_n = int(correct.sum())
            unresolved_r = int((resistant & ~evaluable).sum())
            unresolved_s = int((susceptible & ~evaluable).sum())
            yield_rows.append(
                {
                    "dataset_id": dataset_id,
                    "antibiotic": antibiotic,
                    "phenotype_linked_binary_n": linked,
                    "qc_pass_n": int(group["basic_qc_status"].eq("PASS").sum()),
                    "target_callable_n": int(group["callable"].eq("yes").sum()),
                    "primary_evaluable_n": evaluable_n,
                    "correct_n": correct_n,
                    "incorrect_n": evaluable_n - correct_n,
                    "unresolved_n": linked - evaluable_n,
                    "resistant_total_n": int(resistant.sum()),
                    "susceptible_total_n": int(susceptible.sum()),
                    "tp": counts["tp"],
                    "tn": counts["tn"],
                    "fp": counts["fp"],
                    "fn": counts["fn"],
                    "evaluable_fraction": safe_rate(evaluable_n, linked),
                    "actionable_correct_yield": safe_rate(correct_n, linked),
                    "resistant_detection_yield": safe_rate(int(counts["tp"]), int(resistant.sum())),
                    "susceptible_correct_yield": safe_rate(int(counts["tn"]), int(susceptible.sum())),
                }
            )

            sensitivity_observed = float(counts["sensitivity"])
            specificity_observed = float(counts["specificity"])
            sensitivity_lower = safe_rate(int(counts["tp"]), int(resistant.sum()))
            sensitivity_upper = safe_rate(int(counts["tp"]) + unresolved_r, int(resistant.sum()))
            specificity_lower = safe_rate(int(counts["tn"]), int(susceptible.sum()))
            specificity_upper = safe_rate(int(counts["tn"]) + unresolved_s, int(susceptible.sum()))
            for metric, observed, lower, upper, denominator, resolved, unresolved in (
                (
                    "sensitivity",
                    sensitivity_observed,
                    sensitivity_lower,
                    sensitivity_upper,
                    int(resistant.sum()),
                    int(counts["tp"]) + int(counts["fn"]),
                    unresolved_r,
                ),
                (
                    "specificity",
                    specificity_observed,
                    specificity_lower,
                    specificity_upper,
                    int(susceptible.sum()),
                    int(counts["tn"]) + int(counts["fp"]),
                    unresolved_s,
                ),
                (
                    "balanced_accuracy",
                    float(counts["balanced_accuracy"]),
                    (sensitivity_lower + specificity_lower) / 2,
                    (sensitivity_upper + specificity_upper) / 2,
                    linked,
                    evaluable_n,
                    linked - evaluable_n,
                ),
            ):
                bound_rows.append(
                    {
                        "dataset_id": dataset_id,
                        "antibiotic": antibiotic,
                        "metric": metric,
                        "evaluable_estimate": observed,
                        "logical_lower_bound": lower,
                        "logical_upper_bound": upper,
                        "phenotype_linked_denominator": denominator,
                        "resolved_n": resolved,
                        "unresolved_n": unresolved,
                        "bound_definition": "unresolved isolates assigned to the least or most favourable compatible prediction",
                    }
                )
    return pd.DataFrame(yield_rows), pd.DataFrame(bound_rows)


def selective_abstention(predictions: pd.DataFrame, manifold: pd.DataFrame) -> pd.DataFrame:
    distances = manifold[["dataset_id", "isolate_id", "development_5nn_percentile"]].copy()
    merged = predictions.merge(distances, on=["dataset_id", "isolate_id"], how="inner")
    merged = merged[
        merged["analysis_status"].eq("PRIMARY") & merged["phenotype"].isin(["R", "S"])
    ].copy()
    rows: list[dict[str, object]] = []
    for dataset_id in DATASET_ORDER:
        for antibiotic in DRUG_ORDER:
            group = merged[
                merged["dataset_id"].eq(dataset_id) & merged["antibiotic"].eq(antibiotic)
            ].copy()
            if group.empty:
                continue
            for cutoff in ABSTENTION_CUTOFFS:
                accepted = group[group["development_5nn_percentile"].le(cutoff)].copy()
                counts = confusion(accepted)
                gate = (
                    counts["n_resistant"] >= 10
                    and counts["n_susceptible"] >= 10
                    and float(counts["sensitivity"]) >= 0.90
                    and float(counts["specificity"]) >= 0.90
                    and float(counts["false_susceptible_rate"]) <= 0.10
                )
                rows.append(
                    {
                        "dataset_id": dataset_id,
                        "antibiotic": antibiotic,
                        "development_percentile_cutoff": cutoff,
                        "threshold_source": "HpGP empirical distribution of phenotype-blind 5-nearest-neighbour distance",
                        "original_primary_n": int(len(group)),
                        "accepted_n": int(len(accepted)),
                        "rejected_n": int(len(group) - len(accepted)),
                        "retained_coverage": safe_rate(len(accepted), len(group)),
                        **counts,
                        "passes_frozen_safety_gate": "yes" if gate else "no",
                    }
                )
    return pd.DataFrame(rows)


def external_differences(
    predictions: pd.DataFrame, bootstrap: pd.DataFrame
) -> pd.DataFrame:
    primary = predictions[
        predictions["analysis_status"].eq("PRIMARY")
        & predictions["phenotype"].isin(["R", "S"])
        & predictions["antibiotic"].eq("levofloxacin")
    ].copy()
    ningxia = confusion(primary[primary["dataset_id"].eq("CHINA_NINGXIA_2022")])
    reads = confusion(primary[primary["dataset_id"].eq("ZENODO_10369064")])
    definitions = [
        ("sensitivity", int(ningxia["tp"]), int(ningxia["n_resistant"]), int(reads["tp"]), int(reads["n_resistant"])),
        ("specificity", int(ningxia["tn"]), int(ningxia["n_susceptible"]), int(reads["tn"]), int(reads["n_susceptible"])),
        ("false_susceptible_rate", int(ningxia["fn"]), int(ningxia["n_resistant"]), int(reads["fn"]), int(reads["n_resistant"])),
    ]
    rows: list[dict[str, object]] = []
    for metric, success_a, total_a, success_b, total_b in definitions:
        estimate = safe_rate(success_a, total_a) - safe_rate(success_b, total_b)
        low, high = newcombe_difference_interval(success_a, total_a, success_b, total_b)
        _, p_value = fisher_exact(
            [[success_a, total_a - success_a], [success_b, total_b - success_b]],
            alternative="two-sided",
        )
        rows.append(
            {
                "comparison": "Ningxia minus read cohort",
                "antibiotic": "levofloxacin",
                "metric": metric,
                "ningxia_estimate": safe_rate(success_a, total_a),
                "read_cohort_estimate": safe_rate(success_b, total_b),
                "absolute_difference": estimate,
                "difference_ci_low": low,
                "difference_ci_high": high,
                "ci_method": "Newcombe independent-proportion difference interval",
                "fisher_exact_p": p_value,
            }
        )

    boot = bootstrap[bootstrap["antibiotic"].eq("levofloxacin")].copy()
    a = boot[boot["dataset_id"].eq("CHINA_NINGXIA_2022")].sort_values("replicate")
    b = boot[boot["dataset_id"].eq("ZENODO_10369064")].sort_values("replicate")
    differences = a["balanced_accuracy"].to_numpy(float) - b["balanced_accuracy"].to_numpy(float)
    differences = differences[np.isfinite(differences)]
    rows.append(
        {
            "comparison": "Ningxia minus read cohort",
            "antibiotic": "levofloxacin",
            "metric": "balanced_accuracy",
            "ningxia_estimate": float(ningxia["balanced_accuracy"]),
            "read_cohort_estimate": float(reads["balanced_accuracy"]),
            "absolute_difference": float(ningxia["balanced_accuracy"]) - float(reads["balanced_accuracy"]),
            "difference_ci_low": float(np.quantile(differences, 0.025)),
            "difference_ci_high": float(np.quantile(differences, 0.975)),
            "ci_method": "empirical difference of independent 2,000 near-clone-group bootstrap distributions",
            "fisher_exact_p": float("nan"),
        }
    )
    return pd.DataFrame(rows)


def key_summary(
    yield_table: pd.DataFrame,
    bounds: pd.DataFrame,
    selective: pd.DataFrame,
    differences: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for dataset_id, antibiotic, metric in (
        ("CHINA_NINGXIA_2022", "clarithromycin", "actionable_correct_yield"),
        ("CHINA_NINGXIA_2022", "levofloxacin", "actionable_correct_yield"),
    ):
        row = yield_table[
            yield_table["dataset_id"].eq(dataset_id)
            & yield_table["antibiotic"].eq(antibiotic)
        ].iloc[0]
        rows.append(
            {
                "finding": f"{dataset_id}_{antibiotic}_{metric}",
                "estimate": row[metric],
                "detail": f"{int(row['correct_n'])}/{int(row['phenotype_linked_binary_n'])} phenotype-linked isolates yielded a correct primary result",
            }
        )
    for metric in ("sensitivity", "specificity"):
        row = bounds[
            bounds["dataset_id"].eq("CHINA_NINGXIA_2022")
            & bounds["antibiotic"].eq("clarithromycin")
            & bounds["metric"].eq(metric)
        ].iloc[0]
        rows.append(
            {
                "finding": f"CHINA_NINGXIA_2022_clarithromycin_{metric}_logical_bounds",
                "estimate": row["evaluable_estimate"],
                "detail": f"logical bounds {row['logical_lower_bound']:.6f}-{row['logical_upper_bound']:.6f}",
            }
        )
    ningxia_lvx = selective[
        selective["dataset_id"].eq("CHINA_NINGXIA_2022")
        & selective["antibiotic"].eq("levofloxacin")
        & selective["n_resistant"].ge(10)
        & selective["n_susceptible"].ge(10)
    ].copy()
    best = ningxia_lvx.sort_values(
        ["false_susceptible_rate", "retained_coverage"], ascending=[True, False]
    ).iloc[0]
    rows.append(
        {
            "finding": "CHINA_NINGXIA_2022_levofloxacin_best_eligible_manifold_abstention",
            "estimate": best["false_susceptible_rate"],
            "detail": f"cutoff={best['development_percentile_cutoff']:.2f}; retained_coverage={best['retained_coverage']:.6f}; safety_gate={best['passes_frozen_safety_gate']}",
        }
    )
    for row in differences.itertuples(index=False):
        rows.append(
            {
                "finding": f"external_levofloxacin_difference_{row.metric}",
                "estimate": row.absolute_difference,
                "detail": f"95% interval {row.difference_ci_low:.6f}-{row.difference_ci_high:.6f}; Fisher P={row.fisher_exact_p if math.isfinite(row.fisher_exact_p) else 'NE'}",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    output = root / "results/extended_analysis"
    output.mkdir(parents=True, exist_ok=True)
    predictions = pd.read_csv(root / "results/external_validation/sample_level_predictions.csv")
    manifold = pd.read_csv(root / "results/extended_analysis/development_manifold_distance.csv")
    bootstrap = pd.read_csv(root / "results/extended_analysis/external_bootstrap_distributions.csv")

    yield_table, bounds = end_to_end_tables(predictions)
    selective = selective_abstention(predictions, manifold)
    differences = external_differences(predictions, bootstrap)
    summary = key_summary(yield_table, bounds, selective, differences)

    yield_table.to_csv(output / "coverage_aware_yield.csv", index=False)
    bounds.to_csv(output / "callability_identification_bounds.csv", index=False)
    selective.to_csv(output / "manifold_abstention_metrics.csv", index=False)
    differences.to_csv(output / "external_performance_differences.csv", index=False)
    summary.to_csv(output / "coverage_selective_summary.csv", index=False)


if __name__ == "__main__":
    main()
