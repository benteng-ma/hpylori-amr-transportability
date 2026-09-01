#!/usr/bin/env python3
"""Generate explanatory and robustness analyses for the expanded manuscript.

The frozen marker catalogues, phenotype labels, diagnostic denominators, and
transportability gates are never refit here.  All post-freeze analyses are
descriptive, error-decomposing, or explicitly labelled exploratory.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import fisher_exact, mannwhitneyu
from sklearn.metrics import pairwise_distances, roc_auc_score


COHORT_ORDER = ["HPGP_GLOBAL", "CHINA_NINGXIA_2022", "ZENODO_10369064"]
EXTERNAL_COHORTS = ["CHINA_NINGXIA_2022", "ZENODO_10369064"]
DRUGS = ["clarithromycin", "levofloxacin"]
FROZEN_LVX = {"N87K", "N87I", "D91G", "D91N", "D91Y", "A88V", "A88P"}
FROZEN_CLR = {"A2142G", "A2142C", "A2143G"}
AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",
}


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "pass"}


def to_numeric(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    result = frame.copy()
    for column in columns:
        if column in result:
            result[column] = pd.to_numeric(result[column], errors="coerce")
    return result


def outcome(phenotype: str, prediction: str) -> str:
    mapping = {
        ("R", "R"): "TP",
        ("S", "S"): "TN",
        ("S", "R"): "FP",
        ("R", "S"): "FN",
    }
    return mapping.get((phenotype, prediction), "NOT_EVALUABLE")


def frozen_mutations(antibiotic: str, marker_summary: str) -> list[str]:
    """Return frozen resistance markers from heterogeneous caller summaries."""
    if antibiotic == "clarithromycin" and marker_summary.startswith("A2142:"):
        fractions: dict[str, float] = {}
        for token in marker_summary.split(";"):
            name, value = token.split(":", 1)
            fractions[name] = float(value)
        markers = []
        if fractions.get("A2142", 0.0) >= 0.20:
            markers.append("A2142 resistant allele")
        if fractions.get("A2143", 0.0) >= 0.20:
            markers.append("A2143G")
        return markers
    candidates = [token.strip() for token in marker_summary.split(";") if token.strip()]
    frozen = FROZEN_CLR if antibiotic == "clarithromycin" else FROZEN_LVX
    return [token for token in candidates if token in frozen]


def compute_metrics(frame: pd.DataFrame) -> dict[str, float]:
    phenotype = frame["phenotype"].to_numpy()
    prediction = frame["prediction"].to_numpy()
    tp = int(np.sum((phenotype == "R") & (prediction == "R")))
    tn = int(np.sum((phenotype == "S") & (prediction == "S")))
    fp = int(np.sum((phenotype == "S") & (prediction == "R")))
    fn = int(np.sum((phenotype == "R") & (prediction == "S")))
    sensitivity = tp / (tp + fn) if tp + fn else math.nan
    specificity = tn / (tn + fp) if tn + fp else math.nan
    return {
        "n": len(frame), "tp": tp, "tn": tn, "fp": fp, "fn": fn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "false_susceptible_rate": fn / (tp + fn) if tp + fn else math.nan,
        "false_resistant_rate": fp / (tn + fp) if tn + fp else math.nan,
        "balanced_accuracy": (sensitivity + specificity) / 2
        if not (math.isnan(sensitivity) or math.isnan(specificity)) else math.nan,
    }


def bh_adjust(values: list[float]) -> list[float]:
    array = np.asarray(values, dtype=float)
    order = np.argsort(array)
    adjusted = np.empty(len(array), dtype=float)
    running = 1.0
    for rank_index in range(len(array) - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, array[original_index] * len(array) / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted.tolist()


def marker_architecture(predictions: pd.DataFrame, output: Path) -> None:
    eligible = predictions[
        predictions["phenotype"].isin(["R", "S"])
        & predictions["prediction"].isin(["R", "S"])
    ].copy()
    eligible["outcome"] = [outcome(p, g) for p, g in zip(eligible["phenotype"], eligible["prediction"])]
    eligible["frozen_markers"] = [
        frozen_mutations(drug, summary)
        for drug, summary in zip(eligible["antibiotic"], eligible["marker_summary"])
    ]
    eligible["marker_combination"] = [";".join(markers) if markers else "WT" for markers in eligible["frozen_markers"]]
    rows = []
    for record in eligible.itertuples(index=False):
        markers = record.frozen_markers or ["WT"]
        for marker in markers:
            rows.append({
                "dataset_id": record.dataset_id,
                "isolate_id": record.isolate_id,
                "antibiotic": record.antibiotic,
                "phenotype": record.phenotype,
                "prediction": record.prediction,
                "outcome": record.outcome,
                "mutation": marker,
            })
    sample_markers = pd.DataFrame(rows)
    sample_markers.to_csv(output / "mutation_spectrum_samples.csv", index=False)
    spectrum = (
        sample_markers.groupby(
            ["dataset_id", "antibiotic", "phenotype", "prediction", "outcome", "mutation"],
            dropna=False,
        ).size().rename("n").reset_index()
    )
    spectrum.to_csv(output / "mutation_spectrum_summary.csv", index=False)
    combinations = (
        eligible.groupby(
            ["dataset_id", "antibiotic", "phenotype", "prediction", "outcome", "marker_combination"],
            dropna=False,
        ).size().rename("n").reset_index()
    )
    combinations.to_csv(output / "marker_combination_summary.csv", index=False)

    prevalence_rows = []
    for (dataset, drug), part in eligible.groupby(["dataset_id", "antibiotic"], sort=False):
        prevalence_rows.append({
            "dataset_id": dataset,
            "antibiotic": drug,
            "n": len(part),
            "phenotypic_resistance_prevalence": float((part["phenotype"] == "R").mean()),
            "marker_resistance_prevalence": float((part["prediction"] == "R").mean()),
            "prevalence_gap_marker_minus_phenotype": float(
                (part["prediction"] == "R").mean() - (part["phenotype"] == "R").mean()
            ),
        })
    pd.DataFrame(prevalence_rows).to_csv(output / "marker_phenotype_prevalence_shift.csv", index=False)

    lineage = (
        eligible.groupby(["dataset_id", "antibiotic", "lineage_recomputed"], dropna=False)
        .agg(
            n=("isolate_id", "size"),
            phenotype_resistant=("phenotype", lambda x: int((x == "R").sum())),
            marker_resistant=("prediction", lambda x: int((x == "R").sum())),
            errors=("outcome", lambda x: int(x.isin(["FN", "FP"]).sum())),
        ).reset_index()
    )
    lineage["phenotype_resistance_prevalence"] = lineage["phenotype_resistant"] / lineage["n"]
    lineage["marker_resistance_prevalence"] = lineage["marker_resistant"] / lineage["n"]
    lineage["error_rate"] = lineage["errors"] / lineage["n"]
    lineage.to_csv(output / "marker_prevalence_by_lineage.csv", index=False)


def qc_callability(predictions: pd.DataFrame, qc: pd.DataFrame, output: Path) -> None:
    numeric_columns = [
        "contigs", "assembly_size_bp", "n50_bp", "longest_contig_bp", "gc_percent",
        "ambiguous_bases", "completeness_percent", "contamination_percent",
    ]
    qc = to_numeric(qc, numeric_columns)
    joined = predictions.merge(qc, on=["dataset_id", "isolate_id"], how="left", suffixes=("", "_qc"))
    joined["callable_binary"] = joined["callable"].map(truthy).astype(int)
    joined["log10_n50_bp"] = np.log10(joined["n50_bp"].where(joined["n50_bp"] > 0))
    keep = [
        "dataset_id", "isolate_id", "antibiotic", "callable", "callable_binary",
        "final_qc_status", "contigs", "assembly_size_bp", "n50_bp", "log10_n50_bp",
        "longest_contig_bp", "completeness_percent", "contamination_percent",
    ]
    joined[keep].to_csv(output / "qc_callability_samples.csv", index=False)

    association_rows = []
    metrics = ["contigs", "log10_n50_bp", "completeness_percent", "contamination_percent"]
    for (dataset, drug), part in joined.groupby(["dataset_id", "antibiotic"], sort=False):
        for metric in metrics:
            available = part[["callable_binary", metric]].dropna()
            yes = available.loc[available["callable_binary"] == 1, metric].to_numpy(dtype=float)
            no = available.loc[available["callable_binary"] == 0, metric].to_numpy(dtype=float)
            row = {
                "dataset_id": dataset, "antibiotic": drug, "metric": metric,
                "n_callable": len(yes), "n_uncallable": len(no),
                "callable_median": float(np.median(yes)) if len(yes) else math.nan,
                "uncallable_median": float(np.median(no)) if len(no) else math.nan,
                "mann_whitney_u": math.nan, "mann_whitney_p": math.nan,
                "cliffs_delta_callable_minus_uncallable": math.nan,
                "auc_metric_predicting_callable": math.nan,
            }
            if len(yes) >= 3 and len(no) >= 3:
                test = mannwhitneyu(yes, no, alternative="two-sided")
                row["mann_whitney_u"] = float(test.statistic)
                row["mann_whitney_p"] = float(test.pvalue)
                row["cliffs_delta_callable_minus_uncallable"] = float(
                    2 * test.statistic / (len(yes) * len(no)) - 1
                )
                labels = available["callable_binary"].to_numpy(dtype=int)
                values = available[metric].to_numpy(dtype=float)
                row["auc_metric_predicting_callable"] = float(roc_auc_score(labels, values))
            association_rows.append(row)
    pd.DataFrame(association_rows).to_csv(output / "qc_callability_associations.csv", index=False)


def manifold_analysis(predictions: pd.DataFrame, lineages: pd.DataFrame, output: Path) -> None:
    pc_columns = [f"PC{i}" for i in range(1, 11)]
    lineages = to_numeric(lineages, pc_columns)
    development = lineages[lineages["dataset_id"] == "HPGP_GLOBAL"].copy()
    mean = development[pc_columns].mean()
    sd = development[pc_columns].std(ddof=1).replace(0, 1)
    scaled = (lineages[pc_columns] - mean) / sd
    dev_scaled = scaled.loc[development.index]
    distances = pairwise_distances(scaled.to_numpy(), dev_scaled.to_numpy(), metric="euclidean")
    dev_positions = {index: position for position, index in enumerate(development.index)}
    for row_position, original_index in enumerate(lineages.index):
        if original_index in dev_positions:
            distances[row_position, dev_positions[original_index]] = np.inf
    k = 5
    nearest = np.partition(distances, k - 1, axis=1)[:, :k]
    lineages = lineages.copy()
    lineages["development_5nn_distance"] = nearest.mean(axis=1)
    hp_scores = lineages.loc[lineages["dataset_id"] == "HPGP_GLOBAL", "development_5nn_distance"].to_numpy()
    lineages["development_5nn_percentile"] = [float(np.mean(hp_scores <= score)) for score in lineages["development_5nn_distance"]]
    lineages.to_csv(output / "development_manifold_distance.csv", index=False)

    summary_rows = []
    hp = lineages[lineages["dataset_id"] == "HPGP_GLOBAL"]["development_5nn_distance"].to_numpy()
    clusters = sorted(lineages["lineage_recomputed"].unique())
    hp_counts = lineages[lineages["dataset_id"] == "HPGP_GLOBAL"]["lineage_recomputed"].value_counts().reindex(clusters, fill_value=0).to_numpy(dtype=float)
    hp_distribution = hp_counts / hp_counts.sum()
    for dataset in COHORT_ORDER:
        values = lineages[lineages["dataset_id"] == dataset]["development_5nn_distance"].to_numpy()
        counts = lineages[lineages["dataset_id"] == dataset]["lineage_recomputed"].value_counts().reindex(clusters, fill_value=0).to_numpy(dtype=float)
        distribution = counts / counts.sum()
        row = {
            "dataset_id": dataset, "n": len(values),
            "median_5nn_distance": float(np.median(values)),
            "q1_5nn_distance": float(np.quantile(values, 0.25)),
            "q3_5nn_distance": float(np.quantile(values, 0.75)),
            "median_development_percentile": float(np.median(
                lineages[lineages["dataset_id"] == dataset]["development_5nn_percentile"]
            )),
            "mann_whitney_vs_hpgp_p": math.nan,
            "cliffs_delta_vs_hpgp": math.nan,
            "lineage_jensen_shannon_distance_vs_hpgp": float(jensenshannon(distribution, hp_distribution, base=2)),
        }
        if dataset != "HPGP_GLOBAL":
            test = mannwhitneyu(values, hp, alternative="two-sided")
            row["mann_whitney_vs_hpgp_p"] = float(test.pvalue)
            row["cliffs_delta_vs_hpgp"] = float(2 * test.statistic / (len(values) * len(hp)) - 1)
        summary_rows.append(row)
    pd.DataFrame(summary_rows).to_csv(output / "development_manifold_summary.csv", index=False)

    scored = predictions.merge(
        lineages[["dataset_id", "isolate_id", "development_5nn_distance", "development_5nn_percentile"]],
        on=["dataset_id", "isolate_id"], how="left",
    )
    scored = scored[
        scored["dataset_id"].isin(EXTERNAL_COHORTS)
        & scored["phenotype"].isin(["R", "S"])
        & scored["prediction"].isin(["R", "S"])
    ].copy()
    scored["outcome"] = [outcome(p, g) for p, g in zip(scored["phenotype"], scored["prediction"])]
    scored["error"] = scored["outcome"].isin(["FN", "FP"])
    scored.to_csv(output / "external_error_manifold_samples.csv", index=False)
    error_rows = []
    for (dataset, drug), part in scored.groupby(["dataset_id", "antibiotic"], sort=False):
        correct = part.loc[~part["error"], "development_5nn_distance"].dropna().to_numpy()
        errors = part.loc[part["error"], "development_5nn_distance"].dropna().to_numpy()
        row = {
            "dataset_id": dataset, "antibiotic": drug,
            "n_correct": len(correct), "n_error": len(errors),
            "median_correct": float(np.median(correct)) if len(correct) else math.nan,
            "median_error": float(np.median(errors)) if len(errors) else math.nan,
            "mann_whitney_p": math.nan, "cliffs_delta_error_minus_correct": math.nan,
        }
        if len(correct) >= 3 and len(errors) >= 3:
            test = mannwhitneyu(errors, correct, alternative="two-sided")
            row["mann_whitney_p"] = float(test.pvalue)
            row["cliffs_delta_error_minus_correct"] = float(
                2 * test.statistic / (len(errors) * len(correct)) - 1
            )
        error_rows.append(row)
    pd.DataFrame(error_rows).to_csv(output / "external_error_manifold_associations.csv", index=False)


def parse_missense(effect: str) -> str | None:
    match = re.search(r"p\.([A-Z][a-z]{2}|Ter)(\d+)([A-Z][a-z]{2}|Ter)", effect)
    if not match:
        return None
    ref, position, alt = match.groups()
    if ref not in AA3_TO_1 or alt not in AA3_TO_1:
        return None
    return f"{AA3_TO_1[ref]}{position}{AA3_TO_1[alt]}"


def read_gyrA_variants(tab_path: Path) -> list[str]:
    variants: set[str] = set()
    with tab_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("GENE") != "gyrA" or "missense_variant" not in row.get("EFFECT", ""):
                continue
            variant = parse_missense(row["EFFECT"])
            if variant:
                variants.add(variant)
    return sorted(variants, key=lambda value: (int(re.search(r"\d+", value).group()), value))


def gyra_residual_scan(root: Path, predictions: pd.DataFrame, qc: pd.DataFrame, output: Path) -> None:
    samples = qc[
        (qc["dataset_id"] == "CHINA_NINGXIA_2022")
        & (qc["final_qc_status"] == "PASS")
    ][["dataset_id", "isolate_id"]].copy()
    variant_map: dict[str, list[str]] = {}
    for isolate in samples["isolate_id"]:
        tab = root / "data/phylogeny/snippy" / f"CHINA_NINGXIA_2022__{isolate}" / "snps.tab"
        variant_map[isolate] = read_gyrA_variants(tab)

    audit = read_csv(root / "results/external_validation/ningxia_published_genotype_audit_samples.csv")
    audit = audit[audit["antibiotic"] == "levofloxacin"].copy()
    lvx = predictions[
        (predictions["dataset_id"] == "CHINA_NINGXIA_2022")
        & (predictions["antibiotic"] == "levofloxacin")
        & predictions["phenotype"].isin(["R", "S"])
        & predictions["prediction"].isin(["R", "S"])
    ].copy()
    lvx["outcome"] = [outcome(p, g) for p, g in zip(lvx["phenotype"], lvx["prediction"])]
    lvx["full_gyrA_missense_variants"] = lvx["isolate_id"].map(lambda value: ";".join(variant_map.get(value, [])))
    lvx = lvx.merge(
        audit[[
            "isolate_id", "published_genotype_call", "recalled_frozen_panel_call",
            "recalled_target_substitutions", "snippy_qrdr_substitutions",
            "snippy_frozen_panel_call", "calls_concordant",
        ]], on="isolate_id", how="left",
    )
    lvx["false_susceptible_mechanism"] = "not_false_susceptible"
    fn = lvx["outcome"] == "FN"
    lvx.loc[fn & lvx["recalled_target_substitutions"].str.contains("N87Y", regex=False), "false_susceptible_mechanism"] = "off_panel_N87Y"
    lvx.loc[fn & (lvx["recalled_target_substitutions"] == ""), "false_susceptible_mechanism"] = "QRDR_wild_type_in_deposited_assembly"
    remaining = fn & (lvx["false_susceptible_mechanism"] == "not_false_susceptible")
    lvx.loc[remaining, "false_susceptible_mechanism"] = "other_off_panel_QRDR"
    lvx.to_csv(output / "ningxia_lvx_error_audit_enriched.csv", index=False)
    (
        lvx[lvx["outcome"] == "FN"].groupby("false_susceptible_mechanism")
        .size().rename("n").reset_index()
        .assign(proportion=lambda frame: frame["n"] / frame["n"].sum())
        .to_csv(output / "ningxia_lvx_false_susceptible_decomposition.csv", index=False)
    )

    all_variants = sorted({variant for variants in variant_map.values() for variant in variants}, key=lambda value: (int(re.search(r"\d+", value).group()), value))
    association_rows = []
    for variant in all_variants:
        present = lvx["full_gyrA_missense_variants"].map(lambda value: variant in value.split(";") if value else False)
        resistant = lvx["phenotype"] == "R"
        r_with = int((present & resistant).sum())
        r_without = int((~present & resistant).sum())
        s_with = int((present & ~resistant).sum())
        s_without = int((~present & ~resistant).sum())
        odds_ratio, p_value = fisher_exact([[r_with, r_without], [s_with, s_without]], alternative="two-sided")
        association_rows.append({
            "variant": variant,
            "position": int(re.search(r"\d+", variant).group()),
            "frozen_catalogue_marker": variant in FROZEN_LVX,
            "resistant_with": r_with, "resistant_without": r_without,
            "susceptible_with": s_with, "susceptible_without": s_without,
            "false_susceptible_with": int((present & (lvx["outcome"] == "FN")).sum()),
            "true_positive_with": int((present & (lvx["outcome"] == "TP")).sum()),
            "odds_ratio": float(odds_ratio), "fisher_p": float(p_value),
        })
    associations = pd.DataFrame(association_rows)
    if len(associations):
        associations["bh_q"] = bh_adjust(associations["fisher_p"].tolist())
        associations = associations.sort_values(["bh_q", "fisher_p", "position"])
    associations.to_csv(output / "ningxia_gyrA_missense_associations.csv", index=False)

    sample_variant_rows = []
    for isolate, variants in variant_map.items():
        for variant in variants:
            sample_variant_rows.append({"isolate_id": isolate, "variant": variant})
    pd.DataFrame(sample_variant_rows, columns=["isolate_id", "variant"]).to_csv(
        output / "ningxia_gyrA_missense_samples.csv", index=False
    )

    expanded = lvx.copy()
    expanded["prediction"] = np.where(
        (expanded["prediction"] == "S")
        & expanded["recalled_target_substitutions"].str.contains("N87Y", regex=False),
        "R", expanded["prediction"],
    )
    frozen_metrics = compute_metrics(lvx)
    expanded_metrics = compute_metrics(expanded)
    pd.DataFrame([
        {"panel": "frozen_catalogue", **frozen_metrics},
        {"panel": "exploratory_add_N87Y", **expanded_metrics},
    ]).to_csv(output / "exploratory_N87Y_rescue.csv", index=False)


def robustness(predictions: pd.DataFrame, metrics: pd.DataFrame, output: Path) -> None:
    rng = np.random.default_rng(20260831)
    bootstrap_rows = []
    for dataset in EXTERNAL_COHORTS:
        for drug in DRUGS:
            part = predictions[
                (predictions["dataset_id"] == dataset)
                & (predictions["antibiotic"] == drug)
                & predictions["phenotype"].isin(["R", "S"])
                & predictions["prediction"].isin(["R", "S"])
            ].copy()
            if not len(part):
                continue
            units = part["bootstrap_unit"].unique()
            by_unit = {unit: part[part["bootstrap_unit"] == unit] for unit in units}
            for replicate in range(2000):
                sampled_units = rng.choice(units, size=len(units), replace=True)
                sampled = pd.concat([by_unit[unit] for unit in sampled_units], ignore_index=True)
                values = compute_metrics(sampled)
                bootstrap_rows.append({
                    "dataset_id": dataset, "antibiotic": drug,
                    "replicate": replicate + 1, **values,
                })
    pd.DataFrame(bootstrap_rows).to_csv(output / "external_bootstrap_distributions.csv", index=False)

    ningxia = predictions[
        (predictions["dataset_id"] == "CHINA_NINGXIA_2022")
        & (predictions["antibiotic"] == "levofloxacin")
        & predictions["phenotype"].isin(["R", "S"])
        & predictions["prediction"].isin(["R", "S"])
    ].copy()
    influence_rows = [{"dropped_isolate": "NONE_BASELINE", **compute_metrics(ningxia)}]
    for isolate in ningxia["isolate_id"]:
        influence_rows.append({
            "dropped_isolate": isolate,
            **compute_metrics(ningxia[ningxia["isolate_id"] != isolate]),
        })
    pd.DataFrame(influence_rows).to_csv(output / "ningxia_lvx_leave_one_out_influence.csv", index=False)

    cohort_metrics = to_numeric(metrics[metrics["stratum_type"] == "COHORT"].copy(), [
        "n_resistant", "n_susceptible", "sensitivity", "specificity", "false_susceptible_rate",
    ])
    prevalence_rows = []
    for _, row in cohort_metrics[cohort_metrics["dataset_id"].isin(EXTERNAL_COHORTS)].iterrows():
        sensitivity = row["sensitivity"]
        specificity = row["specificity"]
        if pd.isna(sensitivity) or pd.isna(specificity):
            continue
        for prevalence in np.linspace(0.05, 0.75, 71):
            ppv_denominator = sensitivity * prevalence + (1 - specificity) * (1 - prevalence)
            npv_denominator = specificity * (1 - prevalence) + (1 - sensitivity) * prevalence
            prevalence_rows.append({
                "dataset_id": row["dataset_id"], "antibiotic": row["antibiotic"],
                "observed_n_resistant": row["n_resistant"], "observed_n_susceptible": row["n_susceptible"],
                "assumed_prevalence": prevalence,
                "ppv": sensitivity * prevalence / ppv_denominator if ppv_denominator else math.nan,
                "npv": specificity * (1 - prevalence) / npv_denominator if npv_denominator else math.nan,
            })
    pd.DataFrame(prevalence_rows).to_csv(output / "predictive_values_by_assumed_prevalence.csv", index=False)

    grid_rows = []
    for drug in DRUGS:
        drug_metrics = cohort_metrics[
            (cohort_metrics["antibiotic"] == drug)
            & cohort_metrics["dataset_id"].isin(EXTERNAL_COHORTS)
        ]
        for minimum in [2, 5, 10, 15, 20]:
            eligible = drug_metrics[
                (drug_metrics["n_resistant"] >= minimum)
                & (drug_metrics["n_susceptible"] >= minimum)
            ]
            for threshold in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]:
                if len(eligible) < 2:
                    label = "INSUFFICIENT_DATA"
                elif (eligible["false_susceptible_rate"] > threshold).any():
                    label = "HIGH_FALSE_SUSCEPTIBLE_RISK"
                else:
                    label = "PASSES_EXPLORATORY_GRID"
                grid_rows.append({
                    "antibiotic": drug, "minimum_resistant_and_susceptible_per_external_cohort": minimum,
                    "false_susceptible_gate": threshold, "eligible_external_cohorts": len(eligible),
                    "classification": label,
                })
    pd.DataFrame(grid_rows).to_csv(output / "transport_gate_robustness_grid.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=Path("results/extended_analysis"))
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output if args.output.is_absolute() else root / args.output
    output.mkdir(parents=True, exist_ok=True)

    predictions = read_csv(root / "results/external_validation/sample_level_predictions.csv")
    qc = read_csv(root / "results/qc/assembly_qc_with_checkm2.csv")
    lineages = read_csv(root / "results/lineage_validation/lineage_assignments.csv")
    metrics = read_csv(root / "results/external_validation/frozen_panel_metrics.csv")

    marker_architecture(predictions, output)
    qc_callability(predictions, qc, output)
    manifold_analysis(predictions, lineages, output)
    gyra_residual_scan(root, predictions, qc, output)
    robustness(predictions, metrics, output)


if __name__ == "__main__":
    main()
