#!/usr/bin/env python3
"""Build reviewer-readable supplementary tables for the Word/PDF supplement.

The complete row-level records remain in Supplementary_Data_S1-S20.xlsx.  This
builder deliberately prints compact, decision-relevant summaries in the
paginated supplement while deriving every displayed value from the same
versioned CSV outputs used to create the workbook.  Worksheet S12 is a
machine-readable figure-source manifest and is intentionally not printed as a
reader-facing result table.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


COHORT = {
    "HPGP_GLOBAL": "HpGP",
    "CHINA_NINGXIA_2022": "Ningxia",
    "ZENODO_10369064": "Read cohort",
    "MULTI_COHORT": "Multiple cohorts",
}
DRUG = {"clarithromycin": "Clarithromycin", "levofloxacin": "Levofloxacin"}


def read_csv(root: Path, relative: str) -> list[dict[str, str]]:
    with (root / relative).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def clean(value: object) -> str:
    if value is None:
        return "NE"
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return "NE"
    return text.replace("|", "/").replace("\n", " ")


def number(value: object, digits: int = 2) -> str:
    try:
        parsed = float(str(value))
    except (TypeError, ValueError):
        return clean(value)
    if not math.isfinite(parsed):
        return "NE"
    if abs(parsed - round(parsed)) < 1e-12:
        return str(int(round(parsed)))
    return f"{parsed:.{digits}f}".rstrip("0").rstrip(".")


def percent(value: object, digits: int = 1) -> str:
    try:
        parsed = float(str(value))
    except (TypeError, ValueError):
        return "NE"
    if not math.isfinite(parsed):
        return "NE"
    return f"{100 * parsed:.{digits}f}%"


def pvalue(value: object) -> str:
    try:
        parsed = float(str(value))
    except (TypeError, ValueError):
        return "NE"
    if not math.isfinite(parsed):
        return "NE"
    return f"{parsed:.2e}" if parsed < 0.001 else f"{parsed:.3f}"


def pct_ci(row: dict[str, str], metric: str) -> str:
    estimate = percent(row.get(metric))
    low = percent(row.get(f"{metric}_ci_low"))
    high = percent(row.get(f"{metric}_ci_high"))
    return estimate if low == "NE" else f"{estimate} ({low}-{high})"


def label(value: str) -> str:
    return clean(value).replace("_", " ").title()


def cohort(value: str) -> str:
    return COHORT.get(value, clean(value))


def drug(value: str) -> str:
    return DRUG.get(value, clean(value))


def drug_short(value: str) -> str:
    return {"clarithromycin": "CLR", "levofloxacin": "LVX"}.get(value, clean(value))


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def interval(values: list[float]) -> str:
    if not values:
        return "NE"
    return (
        f"{percent(quantile(values, 0.5))} "
        f"({percent(quantile(values, 0.025))}-{percent(quantile(values, 0.975))})"
    )


def md_table(title: str, headers: list[str], rows: list[list[object]], note: str | None = None) -> str:
    lines = [f"### {title}", "", "| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(clean(value) for value in row) + " |")
    if note:
        lines.extend(["", f"*Note:* {note}"])
    return "\n".join(lines)


def model_name(value: str) -> str:
    return {
        "LINEAGE_ONLY_LOGISTIC": "Lineage only",
        "MUTATION_ONLY_LOGISTIC": "Mutation only",
        "MUTATION_PLUS_LINEAGE_LOGISTIC": "Mutation + lineage",
    }.get(value, label(value))


def build(root: Path, output: Path) -> None:
    blocks: list[str] = [
        "## Supplementary tables",
        "",
        "The tables below provide reviewer-readable results within this Word/PDF supplement. Complete row-level and replicate-level data remain in Supplementary Data S1-S20.xlsx; the workbook is the machine-readable companion to, not a replacement for, the tables printed here. Worksheet S12 is a provenance manifest for figure source-data files and is retained only in the machine-readable workbook and repository; because it contains file inventory metadata rather than scientific results, it is not reproduced as a paginated table. NE, not estimable; CLR, clarithromycin; LVX, levofloxacin; FSR, false-susceptible rate; BA, balanced accuracy; QC, quality control.",
    ]

    # S1: cohort-level crosswalk summary.
    rows = read_csv(root, "results/external_validation/isolate_sequence_phenotype_crosswalk.csv")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["dataset_id"]].append(row)
    shown = []
    for dataset_id in ("HPGP_GLOBAL", "CHINA_NINGXIA_2022", "ZENODO_10369064"):
        items = grouped[dataset_id]
        countries = sorted({clean(item.get("country")) for item in items if clean(item.get("country")) != "NE"})
        geography = countries[0] if len(countries) == 1 else f"Multiple ({len(countries)})"
        years = sorted({int(float(item["collection_year"])) for item in items if clean(item.get("collection_year")) != "NE"})
        year_range = "NE" if not years else (str(years[0]) if years[0] == years[-1] else f"{years[0]}-{years[-1]}")
        shown.append([
            cohort(dataset_id), geography, len(items),
            sum("PASS" in item.get("final_qc_status", "").upper() for item in items),
            sum(item.get("clarithromycin_phenotype") in {"R", "S", "I"} for item in items),
            sum(item.get("levofloxacin_phenotype") in {"R", "S", "I"} for item in items),
            year_range,
        ])
    blocks.append(md_table(
        "Table S1. Isolate-sequence-phenotype-prediction crosswalk: cohort summary",
        ["Cohort", "Geography", "Materialized", "Final-QC pass", "CLR phenotype", "LVX phenotype", "Collection years"],
        shown,
        "The complete 526-row isolate-level crosswalk, accessions, checksums, QC fields, phenotypes, marker calls and predictions are provided in workbook sheet S01_Crosswalk.",
    ))

    # S2: all panel callability rows.
    rows = read_csv(root, "results/qc/panel_callability.csv")
    shown = [[cohort(r["dataset_id"]), drug(r["antibiotic"]), r["n_phenotype_linked"], r["n_callable"], r["n_uncallable"], percent(r["callability"])] for r in rows]
    blocks.append(md_table("Table S2. Complete genome quality and target callability", ["Cohort", "Drug", "Phenotype-linked", "Callable", "Uncallable", "Callability"], shown))

    # S3: lineage composition and near-clone pairs.
    lineage_rows = read_csv(root, "results/lineage_validation/lineage_assignments.csv")
    counts = Counter((r["dataset_id"], r["lineage_recomputed"]) for r in lineage_rows)
    totals = Counter(r["dataset_id"] for r in lineage_rows)
    shown = [[cohort(ds), lineage, count, percent(count / totals[ds])] for (ds, lineage), count in sorted(counts.items())]
    blocks.append(md_table("Table S3a. Phenotype-blind core-SNP lineage composition", ["Cohort", "Fixed SNP cluster", "n", "Within-cohort proportion"], shown))
    pairs = read_csv(root, "results/qc/pairwise_relatedness_candidates.csv")
    shown = [[f"{r['isolate_a']} / {r['isolate_b']}", number(r["mash_distance"], 6), number(r["skani_ani_percent"], 3), percent(r["skani_align_fraction_a"]), percent(r["skani_align_fraction_b"])] for r in pairs]
    blocks.append(md_table("Table S3b. Near-clone pairs meeting the frozen conjunction", ["Isolate pair", "Mash distance", "ANI (%)", "Aligned A", "Aligned B"], shown, "All four pairs were within HpGP; no near-clone component crossed an internal/external cohort boundary."))

    # S4: prediction accounting.
    rows = read_csv(root, "results/external_validation/sample_level_predictions.csv")
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["dataset_id"], row["antibiotic"])].append(row)
    shown = []
    for (dataset_id, antibiotic), items in sorted(grouped.items()):
        primary = [r for r in items if r.get("analysis_status") == "PRIMARY"]
        shown.append([cohort(dataset_id), drug(antibiotic), len(items), len(primary), sum(r.get("prediction") == "R" for r in primary), sum(r.get("prediction") == "S" for r in primary), sum(r.get("correct") == "yes" for r in primary), sum(r.get("correct") == "no" for r in primary)])
    blocks.append(md_table("Table S4. Frozen-panel sample-prediction accounting", ["Cohort", "Drug", "All records", "Primary", "Pred R", "Pred S", "Correct", "Errors"], shown, "The full 1,052-row sample-level table is provided in workbook sheet S04_Predictions."))

    # S5: primary cohort performance.
    rows = [r for r in read_csv(root, "results/external_validation/frozen_panel_metrics.csv") if r.get("stratum_type") == "COHORT"]
    shown = [[cohort(r["dataset_id"]), drug(r["antibiotic"]), f"{r['n']} ({r['n_resistant']}/{r['n_susceptible']})", pct_ci(r, "sensitivity"), pct_ci(r, "specificity"), pct_ci(r, "false_susceptible_rate"), percent(r["balanced_accuracy"])] for r in rows]
    blocks.append(md_table("Table S5. Frozen catalogue performance by cohort", ["Cohort", "Drug", "n (R/S)", "Sensitivity (95% CI)", "Specificity (95% CI)", "FSR (95% CI)", "BA"], shown))

    # S6: validation-design benchmark.
    rows = read_csv(root, "results/lineage_validation/leakage_benchmark_summary.csv")
    for antibiotic in ("clarithromycin", "levofloxacin"):
        shown = [[model_name(r["model"]), label(r["split_type"]), r["folds"], percent(r["sensitivity_mean"]), percent(r["specificity_mean"]), percent(r["false_susceptible_rate_mean"]), percent(r["balanced_accuracy_mean"])] for r in rows if r["antibiotic"] == antibiotic]
        suffix = "a" if antibiotic == "clarithromycin" else "b"
        blocks.append(md_table(f"Table S6{suffix}. Validation-design benchmark: {drug(antibiotic)}", ["Model", "Validation design", "Folds", "Sensitivity", "Specificity", "FSR", "BA"], shown))

    # S7: phenotype and breakpoint sensitivity scenarios.
    rows = [r for r in read_csv(root, "results/sensitivity/phenotype_sensitivity_metrics.csv") if r.get("stratum_type") == "COHORT"]
    for antibiotic in ("clarithromycin", "levofloxacin"):
        shown = [[label(r["scenario"]), cohort(r["dataset_id"]), r["n"], percent(r["sensitivity"]), percent(r["specificity"]), percent(r["false_susceptible_rate"]), percent(r["balanced_accuracy"])] for r in rows if r["antibiotic"] == antibiotic]
        suffix = "a" if antibiotic == "clarithromycin" else "b"
        blocks.append(md_table(f"Table S7{suffix}. Phenotype and breakpoint sensitivity: {drug(antibiotic)}", ["Scenario", "Cohort", "n", "Sensitivity", "Specificity", "FSR", "BA"], shown))

    # S8: negative-control inventory and key outputs.
    status = read_csv(root, "results/negative_controls/negative_control_status.csv")
    shown = [[r["control_id"], label(r["control"]), label(r["status"]), r["reason"]] for r in status]
    blocks.append(md_table("Table S8a. Negative-control completion inventory", ["ID", "Control", "Status", "Interpretation/implementation"], shown))
    clone = read_csv(root, "results/negative_controls/clone_thinned_panel_metrics.csv")
    shown = [[cohort(r["dataset_id"]), drug(r["antibiotic"]), r["n"], percent(r["sensitivity"]), percent(r["specificity"]), percent(r["false_susceptible_rate"]), percent(r["balanced_accuracy"])] for r in clone]
    blocks.append(md_table("Table S8b. Near-clone-thinned performance", ["Cohort", "Drug", "n", "Sensitivity", "Specificity", "FSR", "BA"], shown))
    concordance = read_csv(root, "results/negative_controls/raw_assembly_concordance_summary.csv")
    shown = [[drug(r["antibiotic"]), r["n_raw_callable"], r["n_support_assembly_callable"], r["n_both_callable"], r["n_concordant"], r["n_discordant"], percent(r["concordance_when_both_callable"])] for r in concordance]
    blocks.append(md_table("Table S8c. Raw-read versus support-assembly concordance", ["Drug", "Raw callable", "Assembly callable", "Both callable", "Concordant", "Discordant", "Agreement"], shown))
    control_groups: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for path, name in (("results/negative_controls/random_snp_panel_metrics.csv", "Random SNP panel"), ("results/negative_controls/label_permutation_metrics.csv", "Label permutation")):
        for row in read_csv(root, path):
            try:
                value = float(row["balanced_accuracy"])
            except (TypeError, ValueError):
                continue
            if math.isfinite(value):
                control_groups[(name, row["antibiotic"], row["split_type"])].append(value)
    shown = [[name, drug(antibiotic), label(split), len(values), interval(values)] for (name, antibiotic, split), values in sorted(control_groups.items())]
    blocks.append(md_table("Table S8d. Random-panel and label-permutation balanced-accuracy distributions", ["Control", "Drug", "Design", "Replicates", "Median (2.5th-97.5th percentile)"], shown))

    # S9: MIC lower-bound summaries.
    rows = read_csv(root, "results/mic/ningxia_mic_lower_bound_summary.csv")
    shown = [[drug(r["antibiotic"]), r["marker_prediction"], r["n"], r["n_right_censored"], number(r["mic_median_lower_bound"]), f"{number(r['mic_q1_lower_bound'])}-{number(r['mic_q3_lower_bound'])}", pvalue(r["mann_whitney_p"])] for r in rows]
    blocks.append(md_table("Table S9. Ningxia MIC lower-bound summaries", ["Drug", "Marker prediction", "n", "Right-censored", "Median (mg/L)", "IQR (mg/L)", "Mann-Whitney P"], shown))

    # S10: source-genotype and independent-caller audit.
    rows = read_csv(root, "results/external_validation/ningxia_published_genotype_audit_summary.csv")
    shown = []
    for r in rows:
        if r["row_type"] == "PHENOTYPE_PERFORMANCE":
            result = f"Sens {percent(r['sensitivity'])}; spec {percent(r['specificity'])}; FSR {percent(r['false_susceptible_rate'])}"
        else:
            result = f"{r['n_agree']}/{r['n']} agree ({percent(r['call_concordance'])})"
        shown.append([label(r["row_type"]), drug_short(r["antibiotic"]), label(r["estimate_source"]), label(r["scope"]), r["n"], result])
    blocks.append(md_table("Table S10. Ningxia source-genotype and independent-caller audit", ["Audit", "Drug", "Estimate/call source", "Scope", "n", "Result"], shown))

    # S11: frozen classification.
    rows = read_csv(root, "results/external_validation/transportability_classification.csv")
    shown = [[drug(r["antibiotic"]), label(r["transportability_label"]), r["reason"], r["geographic_scope"]] for r in rows]
    blocks.append(md_table("Table S11. Frozen transportability classification", ["Drug", "Classification", "Frozen-rule reason", "Scope"], shown))

    # S13: marker architecture and prevalence shifts.
    rows = read_csv(root, "results/extended_analysis/marker_phenotype_prevalence_shift.csv")
    shown = [[cohort(r["dataset_id"]), drug(r["antibiotic"]), r["n"], percent(r["phenotypic_resistance_prevalence"]), percent(r["marker_resistance_prevalence"]), f"{100 * float(r['prevalence_gap_marker_minus_phenotype']):+.1f} pp"] for r in rows]
    blocks.append(md_table("Table S13a. Phenotype-marker prevalence shift", ["Cohort", "Drug", "n", "Phenotypic R", "Marker R", "Marker minus phenotype"], shown))
    rows = read_csv(root, "results/extended_analysis/marker_prevalence_by_lineage.csv")
    shown = [[cohort(r["dataset_id"]), drug_short(r["antibiotic"]), r["lineage_recomputed"].replace("SNP_CLUSTER_", "C"), r["n"], percent(r["phenotype_resistance_prevalence"]), percent(r["marker_resistance_prevalence"]), percent(r["error_rate"])] for r in rows]
    blocks.append(md_table("Table S13b. Marker prevalence and error by fixed SNP cluster", ["Cohort", "Drug", "Cluster", "n", "Phenotypic R", "Marker R", "Error rate"], shown, "C01-C08 correspond to the fixed labels SNP_CLUSTER_01-SNP_CLUSTER_08."))

    # S14: QC-callability associations.
    rows = read_csv(root, "results/extended_analysis/qc_callability_associations.csv")
    shown = [[cohort(r["dataset_id"]), drug_short(r["antibiotic"]), label(r["metric"]), f"{r['n_callable']}/{clean(r['n_uncallable'])}", number(r["callable_median"], 3), number(r["uncallable_median"], 3), pvalue(r["mann_whitney_p"]), number(r["cliffs_delta_callable_minus_uncallable"], 3)] for r in rows]
    blocks.append(md_table("Table S14. Assembly-quality associations with target callability", ["Cohort", "Drug", "Metric", "Callable/uncallable", "Median callable", "Median uncallable", "P", "Cliff's delta"], shown))

    # S15: manifold and error-distance analyses.
    rows = read_csv(root, "results/extended_analysis/development_manifold_summary.csv")
    shown = [[cohort(r["dataset_id"]), r["n"], number(r["median_5nn_distance"], 3), number(r["median_development_percentile"], 3), pvalue(r["mann_whitney_vs_hpgp_p"]), number(r["cliffs_delta_vs_hpgp"], 3), number(r["lineage_jensen_shannon_distance_vs_hpgp"], 3)] for r in rows]
    blocks.append(md_table("Table S15a. Development-manifold and lineage-composition shift", ["Cohort", "n", "Median 5-NN distance", "Median HpGP percentile", "P vs HpGP", "Cliff's delta", "Lineage JS distance"], shown))
    rows = read_csv(root, "results/extended_analysis/external_error_manifold_associations.csv")
    shown = [[cohort(r["dataset_id"]), drug_short(r["antibiotic"]), f"{r['n_correct']}/{r['n_error']}", number(r["median_correct"], 3), number(r["median_error"], 3), pvalue(r["mann_whitney_p"]), number(r["cliffs_delta_error_minus_correct"], 3)] for r in rows]
    blocks.append(md_table("Table S15b. Development-manifold distance by correct versus incorrect call", ["Cohort", "Drug", "Correct/error n", "Median correct", "Median error", "P", "Cliff's delta"], shown))

    # S16: residual mechanism audit.
    rows = read_csv(root, "results/extended_analysis/ningxia_lvx_false_susceptible_decomposition.csv")
    shown = [[label(r["false_susceptible_mechanism"]), r["n"], percent(r["proportion"])] for r in rows]
    blocks.append(md_table("Table S16a. Ningxia levofloxacin false-susceptible decomposition", ["Mechanism category", "n", "Proportion"], shown))
    rows = read_csv(root, "results/extended_analysis/ningxia_gyrA_missense_associations.csv")
    def assoc_key(row: dict[str, str]) -> tuple[float, float]:
        try:
            return (float(row["bh_q"]), float(row["fisher_p"]))
        except (TypeError, ValueError):
            return (float("inf"), float("inf"))
    shown = [[r["variant"], r["frozen_catalogue_marker"], f"{r['resistant_with']}/{r['resistant_without']}", f"{r['susceptible_with']}/{r['susceptible_without']}", f"{r['false_susceptible_with']}/{r['true_positive_with']}", pvalue(r["fisher_p"]), pvalue(r["bh_q"])] for r in sorted(rows, key=assoc_key)[:12]]
    blocks.append(md_table("Table S16b. Twelve smallest adjusted P values in the full-gyrA missense residual scan", ["Variant", "Frozen marker", "R with/without", "S with/without", "FN/TP with", "Fisher P", "BH q"], shown, "All 140 tested substitution rows are retained in workbook sheet S16_Residual_Mechanisms; none passed false-discovery correction."))
    rows = read_csv(root, "results/extended_analysis/exploratory_N87Y_rescue.csv")
    shown = [[label(r["panel"]), r["n"], f"{r['tp']}/{r['tn']}/{r['fp']}/{r['fn']}", percent(r["sensitivity"]), percent(r["specificity"]), percent(r["false_susceptible_rate"]), percent(r["balanced_accuracy"])] for r in rows]
    blocks.append(md_table("Table S16c. Exploratory add-only N87Y scenario", ["Panel", "n", "TP/TN/FP/FN", "Sensitivity", "Specificity", "FSR", "BA"], shown))

    # S17: stability, predictive-value and decision-grid summaries.
    rows = read_csv(root, "results/extended_analysis/external_bootstrap_distributions.csv")
    grouped_metrics: dict[tuple[str, str], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for r in rows:
        for metric in ("sensitivity", "specificity", "false_susceptible_rate", "balanced_accuracy"):
            try:
                value = float(r[metric])
            except (TypeError, ValueError):
                continue
            if math.isfinite(value):
                grouped_metrics[(r["dataset_id"], r["antibiotic"])][metric].append(value)
    shown = [[cohort(ds), drug_short(ab), len(metrics["balanced_accuracy"]), interval(metrics["sensitivity"]), interval(metrics["specificity"]), interval(metrics["false_susceptible_rate"]), interval(metrics["balanced_accuracy"])] for (ds, ab), metrics in sorted(grouped_metrics.items())]
    blocks.append(md_table("Table S17a. Near-clone-group bootstrap distributions", ["Cohort", "Drug", "Replicates", "Sensitivity", "Specificity", "FSR", "BA"], shown, "Cells report median (2.5th-97.5th percentile)."))
    rows = read_csv(root, "results/extended_analysis/ningxia_lvx_leave_one_out_influence.csv")
    leave_one_out = [r for r in rows if r["dropped_isolate"] != "NONE_BASELINE"]
    shown = [["Baseline", percent(rows[0]["false_susceptible_rate"]), percent(rows[0]["balanced_accuracy"])], ["Leave-one-out range", f"{percent(min(float(r['false_susceptible_rate']) for r in leave_one_out))}-{percent(max(float(r['false_susceptible_rate']) for r in leave_one_out))}", f"{percent(min(float(r['balanced_accuracy']) for r in leave_one_out))}-{percent(max(float(r['balanced_accuracy']) for r in leave_one_out))}"]]
    blocks.append(md_table("Table S17b. Ningxia levofloxacin leave-one-out influence", ["Analysis", "FSR", "BA"], shown))
    rows = read_csv(root, "results/extended_analysis/predictive_values_by_assumed_prevalence.csv")
    selected = []
    for r in rows:
        prevalence = float(r["assumed_prevalence"])
        if any(abs(prevalence - target) < 1e-9 for target in (0.10, 0.25, 0.50)):
            selected.append([cohort(r["dataset_id"]), drug_short(r["antibiotic"]), percent(prevalence), percent(r["ppv"]), percent(r["npv"])])
    blocks.append(md_table("Table S17c. Predictive values under selected assumed prevalences", ["Cohort", "Drug", "Assumed resistance prevalence", "PPV", "NPV"], selected, "These are mathematical prevalence scenarios, not empirical recalibration."))
    rows = read_csv(root, "results/extended_analysis/transport_gate_robustness_grid.csv")
    counts = Counter((r["antibiotic"], r["classification"]) for r in rows)
    shown = [[drug(ab), label(classification), count] for (ab, classification), count in sorted(counts.items())]
    blocks.append(md_table("Table S17d. Transport-gate robustness-grid classifications", ["Drug", "Classification", "Grid cells"], shown, "The prespecified rule remained the primary classification; this grid is a post-freeze robustness display."))

    # S18: four-domain transport-shift decomposition.
    rows = read_csv(root, "results/transport_shift/three_layer_transport_summary.csv")
    shown = [[label(r["layer"]), cohort(r["dataset_id"]), drug_short(r["antibiotic"]) if r["antibiotic"] != "not_applicable" else "-", label(r["measure"]), percent(r["estimate"]) if float(r["estimate"]) <= 1 else number(r["estimate"]), r["detail"]] for r in rows]
    blocks.append(md_table("Table S18a. Four-domain transport-shift summary", ["Domain", "Cohort", "Drug", "Measure", "Estimate", "Detail"], shown))
    rows = read_csv(root, "results/transport_shift/dominant_lineage_comparisons.csv")
    shown = [[label(r["metric"]), cohort(r["dataset_a"]), cohort(r["dataset_b"]), f"{r['dataset_a_successes']}/{r['dataset_a_failures']}", f"{r['dataset_b_successes']}/{r['dataset_b_failures']}", number(r["odds_ratio"], 3), pvalue(r["fisher_p"])] for r in rows]
    blocks.append(md_table("Table S18b. Dominant-lineage conditional comparisons", ["Metric", "Dataset A", "Dataset B", "A success/fail", "B success/fail", "Odds ratio", "Fisher P"], shown))
    rows = read_csv(root, "results/transport_shift/cohort_discriminator_summary.csv")
    shown = [[cohort(r["external_dataset"]), f"{r['n_hpGP']}/{r['n_external']}", r["features"], percent(r["pooled_oof_auc"]), f"{percent(r['group_bootstrap_ci_low'])}-{percent(r['group_bootstrap_ci_high'])}", pvalue(r["permutation_p"])] for r in rows]
    blocks.append(md_table("Table S18c. Phenotype-blind cohort-discriminator audit", ["External cohort", "HpGP/external n", "Features", "OOF AUC", "Group-bootstrap 95% CI", "Permutation P"], shown))
    rows = read_csv(root, "results/transport_shift/ningxia_resistant_mic_summary.csv")
    shown = [[r["outcome"], r["n"], number(r["median_lower_bound_mic_mg_L"]), f"{number(r['q1_lower_bound_mic_mg_L'])}-{number(r['q3_lower_bound_mic_mg_L'])}", r["right_censored_n"]] for r in rows]
    blocks.append(md_table("Table S18d. Ningxia levofloxacin resistant-MIC severity", ["Outcome", "n", "Median lower-bound MIC (mg/L)", "IQR (mg/L)", "Right-censored"], shown))

    # S19: 23S target-recovery sensitivity.
    rows = read_csv(root, "results/extended_analysis/23s_callability_sensitivity_summary.csv")
    shown = [[label(r["section"]), label(r["metric"]), r["value"], r["denominator"], r["notes"]] for r in rows]
    blocks.append(md_table("Table S19a. Frozen 23S configuration and recovery-class summary", ["Section", "Metric", "Value", "Denominator", "Note"], shown))
    grid = read_csv(root, "results/extended_analysis/23s_callability_sensitivity_grid.csv")
    grouped = defaultdict(list)
    for r in grid:
        grouped[r["blast_task"]].append(r)
    shown = []
    for task, items in sorted(grouped.items()):
        identities = [float(r["minimum_identity"]) for r in items]
        coverages = [float(r["minimum_query_coverage"]) for r in items]
        callable_counts = sorted({int(float(r["n_callable_final_qc"])) for r in items})
        shown.append([task, len(items), f"{percent(min(identities))}-{percent(max(identities))}", f"{percent(min(coverages))}-{percent(max(coverages))}", ", ".join(map(str, callable_counts))])
    blocks.append(md_table("Table S19b. 23S callability sensitivity-grid audit", ["BLAST task", "Grid cells", "Identity range", "Coverage range", "Final-QC callable counts"], shown, "Every grid setting required both resistance-marker bases to be spanned; no threshold relaxation rescued an additional final-QC assembly."))

    # S20: coverage-aware and selective-use stress tests.
    rows = read_csv(root, "results/extended_analysis/coverage_aware_yield.csv")
    shown = [[
        cohort(r["dataset_id"]),
        drug_short(r["antibiotic"]),
        r["phenotype_linked_binary_n"],
        r["qc_pass_n"],
        r["primary_evaluable_n"],
        r["correct_n"],
        r["unresolved_n"],
        percent(r["actionable_correct_yield"]),
    ] for r in rows]
    blocks.append(md_table(
        "Table S20a. End-to-end correct-result yield",
        ["Cohort", "Drug", "Binary phenotype", "QC pass", "Primary evaluable", "Correct", "Unresolved", "Correct yield"],
        shown,
        "Correct yield uses every phenotype-linked binary isolate as the denominator; it is not accuracy conditional on receiving a result.",
    ))

    rows = [r for r in read_csv(root, "results/extended_analysis/callability_identification_bounds.csv") if r["metric"] in {"sensitivity", "specificity"}]
    shown = [[
        cohort(r["dataset_id"]),
        drug_short(r["antibiotic"]),
        label(r["metric"]),
        f"{r['resolved_n']}/{r['phenotype_linked_denominator']}",
        percent(r["evaluable_estimate"]),
        f"{percent(r['logical_lower_bound'])}-{percent(r['logical_upper_bound'])}",
    ] for r in rows]
    blocks.append(md_table(
        "Table S20b. Callability-aware logical performance bounds",
        ["Cohort", "Drug", "Metric", "Resolved/total", "Evaluable estimate", "Logical lower-upper bound"],
        shown,
        "Bounds assign unresolved isolates to the least or most favourable compatible prediction; they are not confidence intervals.",
    ))

    selected_cutoffs = {0.05, 0.10, 0.20, 0.50, 0.90, 1.00}
    rows = []
    for r in read_csv(root, "results/extended_analysis/manifold_abstention_metrics.csv"):
        try:
            cutoff = float(r["development_percentile_cutoff"])
        except (TypeError, ValueError):
            continue
        if r["dataset_id"] == "HPGP_GLOBAL" or r["antibiotic"] != "levofloxacin":
            continue
        if not any(abs(cutoff - value) < 1e-9 for value in selected_cutoffs):
            continue
        rows.append([cohort(r["dataset_id"]), percent(cutoff), f"{r['accepted_n']}/{r['original_primary_n']}", percent(r["retained_coverage"]), percent(r["sensitivity"]), percent(r["specificity"]), percent(r["false_susceptible_rate"]), r["passes_frozen_safety_gate"]])
    blocks.append(md_table(
        "Table S20c. Phenotype-blind development-manifold abstention for levofloxacin",
        ["Cohort", "HpGP distance cutoff", "Accepted/primary", "Coverage", "Sensitivity", "Specificity", "FSR", "Passes gate"],
        rows,
        "Cutoffs were fixed HpGP empirical distance percentiles and did not use external phenotypes; a subset also required at least 10 resistant and 10 susceptible isolates to pass the frozen gate.",
    ))

    rows = read_csv(root, "results/extended_analysis/external_performance_differences.csv")
    shown = [[
        label(r["metric"]),
        percent(r["ningxia_estimate"]),
        percent(r["read_cohort_estimate"]),
        f"{100 * float(r['absolute_difference']):+.1f} pp",
        f"{100 * float(r['difference_ci_low']):+.1f} to {100 * float(r['difference_ci_high']):+.1f} pp",
        pvalue(r["fisher_exact_p"]),
    ] for r in rows]
    blocks.append(md_table(
        "Table S20d. External-cohort levofloxacin performance differences",
        ["Metric", "Ningxia", "Read cohort", "Absolute difference", "95% interval", "Fisher P"],
        shown,
        "Differences are Ningxia minus read cohort. Balanced-accuracy uncertainty uses the near-clone-group bootstrap and has no single Fisher exact P value.",
    ))

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n\n".join(blocks).rstrip() + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output or root / "manuscript/supplementary_tables_generated.md"
    build(root, output)


if __name__ == "__main__":
    main()
