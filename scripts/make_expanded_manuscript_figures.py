#!/usr/bin/env python3
"""Render the expanded eight-figure manuscript and supplementary figure set."""

from __future__ import annotations

import argparse
import csv
import math
import textwrap
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

from make_manuscript_figures import (
    COHORT_COLORS,
    COHORT_LABELS,
    COHORT_ORDER,
    DRUG_LABELS,
    numeric,
    panel_label,
    read_csv,
    save_figure,
    style_axis,
    write_source,
)


OUTCOME_COLORS = {"TP": "#2A9D8F", "TN": "#4C78A8", "FP": "#E9C46A", "FN": "#D64B4B"}
SHORT_COHORT = {"HPGP_GLOBAL": "HpGP", "CHINA_NINGXIA_2022": "Ningxia", "ZENODO_10369064": "Reads"}
SHORT_DRUG = {"clarithromycin": "CLR", "levofloxacin": "LVX"}
EXTERNAL_COHORTS = ["CHINA_NINGXIA_2022", "ZENODO_10369064"]


def outcome(phenotype: str, prediction: str) -> str:
    return {
        ("R", "R"): "TP", ("S", "S"): "TN",
        ("S", "R"): "FP", ("R", "S"): "FN",
    }.get((phenotype, prediction), "NE")


def cohort_color(dataset: str) -> str:
    return COHORT_COLORS.get(dataset, "#666666")


def figure1(root: Path, figures: Path, source: Path) -> None:
    phenotypes = read_csv(root / "metadata/phenotype_manifest.csv")
    qc = numeric(read_csv(root / "results/qc/assembly_qc_with_checkm2.csv"), [
        "contigs", "n50_bp", "completeness_percent", "contamination_percent",
    ])
    callability = numeric(read_csv(root / "results/qc/panel_callability.csv"), [
        "n_phenotype_linked", "n_callable", "n_uncallable", "callability",
    ])
    qc_call = numeric(read_csv(root / "results/extended_analysis/qc_callability_samples.csv"), [
        "contigs", "n50_bp", "log10_n50_bp", "callable_binary", "completeness_percent", "contamination_percent",
    ])

    phenotype_counts = phenotypes.groupby("dataset_id")["isolate_id"].nunique().to_dict()
    materialized = qc.groupby("dataset_id")["isolate_id"].nunique().to_dict()
    passed = qc[qc["final_qc_status"] == "PASS"].groupby("dataset_id")["isolate_id"].nunique().to_dict()
    flow = pd.DataFrame([
        {
            "dataset_id": dataset,
            "phenotype_linked": int(phenotype_counts.get(dataset, 0)),
            "sequence_materialized": int(materialized.get(dataset, 0)),
            "final_qc_pass": int(passed.get(dataset, 0)),
            "sequence_source": "paired reads + SKESA" if dataset == "ZENODO_10369064" else "public assembly",
        }
        for dataset in COHORT_ORDER
    ])
    pheno = (
        phenotypes[phenotypes["antibiotic"].isin(DRUG_LABELS)]
        .groupby(["dataset_id", "antibiotic", "susceptibility_original"], dropna=False)
        .size().rename("n").reset_index()
    )
    attrition = flow.assign(final_qc_fail=lambda frame: frame.sequence_materialized - frame.final_qc_pass)
    write_source(flow, source, "Figure1A_expanded_study_flow.csv")
    write_source(pheno, source, "Figure1B_expanded_phenotypes.csv")
    write_source(attrition, source, "Figure1C_expanded_qc_attrition.csv")
    write_source(callability, source, "Figure1D_expanded_callability.csv")
    write_source(qc_call, source, "Figure1EF_expanded_qc_callability_samples.csv")
    write_source(qc, source, "Figure1G_expanded_checkm2.csv")

    fig = plt.figure(figsize=(13.5, 10.3), constrained_layout=True)
    fig._compact_panel_labels = True
    grid = fig.add_gridspec(3, 3, height_ratios=[1.18, 1, 1])
    ax_flow = fig.add_subplot(grid[0, :])
    ax_pheno = fig.add_subplot(grid[1, 0])
    ax_attr = fig.add_subplot(grid[1, 1])
    ax_call = fig.add_subplot(grid[1, 2])
    ax_clr = fig.add_subplot(grid[2, 0])
    ax_lvx = fig.add_subplot(grid[2, 1])
    ax_check = fig.add_subplot(grid[2, 2])

    ax_flow.set_xlim(0, 12)
    ax_flow.set_ylim(0, 4)
    ax_flow.axis("off")
    for index, row in flow.iterrows():
        y = 3.15 - index * 1.22
        color = cohort_color(row.dataset_id)
        boxes = [
            FancyBboxPatch((0.15, y - 0.38), 2.65, 0.76, boxstyle="round,pad=0.03", facecolor=color, edgecolor="none", alpha=0.94),
            FancyBboxPatch((4.0, y - 0.38), 3.0, 0.76, boxstyle="round,pad=0.03", facecolor="#F1F3F5", edgecolor="#73777A"),
            FancyBboxPatch((8.2, y - 0.38), 3.25, 0.76, boxstyle="round,pad=0.03", facecolor="white", edgecolor=color, linewidth=1.8),
        ]
        for box in boxes:
            ax_flow.add_patch(box)
        ax_flow.annotate("", xy=(3.95, y), xytext=(2.85, y), arrowprops={"arrowstyle": "->", "color": "#4D4D4D"})
        ax_flow.annotate("", xy=(8.15, y), xytext=(7.05, y), arrowprops={"arrowstyle": "->", "color": "#4D4D4D"})
        ax_flow.text(1.48, y, f"{COHORT_LABELS[row.dataset_id]}\nn={row.phenotype_linked}", ha="center", va="center", color="white", fontsize=9, fontweight="bold")
        ax_flow.text(5.5, y, f"{row.sequence_source}\nmaterialized n={row.sequence_materialized}", ha="center", va="center", fontsize=8.5)
        ax_flow.text(9.82, y, f"Frozen QC pass n={row.final_qc_pass}\ntarget callability retained", ha="center", va="center", fontsize=8.5)
    panel_label(ax_flow, "A")

    pivot = pheno.pivot_table(index=["dataset_id", "antibiotic"], columns="susceptibility_original", values="n", fill_value=0)
    pivot = pivot.reindex(pd.MultiIndex.from_product([COHORT_ORDER, DRUG_LABELS]))
    x = np.arange(len(pivot))
    bottom = np.zeros(len(pivot))
    for label, color in [("S", "#4C78A8"), ("I", "#F2CF5B"), ("R", "#D64B4B")]:
        values = pivot[label].to_numpy() if label in pivot else np.zeros(len(pivot))
        ax_pheno.bar(x, values, bottom=bottom, color=color, label=label)
        bottom += values
    ax_pheno.set_xticks(x, [f"{SHORT_COHORT[d]}\n{SHORT_DRUG[a]}" for d, a in pivot.index], fontsize=7)
    ax_pheno.set_ylabel("Phenotype-linked isolates")
    ax_pheno.legend(frameon=False, fontsize=7, ncol=3)
    style_axis(ax_pheno); panel_label(ax_pheno, "B")

    attr = attrition.set_index("dataset_id").reindex(COHORT_ORDER)
    ax_attr.bar(np.arange(3), attr["final_qc_pass"], color=[cohort_color(x) for x in COHORT_ORDER], label="QC pass")
    ax_attr.bar(np.arange(3), attr["final_qc_fail"], bottom=attr["final_qc_pass"], color="#D9D9D9", edgecolor="#777777", label="QC fail")
    for i, row in enumerate(attr.itertuples()):
        ax_attr.text(i, max(4, row.final_qc_pass - max(5, 0.04 * row.final_qc_pass)), f"{int(row.final_qc_pass)}/{int(row.sequence_materialized)}", ha="center", va="top", fontsize=7, color="white", fontweight="bold")
    ax_attr.set_xticks(np.arange(3), [SHORT_COHORT[x] for x in COHORT_ORDER])
    ax_attr.set_ylabel("Materialized genomes")
    ax_attr.legend(frameon=False, fontsize=7)
    style_axis(ax_attr); panel_label(ax_attr, "C")

    callability["percent"] = 100 * callability["callability"]
    for j, drug in enumerate(DRUG_LABELS):
        part = callability[callability["antibiotic"] == drug].set_index("dataset_id").reindex(COHORT_ORDER)
        ax_call.bar(np.arange(3) + (j - 0.5) * 0.34, part["percent"], width=0.32, label=SHORT_DRUG[drug])
    ax_call.axhline(90, color="#777777", linestyle="--", linewidth=1)
    ax_call.set_ylim(0, 105)
    ax_call.set_xticks(np.arange(3), [SHORT_COHORT[x] for x in COHORT_ORDER])
    ax_call.set_ylabel("Callable predictions (%)")
    ax_call.legend(frameon=False, fontsize=7)
    style_axis(ax_call); panel_label(ax_call, "D")

    for axis, drug, letter in [(ax_clr, "clarithromycin", "E"), (ax_lvx, "levofloxacin", "F")]:
        part = qc_call[(qc_call["antibiotic"] == drug) & qc_call["dataset_id"].isin(EXTERNAL_COHORTS)]
        for dataset in EXTERNAL_COHORTS:
            cohort = part[part["dataset_id"] == dataset]
            callable_rows = cohort[cohort["callable_binary"] == 1]
            uncallable_rows = cohort[cohort["callable_binary"] == 0]
            axis.scatter(callable_rows["contigs"], callable_rows["n50_bp"], s=24, color=cohort_color(dataset), alpha=0.72, label=f"{SHORT_COHORT[dataset]} callable")
            axis.scatter(uncallable_rows["contigs"], uncallable_rows["n50_bp"], s=30, facecolor="none", edgecolor=cohort_color(dataset), linewidth=1.1, alpha=0.85, label=f"{SHORT_COHORT[dataset]} uncallable")
        axis.set_yscale("log")
        axis.set_xscale("log")
        axis.set_xlabel("Contigs")
        axis.set_ylabel("Assembly N50 (bp)")
        style_axis(axis); panel_label(axis, letter)
    handles, labels = ax_clr.get_legend_handles_labels()
    ax_clr.legend(handles[:4], labels[:4], frameon=False, fontsize=6, loc="lower left")

    for dataset in COHORT_ORDER:
        part = qc[qc["dataset_id"] == dataset]
        ax_check.scatter(part["contamination_percent"], part["completeness_percent"], s=18, alpha=0.65, color=cohort_color(dataset), label=SHORT_COHORT[dataset])
    ax_check.axhline(90, color="#777777", linestyle="--", linewidth=1)
    ax_check.axvline(5, color="#777777", linestyle="--", linewidth=1)
    ax_check.set_xscale("symlog", linthresh=1)
    ax_check.set_xlim(0, 70)
    ax_check.set_xlabel("CheckM2 contamination (%)")
    ax_check.set_ylabel("CheckM2 completeness (%)")
    ax_check.legend(frameon=False, fontsize=7)
    style_axis(ax_check); panel_label(ax_check, "G")
    save_figure(fig, figures, "Figure1_study_design_data_lineage")


def figure2(root: Path, figures: Path, source: Path) -> None:
    lineages = numeric(read_csv(root / "results/extended_analysis/development_manifold_distance.csv"), [
        "PC1", "PC2", "development_5nn_distance", "development_5nn_percentile",
    ])
    summary = numeric(read_csv(root / "results/extended_analysis/development_manifold_summary.csv"), [
        "lineage_jensen_shannon_distance_vs_hpgp", "median_5nn_distance",
    ])
    errors = numeric(read_csv(root / "results/extended_analysis/external_error_manifold_samples.csv"), ["development_5nn_distance"])
    groups = read_csv(root / "results/qc/near_clone_groups.csv")
    write_source(lineages, source, "Figure2A-C_expanded_manifold.csv")
    write_source(summary, source, "Figure2D_expanded_divergence.csv")
    write_source(errors, source, "Figure2E_expanded_error_distance.csv")
    write_source(groups, source, "Figure2F_expanded_near_clones.csv")

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.2), constrained_layout=True)
    fig._compact_panel_labels = True
    ax = axes[0, 0]
    for dataset in COHORT_ORDER:
        part = lineages[lineages["dataset_id"] == dataset]
        ax.scatter(part["PC1"], part["PC2"], s=17, alpha=0.7, color=cohort_color(dataset), label=SHORT_COHORT[dataset], edgecolors="none")
    ax.set_xlabel("Development-fitted component 1"); ax.set_ylabel("Component 2")
    ax.legend(frameon=False, fontsize=7)
    style_axis(ax); panel_label(ax, "A")

    mix = pd.crosstab(lineages["lineage_recomputed"], lineages["dataset_id"], normalize="columns").reindex(columns=COHORT_ORDER, fill_value=0)
    image = axes[0, 1].imshow(mix.to_numpy(), aspect="auto", cmap="Blues", vmin=0, vmax=max(0.01, mix.to_numpy().max()))
    axes[0, 1].set_xticks(np.arange(3), [SHORT_COHORT[x] for x in COHORT_ORDER])
    axes[0, 1].set_yticks(np.arange(len(mix)), [x.replace("SNP_CLUSTER_", "C") for x in mix.index], fontsize=7)
    fig.colorbar(image, ax=axes[0, 1], label="Within-cohort proportion", shrink=0.72)
    panel_label(axes[0, 1], "B")

    data = [lineages[lineages["dataset_id"] == dataset]["development_5nn_distance"].dropna().to_numpy() for dataset in COHORT_ORDER]
    violin = axes[0, 2].violinplot(data, positions=np.arange(3), showmedians=True, showextrema=False)
    for body, dataset in zip(violin["bodies"], COHORT_ORDER):
        body.set_facecolor(cohort_color(dataset)); body.set_alpha(0.72)
    violin["cmedians"].set_color("#222222")
    axes[0, 2].set_xticks(np.arange(3), [SHORT_COHORT[x] for x in COHORT_ORDER])
    axes[0, 2].set_yscale("log")
    axes[0, 2].set_ylabel("Mean standardized distance to 5 nearest HpGP genomes")
    style_axis(axes[0, 2]); panel_label(axes[0, 2], "C")

    ext_summary = summary[summary["dataset_id"].isin(EXTERNAL_COHORTS)]
    axes[1, 0].bar(np.arange(2), ext_summary["lineage_jensen_shannon_distance_vs_hpgp"], color=[cohort_color(x) for x in ext_summary["dataset_id"]])
    axes[1, 0].set_xticks(np.arange(2), [SHORT_COHORT[x] for x in ext_summary["dataset_id"]])
    axes[1, 0].set_ylim(0, 1)
    axes[1, 0].set_ylabel("Jensen-Shannon distance vs HpGP")
    style_axis(axes[1, 0]); panel_label(axes[1, 0], "D")

    ningxia = errors[(errors["dataset_id"] == "CHINA_NINGXIA_2022") & (errors["antibiotic"] == "levofloxacin")]
    rng = np.random.default_rng(20260831)
    for i, label in enumerate([False, True]):
        values = ningxia[ningxia["error"].str.lower() == str(label).lower()]["development_5nn_distance"].dropna().to_numpy()
        axes[1, 1].scatter(np.full(len(values), i) + rng.normal(0, 0.045, len(values)), values, s=26, alpha=0.72, color=["#4C78A8", "#D64B4B"][i])
        if len(values):
            axes[1, 1].plot([i - 0.18, i + 0.18], [np.median(values)] * 2, color="black", linewidth=2)
    axes[1, 1].set_xticks([0, 1], ["Correct", "FP/FN"])
    axes[1, 1].set_ylabel("5-nearest-neighbour distance")
    axes[1, 1].text(0.98, 0.98, "P=0.968", transform=axes[1, 1].transAxes, ha="right", va="top", fontsize=8)
    style_axis(axes[1, 1]); panel_label(axes[1, 1], "E")

    sizes = groups.groupby("near_clone_group").size()
    counts = Counter(sizes.to_list())
    x_values = sorted(counts)
    axes[1, 2].bar([str(x) for x in x_values], [counts[x] for x in x_values], color="#6B6ECF")
    cross = groups.groupby("near_clone_group")["dataset_id"].nunique()
    axes[1, 2].text(0.98, 0.95, f"Cross-cohort components: {(cross > 1).sum()}", transform=axes[1, 2].transAxes, ha="right", va="top", fontsize=8)
    axes[1, 2].set_xlabel("Near-clone component size"); axes[1, 2].set_ylabel("Components")
    style_axis(axes[1, 2]); panel_label(axes[1, 2], "F")
    save_figure(fig, figures, "Figure2_population_structure_relatedness")


def forest(axis: plt.Axes, rows: pd.DataFrame, metric: str, gate: float | None = None) -> None:
    y = np.arange(len(rows))[::-1]
    labels = []
    for position, (_, row) in zip(y, rows.iterrows()):
        labels.append(f"{SHORT_COHORT.get(row.dataset_id, row.dataset_id)} (n={int(row.n)})")
        if not bool(row.estimable):
            axis.text(0.02, position, f"NE: {int(row.n_resistant)} R / {int(row.n_susceptible)} S", va="center", fontsize=7, color="#777777")
            continue
        value = row[metric]
        low = row[f"{metric}_ci_low"]
        high = row[f"{metric}_ci_high"]
        axis.errorbar(value, position, xerr=[[value - low], [high - value]], fmt="o", color=cohort_color(row.dataset_id), capsize=3)
    if gate is not None:
        axis.axvline(gate, color="#777777", linestyle="--", linewidth=1)
    axis.set_xlim(0, 1.02)
    axis.set_yticks(y, labels)
    axis.set_xlabel(f"{metric.replace('_', ' ').capitalize()} (95% CI)")
    style_axis(axis)


def figure3(root: Path, figures: Path, source: Path) -> None:
    metrics = numeric(read_csv(root / "results/external_validation/frozen_panel_metrics.csv"), [
        "n", "n_resistant", "n_susceptible", "sensitivity", "specificity", "false_susceptible_rate",
        "balanced_accuracy", "mcc", "sensitivity_ci_low", "sensitivity_ci_high",
        "specificity_ci_low", "specificity_ci_high", "false_susceptible_rate_ci_low",
        "false_susceptible_rate_ci_high", "balanced_accuracy_bootstrap_low", "balanced_accuracy_bootstrap_high",
    ])
    cohorts = metrics[(metrics["stratum_type"] == "COHORT") & metrics["dataset_id"].isin(COHORT_ORDER)].copy()
    cohorts["estimable"] = (cohorts["n_resistant"] >= 10) & (cohorts["n_susceptible"] >= 10)
    cohorts["rank"] = cohorts["dataset_id"].map({x: i for i, x in enumerate(COHORT_ORDER)})
    cohorts = cohorts.sort_values(["antibiotic", "rank"])
    write_source(cohorts.drop(columns="rank"), source, "Figure3_expanded_diagnostic_performance.csv")

    fig, axes = plt.subplots(3, 2, figsize=(10.8, 10.5), constrained_layout=True)
    fig._compact_panel_labels = True
    letters = iter("ABCDEF")
    for row, drug in enumerate(DRUG_LABELS):
        part = cohorts[cohorts["antibiotic"] == drug]
        forest(axes[row, 0], part, "sensitivity", 0.90)
        panel_label(axes[row, 0], next(letters))
        forest(axes[row, 1], part, "specificity", 0.90)
        panel_label(axes[row, 1], next(letters))

    estimable = cohorts[cohorts["estimable"]].copy()
    labels = [f"{SHORT_COHORT[d]} {SHORT_DRUG[a]}" for d, a in zip(estimable.dataset_id, estimable.antibiotic)]
    y = np.arange(len(estimable))[::-1]
    for position, (_, row) in zip(y, estimable.iterrows()):
        low = row["balanced_accuracy_bootstrap_low"]
        high = row["balanced_accuracy_bootstrap_high"]
        axes[2, 0].errorbar(row["balanced_accuracy"], position, xerr=[[row["balanced_accuracy"] - low], [high - row["balanced_accuracy"]]], fmt="o", color=cohort_color(row["dataset_id"]), capsize=3)
    axes[2, 0].axvline(0.90, color="#777777", linestyle="--", linewidth=1)
    axes[2, 0].set_xlim(0, 1.02); axes[2, 0].set_yticks(y, labels); axes[2, 0].set_xlabel("Balanced accuracy (bootstrap 95% interval)")
    style_axis(axes[2, 0]); panel_label(axes[2, 0], next(letters))

    for position, (_, row) in zip(y, estimable.iterrows()):
        value = row["false_susceptible_rate"]
        low = row["false_susceptible_rate_ci_low"]
        high = row["false_susceptible_rate_ci_high"]
        axes[2, 1].errorbar(value, position, xerr=[[value - low], [high - value]], fmt="o", color=cohort_color(row["dataset_id"]), capsize=3)
    axes[2, 1].axvline(0.10, color="#777777", linestyle="--", linewidth=1)
    axes[2, 1].set_xlim(0, 0.86); axes[2, 1].set_yticks(y, labels); axes[2, 1].set_xlabel("False-susceptible rate (exact 95% CI)")
    style_axis(axes[2, 1]); panel_label(axes[2, 1], next(letters))
    save_figure(fig, figures, "Figure3_external_performance_frozen_catalogues")


def figure4(root: Path, figures: Path, source: Path) -> None:
    folds = numeric(read_csv(root / "results/lineage_validation/leakage_benchmark_folds.csv"), ["balanced_accuracy", "auroc_probability"])
    random_panel = numeric(read_csv(root / "results/negative_controls/random_snp_panel_metrics.csv"), ["auroc_probability"])
    permutations = numeric(read_csv(root / "results/negative_controls/label_permutation_metrics.csv"), ["auroc_probability"])
    clone = numeric(read_csv(root / "results/negative_controls/clone_thinned_panel_metrics.csv"), ["balanced_accuracy"])
    concordance = numeric(read_csv(root / "results/negative_controls/raw_assembly_concordance_summary.csv"), [
        "n_raw_callable", "n_support_assembly_callable", "n_both_callable", "n_concordant", "concordance_when_both_callable",
    ])
    write_source(folds, source, "Figure4AB_expanded_validation_design.csv")
    write_source(random_panel, source, "Figure4C_expanded_random_panels.csv")
    write_source(permutations, source, "Figure4D_expanded_label_permutations.csv")
    write_source(clone, source, "Figure4E_expanded_clone_thinning.csv")
    write_source(concordance, source, "Figure4F_expanded_raw_assembly_concordance.csv")

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.3), constrained_layout=True)
    fig._compact_panel_labels = True
    splits = ["RANDOM_ISOLATE_SPLIT", "CLONE_GROUPED_SPLIT", "LEAVE_COUNTRY_OUT", "LEAVE_LINEAGE_OUT", "HPGP_TO_EXTERNAL"]
    models = ["MUTATION_ONLY_LOGISTIC", "LINEAGE_ONLY_LOGISTIC", "MUTATION_PLUS_LINEAGE_LOGISTIC"]
    model_labels = {models[0]: "Mutation", models[1]: "Lineage", models[2]: "Mutation + lineage"}
    colors = {models[0]: "#2878B5", models[1]: "#E39C37", models[2]: "#7A5195"}
    rng = np.random.default_rng(20260831)
    for axis, drug, letter in [(axes[0, 0], "clarithromycin", "A"), (axes[0, 1], "levofloxacin", "B")]:
        part = folds[folds["antibiotic"] == drug]
        for model_index, model in enumerate(models):
            for split_index, split in enumerate(splits):
                values = part[(part["model"] == model) & (part["split_type"] == split)]["balanced_accuracy"].dropna().to_numpy()
                if not len(values):
                    continue
                x = split_index + (model_index - 1) * 0.22
                axis.scatter(x + rng.normal(0, 0.025, len(values)), values, s=11, alpha=0.28, color=colors[model])
                axis.errorbar(x, values.mean(), yerr=values.std(ddof=1) if len(values) > 1 else 0, fmt="D", color=colors[model], capsize=3, markersize=4)
        axis.axhline(0.5, color="#AAAAAA", linestyle=":")
        axis.set_xticks(np.arange(5), ["Random", "Clone", "Country", "Lineage", "External"], rotation=25, ha="right")
        axis.set_ylim(0, 1.02); axis.set_ylabel("Balanced accuracy")
        style_axis(axis); panel_label(axis, letter)
    axes[0, 0].legend(
        handles=[Line2D([0], [0], marker="D", linestyle="", color=colors[m], label=model_labels[m]) for m in models],
        frameon=False,
        fontsize=7,
        ncol=3,
        loc="lower left",
        bbox_to_anchor=(0.0, 1.015),
        borderaxespad=0,
        columnspacing=0.9,
        handletextpad=0.35,
    )

    for axis, frame, title, letter in [
        (axes[0, 2], random_panel, "Size/prevalence-matched random SNP panels", "C"),
        (axes[1, 0], permutations, "Training-label permutations", "D"),
    ]:
        data = [frame[frame["antibiotic"] == drug]["auroc_probability"].dropna().to_numpy() for drug in DRUG_LABELS]
        violin = axis.violinplot(data, positions=[0, 1], showmedians=False, showextrema=False)
        for body, color in zip(violin["bodies"], ["#4C78A8", "#D95F02"]):
            body.set_facecolor(color); body.set_alpha(0.68)
        means = [float(np.mean(values)) for values in data]
        axis.scatter([0, 1], means, marker="D", s=28, color="black", zorder=4)
        for position, value in enumerate(means):
            axis.text(position, min(0.96, value + 0.055), f"mean {value:.2f}", ha="center", fontsize=7)
        axis.axhline(0.5, color="#777777", linestyle="--", linewidth=1)
        axis.set_xticks([0, 1], ["CLR", "LVX"]); axis.set_ylim(0, 1)
        axis.set_ylabel("Probability AUROC")
        style_axis(axis); panel_label(axis, letter)

    baseline = numeric(read_csv(root / "results/external_validation/frozen_panel_metrics.csv"), ["balanced_accuracy"])
    base_cohort = baseline[baseline["stratum_type"] == "COHORT"][["dataset_id", "antibiotic", "balanced_accuracy"]].rename(columns={"balanced_accuracy": "baseline_ba"})
    clone_delta = clone.merge(base_cohort, on=["dataset_id", "antibiotic"], how="left")
    clone_delta["delta"] = clone_delta["balanced_accuracy"] - clone_delta["baseline_ba"]
    clone_delta = clone_delta[clone_delta["dataset_id"].isin(COHORT_ORDER)]
    labels = [f"{SHORT_COHORT[d]} {SHORT_DRUG[a]}" for d, a in zip(clone_delta.dataset_id, clone_delta.antibiotic)]
    for position, row in enumerate(clone_delta.itertuples()):
        axes[1, 1].plot([row.baseline_ba, row.balanced_accuracy], [position, position], color="#999999", linewidth=1)
        axes[1, 1].scatter(row.baseline_ba, position, s=28, facecolor="white", edgecolor=cohort_color(row.dataset_id), linewidth=1.3)
        axes[1, 1].scatter(row.balanced_accuracy, position, s=28, color=cohort_color(row.dataset_id))
    axes[1, 1].set_xlim(0.45, 1.02)
    axes[1, 1].set_yticks(np.arange(len(clone_delta)), labels, fontsize=7)
    axes[1, 1].set_xlabel("Balanced accuracy (open: full; filled: clone-thinned)")
    style_axis(axes[1, 1]); panel_label(axes[1, 1], "E")

    x = np.arange(len(concordance))
    axes[1, 2].bar(x - 0.18, concordance["n_both_callable"], width=0.36, color="#BDBDBD", label="Both callable")
    axes[1, 2].bar(x + 0.18, concordance["n_concordant"], width=0.36, color="#3A9D5D", label="Concordant")
    for i, row in enumerate(concordance.itertuples()):
        axes[1, 2].text(i + 0.18, row.n_concordant + 1, f"{row.concordance_when_both_callable:.0%}", ha="center", fontsize=8)
    axes[1, 2].set_xticks(x, [SHORT_DRUG[x] for x in concordance.antibiotic])
    axes[1, 2].set_ylabel("Raw-read cohort isolates")
    axes[1, 2].set_ylim(0, max(concordance["n_concordant"].max(), concordance["n_both_callable"].max()) + 8)
    axes[1, 2].legend(
        frameon=False,
        fontsize=7,
        ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.015),
        borderaxespad=0,
        columnspacing=1.2,
        handletextpad=0.45,
    )
    style_axis(axes[1, 2]); panel_label(axes[1, 2], "F")
    save_figure(fig, figures, "Figure4_validation_leakage_benchmark")


def confusion_matrix(axis: plt.Axes, frame: pd.DataFrame, row_field: str, col_field: str, letter: str) -> None:
    matrix = pd.crosstab(frame[row_field], frame[col_field]).reindex(index=["R", "S"], columns=["R", "S"], fill_value=0)
    image = axis.imshow(matrix.to_numpy(), cmap="Blues", vmin=0)
    axis.set_xticks([0, 1], ["R", "S"]); axis.set_yticks([0, 1], ["R", "S"])
    axis.set_xlabel(col_field.replace("_", " ")); axis.set_ylabel(row_field.replace("_", " "))
    for i in range(2):
        for j in range(2):
            axis.text(j, i, str(int(matrix.iloc[i, j])), ha="center", va="center", fontsize=12, color="white" if matrix.iloc[i, j] > matrix.to_numpy().max() / 2 else "black")
    panel_label(axis, letter)


def figure5(root: Path, figures: Path, source: Path) -> None:
    audit = numeric(read_csv(root / "results/extended_analysis/ningxia_lvx_error_audit_enriched.csv"), ["mic_numeric"])
    decomposition = numeric(read_csv(root / "results/extended_analysis/ningxia_lvx_false_susceptible_decomposition.csv"), ["n", "proportion"])
    variants = numeric(read_csv(root / "results/extended_analysis/ningxia_gyrA_missense_associations.csv"), ["position", "fisher_p", "bh_q"])
    write_source(audit, source, "Figure5A-E_expanded_Ningxia_error_audit.csv")
    write_source(decomposition, source, "Figure5C_expanded_false_susceptible_decomposition.csv")
    write_source(variants, source, "Figure5F_expanded_gyrA_scan.csv")
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.3), constrained_layout=True)
    fig._compact_panel_labels = True

    confusion_matrix(axes[0, 0], audit, "phenotype", "prediction", "A")
    rng = np.random.default_rng(20260831)
    outcome_order = ["TN", "FP", "TP", "FN"]
    for i, name in enumerate(outcome_order):
        values = audit[(audit["outcome"] == name) & audit["mic_numeric"].notna() & (audit["mic_numeric"] > 0)]["mic_numeric"].to_numpy()
        axes[0, 1].scatter(np.full(len(values), i) + rng.normal(0, 0.055, len(values)), values, s=28, alpha=0.75, color=OUTCOME_COLORS[name], edgecolor="white", linewidth=0.3)
        if len(values):
            axes[0, 1].plot([i - 0.18, i + 0.18], [np.median(values)] * 2, color="black", linewidth=2)
    axes[0, 1].set_yscale("log", base=2); axes[0, 1].set_xticks(np.arange(4), outcome_order)
    axes[0, 1].set_ylabel("MIC lower bound (mg/L)")
    style_axis(axes[0, 1]); panel_label(axes[0, 1], "B")

    labels = ["Deposited QRDR\nwild type", "Off-panel\nN87Y"]
    ordered = decomposition.set_index("false_susceptible_mechanism").reindex(["QRDR_wild_type_in_deposited_assembly", "off_panel_N87Y"])
    bars = axes[0, 2].bar(np.arange(2), ordered["n"], color=["#D64B4B", "#E9C46A"])
    for bar, value in zip(bars, ordered["n"]):
        axes[0, 2].text(bar.get_x() + bar.get_width() / 2, value + 0.3, f"{int(value)}/16", ha="center", fontsize=9)
    axes[0, 2].set_xticks(np.arange(2), labels); axes[0, 2].set_ylim(0, 17)
    axes[0, 2].set_ylabel("False-susceptible isolates")
    style_axis(axes[0, 2]); panel_label(axes[0, 2], "C")

    source_rows = audit[audit["published_genotype_call"].isin(["R", "S"]) & audit["recalled_frozen_panel_call"].isin(["R", "S"])]
    confusion_matrix(axes[1, 0], source_rows, "published_genotype_call", "recalled_frozen_panel_call", "D")
    caller_rows = audit[audit["recalled_frozen_panel_call"].isin(["R", "S"]) & audit["snippy_frozen_panel_call"].isin(["R", "S"])]
    confusion_matrix(axes[1, 1], caller_rows, "recalled_frozen_panel_call", "snippy_frozen_panel_call", "E")

    plot = variants[variants["fisher_p"] > 0].copy()
    plot["minus_log10_p"] = -np.log10(plot["fisher_p"])
    colors = np.where(plot["frozen_catalogue_marker"].str.lower() == "true", "#D95F02", "#4C78A8")
    axes[1, 2].scatter(plot["position"], plot["minus_log10_p"], s=18, color=colors, alpha=0.75)
    axes[1, 2].axhline(-math.log10(0.05), color="#777777", linestyle="--", linewidth=1, label="Nominal P=0.05")
    for _, row in plot.nsmallest(5, "fisher_p").iterrows():
        axes[1, 2].annotate(row["variant"], (row["position"], row["minus_log10_p"]), xytext=(2, 3), textcoords="offset points", fontsize=6)
    axes[1, 2].text(0.02, 0.95, "No FDR-significant variant", transform=axes[1, 2].transAxes, ha="left", va="top", fontsize=8)
    axes[1, 2].set_xlabel("gyrA amino-acid position"); axes[1, 2].set_ylabel("−log10 Fisher P")
    axes[1, 2].legend(handles=[
        Line2D([0], [0], marker="o", linestyle="", color="#D95F02", label="Frozen marker"),
        Line2D([0], [0], marker="o", linestyle="", color="#4C78A8", label="Other missense"),
    ], frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0, 0.88))
    style_axis(axes[1, 2]); panel_label(axes[1, 2], "F")
    save_figure(fig, figures, "Figure5_MIC_false_susceptible_sensitivity")


def gate_grid(axis: plt.Axes, grid: pd.DataFrame, drug: str, letter: str) -> None:
    part = grid[grid["antibiotic"] == drug].copy()
    minima = sorted(part["minimum_resistant_and_susceptible_per_external_cohort"].unique())
    thresholds = sorted(part["false_susceptible_gate"].unique())
    codes = {"INSUFFICIENT_DATA": 0, "HIGH_FALSE_SUSCEPTIBLE_RISK": 1, "PASSES_EXPLORATORY_GRID": 2}
    matrix = np.zeros((len(minima), len(thresholds)))
    for i, minimum in enumerate(minima):
        for j, threshold in enumerate(thresholds):
            label = part[(part["minimum_resistant_and_susceptible_per_external_cohort"] == minimum) & (part["false_susceptible_gate"] == threshold)]["classification"].iloc[0]
            matrix[i, j] = codes[label]
    from matplotlib.colors import ListedColormap
    axis.imshow(matrix, aspect="auto", cmap=ListedColormap(["#BDBDBD", "#D95F02", "#3A9D5D"]), vmin=0, vmax=2)
    axis.set_xticks(np.arange(len(thresholds)), [f"{x:.2f}" for x in thresholds], rotation=45, ha="right", fontsize=7)
    axis.set_yticks(np.arange(len(minima)), [str(int(x)) for x in minima])
    axis.set_xlabel("False-susceptible gate"); axis.set_ylabel("Minimum R and S per external cohort")
    frozen_i = minima.index(10); frozen_j = thresholds.index(0.10)
    axis.scatter(frozen_j, frozen_i, marker="s", s=85, facecolor="none", edgecolor="black", linewidth=1.5)
    panel_label(axis, letter)


def figure6(root: Path, figures: Path, source: Path) -> None:
    classification = read_csv(root / "results/external_validation/transportability_classification.csv")
    metrics = numeric(read_csv(root / "results/external_validation/frozen_panel_metrics.csv"), [
        "n", "n_resistant", "n_susceptible", "sensitivity", "specificity", "false_susceptible_rate",
    ])
    grid_data = numeric(read_csv(root / "results/extended_analysis/transport_gate_robustness_grid.csv"), [
        "minimum_resistant_and_susceptible_per_external_cohort", "false_susceptible_gate", "eligible_external_cohorts",
    ])
    write_source(metrics, source, "Figure6A_expanded_boundary_metrics.csv")
    write_source(classification, source, "Figure6B_expanded_classification.csv")
    write_source(grid_data, source, "Figure6CD_expanded_gate_grid.csv")
    eligible = metrics[metrics["stratum_type"].isin(["COHORT", "LINEAGE"])].copy()
    eligible["row_label"] = np.where(eligible["stratum_type"] == "COHORT", eligible["stratum"].map(lambda x: SHORT_COHORT.get(x, x)), eligible["stratum"].str.replace("SNP_CLUSTER_", "Cluster "))
    eligible["estimable"] = ~(
        (eligible["stratum_type"] == "COHORT")
        & eligible["dataset_id"].isin(EXTERNAL_COHORTS)
        & ((eligible["n_resistant"] < 10) | (eligible["n_susceptible"] < 10))
    )
    fig = plt.figure(figsize=(13.5, max(9.2, 0.25 * len(eligible) + 4.2)), constrained_layout=True)
    fig._compact_panel_labels = True
    layout = fig.add_gridspec(2, 2, height_ratios=[1.55, 1], width_ratios=[1.45, 1])
    ax_heat = fig.add_subplot(layout[0, 0]); ax_text = fig.add_subplot(layout[0, 1])
    ax_clr = fig.add_subplot(layout[1, 0]); ax_lvx = fig.add_subplot(layout[1, 1])
    rows, labels = [], []
    for drug in DRUG_LABELS:
        for _, row in eligible[eligible["antibiotic"] == drug].iterrows():
            rows.append([row[x] for x in ["sensitivity", "specificity", "false_susceptible_rate"]] if row["estimable"] else [np.nan] * 3)
            labels.append(f"{SHORT_DRUG[drug]} | {row.row_label}" + (" [NE]" if not row.estimable else ""))
    matrix = np.asarray(rows, dtype=float)
    scores = np.full_like(matrix, np.nan)
    scores[:, :2] = np.clip(0.5 + (matrix[:, :2] - 0.90) / 0.20, 0, 1)
    scores[:, 2] = np.clip(0.5 + (0.10 - matrix[:, 2]) / 0.20, 0, 1)
    cmap = plt.get_cmap("RdYlGn").copy(); cmap.set_bad("#D9D9D9")
    image = ax_heat.imshow(np.ma.masked_invalid(scores), aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax_heat.set_xticks([0, 1, 2], ["Sensitivity", "Specificity", "False-susceptible\nrate"])
    ax_heat.set_yticks(np.arange(len(labels)), labels, fontsize=7)
    for i in range(matrix.shape[0]):
        for j in range(3):
            ax_heat.text(j, i, "NE" if np.isnan(matrix[i, j]) else f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=6.5)
    cb = fig.colorbar(image, ax=ax_heat, shrink=0.72, label="Frozen gate score")
    cb.set_ticks([0, 0.5, 1], labels=["Fail", "Gate", "Pass"])
    panel_label(ax_heat, "A")

    ax_text.axis("off")
    y = 0.95
    for _, row in classification.iterrows():
        ax_text.text(0.02, y, DRUG_LABELS[row.antibiotic], transform=ax_text.transAxes, fontsize=12, fontweight="bold")
        ax_text.text(0.02, y - 0.07, row.transportability_label.replace("_", " "), transform=ax_text.transAxes, fontsize=10, color="#D95F02", fontweight="bold")
        ax_text.text(0.02, y - 0.14, textwrap.fill(row.reason, 46), transform=ax_text.transAxes, fontsize=8, va="top")
        y -= 0.34
    ax_text.text(0.02, 0.12, "Both external cohorts were Chinese;\nno global external-validation claim.", transform=ax_text.transAxes, fontsize=8, color="#555555")
    ax_text.text(0.02, 0.025, "Gate grids: grey = insufficient data; orange = high risk; green = exploratory pass.", transform=ax_text.transAxes, fontsize=6.5, color="#555555")
    panel_label(ax_text, "B")
    gate_grid(ax_clr, grid_data, "clarithromycin", "C")
    gate_grid(ax_lvx, grid_data, "levofloxacin", "D")
    save_figure(fig, figures, "Figure6_transportability_boundaries")


def figure7(root: Path, figures: Path, source: Path) -> None:
    spectrum = numeric(read_csv(root / "results/extended_analysis/mutation_spectrum_summary.csv"), ["n"])
    combinations = numeric(read_csv(root / "results/extended_analysis/marker_combination_summary.csv"), ["n"])
    prevalence = numeric(read_csv(root / "results/extended_analysis/marker_phenotype_prevalence_shift.csv"), [
        "phenotypic_resistance_prevalence", "marker_resistance_prevalence", "prevalence_gap_marker_minus_phenotype",
    ])
    lineage = numeric(read_csv(root / "results/extended_analysis/marker_prevalence_by_lineage.csv"), [
        "n", "phenotype_resistance_prevalence", "marker_resistance_prevalence", "error_rate",
    ])
    sample_markers = read_csv(root / "results/extended_analysis/mutation_spectrum_samples.csv")
    write_source(spectrum, source, "Figure7AB_mutation_spectrum.csv")
    write_source(combinations, source, "Figure7C_marker_combinations.csv")
    write_source(prevalence, source, "Figure7D_prevalence_shift.csv")
    write_source(lineage, source, "Figure7E_lineage_marker_prevalence.csv")
    write_source(sample_markers, source, "Figure7F_outcome_composition.csv")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.4), constrained_layout=True)
    fig._compact_panel_labels = True

    for axis, drug, letter in [(axes[0, 0], "clarithromycin", "A"), (axes[0, 1], "levofloxacin", "B")]:
        part = spectrum[spectrum["antibiotic"] == drug]
        mutation_order = part.groupby("mutation")["n"].sum().sort_values(ascending=False).index.tolist()[:8]
        counts = part[part["mutation"].isin(mutation_order)].groupby(["mutation", "dataset_id"])["n"].sum().unstack(fill_value=0).reindex(index=mutation_order, columns=COHORT_ORDER, fill_value=0)
        denominators = sample_markers[sample_markers["antibiotic"] == drug].drop_duplicates(["dataset_id", "isolate_id"]).groupby("dataset_id")["isolate_id"].size().reindex(COHORT_ORDER, fill_value=1)
        proportions = counts.divide(denominators, axis=1)
        image = axis.imshow(proportions.to_numpy(), aspect="auto", cmap="YlGnBu", vmin=0, vmax=max(0.01, float(proportions.to_numpy().max())))
        axis.set_xticks(np.arange(3), [SHORT_COHORT[x] for x in COHORT_ORDER])
        axis.set_yticks(np.arange(len(proportions)), proportions.index, fontsize=7)
        for i in range(proportions.shape[0]):
            for j in range(proportions.shape[1]):
                axis.text(j, i, f"{proportions.iloc[i, j]:.2f}", ha="center", va="center", fontsize=6.5)
        fig.colorbar(image, ax=axis, shrink=0.68, label="Fraction of callable isolates")
        panel_label(axis, letter)

    top_combo = combinations.groupby(["antibiotic", "marker_combination"])["n"].sum().sort_values(ascending=False).reset_index().head(12)
    labels = [f"{SHORT_DRUG[d]} | {c}" for d, c in zip(top_combo.antibiotic, top_combo.marker_combination)]
    axes[0, 2].barh(np.arange(len(top_combo))[::-1], top_combo["n"], color=["#4C78A8" if d == "clarithromycin" else "#D95F02" for d in top_combo.antibiotic])
    axes[0, 2].set_yticks(np.arange(len(top_combo))[::-1], labels, fontsize=6.5)
    axes[0, 2].set_xlabel("Isolates")
    style_axis(axes[0, 2]); panel_label(axes[0, 2], "C")

    prev = prevalence.copy(); y = np.arange(len(prev))
    labels = [f"{SHORT_COHORT[d]} {SHORT_DRUG[a]}" for d, a in zip(prev.dataset_id, prev.antibiotic)]
    for i, row in enumerate(prev.itertuples()):
        axes[1, 0].plot([row.phenotypic_resistance_prevalence, row.marker_resistance_prevalence], [i, i], color="#999999")
        axes[1, 0].scatter(row.phenotypic_resistance_prevalence, i, color="#D64B4B", s=28)
        axes[1, 0].scatter(row.marker_resistance_prevalence, i, color="#2878B5", s=28)
    axes[1, 0].set_xlim(0, 0.75); axes[1, 0].set_yticks(y, labels, fontsize=7)
    axes[1, 0].set_xlabel("Prevalence")
    axes[1, 0].legend(handles=[Line2D([0], [0], marker="o", linestyle="", color="#D64B4B", label="Phenotype R"), Line2D([0], [0], marker="o", linestyle="", color="#2878B5", label="Marker R")], frameon=False, fontsize=7)
    style_axis(axes[1, 0]); panel_label(axes[1, 0], "D")

    lvx_lineage = lineage[lineage["antibiotic"] == "levofloxacin"].copy()
    lvx_lineage["gap"] = lvx_lineage["marker_resistance_prevalence"] - lvx_lineage["phenotype_resistance_prevalence"]
    pivot_gap = lvx_lineage.pivot_table(index="lineage_recomputed", columns="dataset_id", values="gap", fill_value=np.nan).reindex(columns=COHORT_ORDER)
    im = axes[1, 1].imshow(pivot_gap.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-0.6, vmax=0.6)
    axes[1, 1].set_xticks(np.arange(3), [SHORT_COHORT[x] for x in COHORT_ORDER])
    axes[1, 1].set_yticks(np.arange(len(pivot_gap)), [x.replace("SNP_CLUSTER_", "C") for x in pivot_gap.index], fontsize=7)
    fig.colorbar(im, ax=axes[1, 1], shrink=0.72, label="Prevalence difference")
    panel_label(axes[1, 1], "E")

    unique = sample_markers.drop_duplicates(["dataset_id", "isolate_id", "antibiotic", "outcome"])
    outcome_counts = unique.groupby(["dataset_id", "antibiotic", "outcome"]).size().unstack(fill_value=0)
    outcome_counts = outcome_counts.reindex(pd.MultiIndex.from_product([COHORT_ORDER, DRUG_LABELS]), fill_value=0)
    bottom = np.zeros(len(outcome_counts))
    for name in ["TN", "TP", "FP", "FN"]:
        values = outcome_counts[name].to_numpy() if name in outcome_counts else np.zeros(len(outcome_counts))
        axes[1, 2].bar(np.arange(len(outcome_counts)), values, bottom=bottom, color=OUTCOME_COLORS[name], label=name)
        bottom += values
    axes[1, 2].set_xticks(np.arange(len(outcome_counts)), [f"{SHORT_COHORT[d]}\n{SHORT_DRUG[a]}" for d, a in outcome_counts.index], fontsize=7)
    axes[1, 2].set_ylabel("Callable isolates")
    axes[1, 2].legend(frameon=False, fontsize=7, ncol=2)
    style_axis(axes[1, 2]); panel_label(axes[1, 2], "F")
    save_figure(fig, figures, "Figure7_mutation_architecture_prevalence_shift")


def figure8(root: Path, figures: Path, source: Path) -> None:
    influence = numeric(read_csv(root / "results/extended_analysis/ningxia_lvx_leave_one_out_influence.csv"), ["false_susceptible_rate", "balanced_accuracy"])
    bootstrap = numeric(read_csv(root / "results/extended_analysis/external_bootstrap_distributions.csv"), ["balanced_accuracy", "false_susceptible_rate"])
    predictive = numeric(read_csv(root / "results/extended_analysis/predictive_values_by_assumed_prevalence.csv"), ["assumed_prevalence", "ppv", "npv"])
    write_source(influence, source, "Figure8AB_leave_one_out.csv")
    write_source(bootstrap, source, "Figure8CD_bootstrap_distributions.csv")
    write_source(predictive, source, "Figure8EF_predictive_value_curves.csv")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.3), constrained_layout=True)
    fig._compact_panel_labels = True

    baseline = influence[influence["dropped_isolate"] == "NONE_BASELINE"].iloc[0]
    dropped = influence[influence["dropped_isolate"] != "NONE_BASELINE"].reset_index(drop=True)
    for axis, metric, gate, title, letter in [
        (axes[0, 0], "false_susceptible_rate", 0.10, "Ningxia LVX leave-one-out false-susceptible rate", "A"),
        (axes[0, 1], "balanced_accuracy", 0.90, "Ningxia LVX leave-one-out balanced accuracy", "B"),
    ]:
        values = dropped[metric].to_numpy()
        axis.plot(np.arange(len(values)), values, marker="o", markersize=3, linewidth=0.8, color="#D95F02")
        axis.axhline(baseline[metric], color="#222222", linewidth=1.2, label="Full cohort")
        axis.axhline(gate, color="#777777", linestyle="--", linewidth=1, label="Frozen gate")
        axis.set_xlabel("Dropped isolate (ordered by identifier)"); axis.set_ylabel(metric.replace("_", " "))
        axis.legend(frameon=False, fontsize=7)
        style_axis(axis); panel_label(axis, letter)

    estimable_boot = bootstrap[bootstrap["n"].astype(float) >= 20].copy()
    labels, data_ba, data_fsr, colors = [], [], [], []
    for dataset, drug in [("CHINA_NINGXIA_2022", "levofloxacin"), ("ZENODO_10369064", "clarithromycin"), ("ZENODO_10369064", "levofloxacin")]:
        part = estimable_boot[(estimable_boot["dataset_id"] == dataset) & (estimable_boot["antibiotic"] == drug)]
        labels.append(f"{SHORT_COHORT[dataset]} {SHORT_DRUG[drug]}")
        data_ba.append(part["balanced_accuracy"].dropna().to_numpy()); data_fsr.append(part["false_susceptible_rate"].dropna().to_numpy())
        colors.append(cohort_color(dataset))
    for axis, datasets, gate, title, letter in [
        (axes[0, 2], data_ba, 0.90, "External bootstrap balanced accuracy", "C"),
        (axes[1, 0], data_fsr, 0.10, "External bootstrap false-susceptible rate", "D"),
    ]:
        violin = axis.violinplot(datasets, positions=np.arange(len(datasets)), showmedians=True, showextrema=False)
        for body, color in zip(violin["bodies"], colors):
            body.set_facecolor(color); body.set_alpha(0.68)
        violin["cmedians"].set_color("black")
        axis.axhline(gate, color="#777777", linestyle="--", linewidth=1)
        axis.set_xticks(np.arange(len(labels)), labels, rotation=25, ha="right", fontsize=7)
        axis.set_ylim(0, 1)
        axis.set_ylabel("Balanced accuracy" if letter == "C" else "False-susceptible rate")
        style_axis(axis); panel_label(axis, letter)

    curves = predictive[
        (predictive["dataset_id"].isin(EXTERNAL_COHORTS))
        & ~((predictive["dataset_id"] == "CHINA_NINGXIA_2022") & (predictive["antibiotic"] == "clarithromycin"))
    ]
    for axis, metric, title, letter in [(axes[1, 1], "ppv", "PPV across assumed resistance prevalence", "E"), (axes[1, 2], "npv", "NPV across assumed resistance prevalence", "F")]:
        for (dataset, drug), part in curves.groupby(["dataset_id", "antibiotic"], sort=False):
            axis.plot(part["assumed_prevalence"], part[metric], color=cohort_color(dataset), linestyle="-" if drug == "levofloxacin" else "--", label=f"{SHORT_COHORT[dataset]} {SHORT_DRUG[drug]}")
        axis.set_xlim(0.05, 0.75); axis.set_ylim(0, 1.02)
        axis.set_xlabel("Assumed resistance prevalence"); axis.set_ylabel(metric.upper())
        axis.legend(frameon=False, fontsize=7)
        style_axis(axis); panel_label(axis, letter)
    save_figure(fig, figures, "Figure8_robustness_predictive_values")


def supplementary_figures(root: Path, figures: Path, source: Path) -> None:
    figures.mkdir(parents=True, exist_ok=True)
    associations = numeric(read_csv(root / "results/extended_analysis/qc_callability_associations.csv"), [
        "cliffs_delta_callable_minus_uncallable", "mann_whitney_p",
    ])
    fig, ax = plt.subplots(figsize=(9.2, 6.2), constrained_layout=True)
    show = associations[associations["cliffs_delta_callable_minus_uncallable"].notna()].copy()
    show["label"] = [f"{SHORT_COHORT.get(d, d)} {SHORT_DRUG[a]} | {m}" for d, a, m in zip(show.dataset_id, show.antibiotic, show.metric)]
    colors = [cohort_color(x) for x in show.dataset_id]
    ax.barh(np.arange(len(show)), show["cliffs_delta_callable_minus_uncallable"], color=colors)
    ax.axvline(0, color="black", linewidth=1)
    ax.set_yticks(np.arange(len(show)), show["label"], fontsize=7)
    ax.set_xlabel("Cliff's delta: callable minus uncallable")
    style_axis(ax)
    save_figure(fig, figures, "Supplementary_Figure_S1_QC_callability_effects")

    variants = numeric(read_csv(root / "results/extended_analysis/ningxia_gyrA_missense_associations.csv"), ["position", "fisher_p", "bh_q"])
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), constrained_layout=True)
    axes[0].scatter(variants["position"], -np.log10(variants["fisher_p"].clip(lower=1e-300)), c=np.where(variants["frozen_catalogue_marker"].str.lower() == "true", "#D95F02", "#4C78A8"), s=22)
    axes[0].axhline(-math.log10(0.05), color="#777777", linestyle="--")
    axes[0].set_xlabel("gyrA amino-acid position"); axes[0].set_ylabel("−log10 Fisher P")
    style_axis(axes[0]); panel_label(axes[0], "A")
    top = variants.nsmallest(20, "fisher_p").sort_values("fisher_p", ascending=False)
    axes[1].barh(np.arange(len(top)), -np.log10(top["fisher_p"].clip(lower=1e-300)), color=np.where(top["frozen_catalogue_marker"].str.lower() == "true", "#D95F02", "#4C78A8"))
    axes[1].set_yticks(np.arange(len(top)), top["variant"], fontsize=7); axes[1].set_xlabel("−log10 Fisher P")
    style_axis(axes[1]); panel_label(axes[1], "B")
    save_figure(fig, figures, "Supplementary_Figure_S2_gyrA_residual_scan")

    sensitivity = numeric(read_csv(root / "results/sensitivity/phenotype_sensitivity_metrics.csv"), ["false_susceptible_rate", "balanced_accuracy"])
    primary = sensitivity[(sensitivity["stratum_type"] == "COHORT") & sensitivity["dataset_id"].isin(EXTERNAL_COHORTS)]
    scenarios = ["PRIMARY_ORIGINAL_I_EXCLUDED", "EXCLUDE_BORDERLINE_MIC", "RECOMPUTED_I_EXCLUDED", "I_AS_S", "I_AS_R"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), constrained_layout=True)
    for dataset, drug, label, color in [
        ("CHINA_NINGXIA_2022", "levofloxacin", "Ningxia LVX", "#D95F02"),
        ("ZENODO_10369064", "clarithromycin", "Reads CLR", "#3A9D5D"),
        ("ZENODO_10369064", "levofloxacin", "Reads LVX", "#2878B5"),
    ]:
        part = primary[(primary.dataset_id == dataset) & (primary.antibiotic == drug) & primary.scenario.isin(scenarios)].set_index("scenario").reindex(scenarios)
        axes[0].plot(np.arange(len(scenarios)), part["false_susceptible_rate"], marker="o", color=color, label=label)
        axes[1].plot(np.arange(len(scenarios)), part["balanced_accuracy"], marker="o", color=color, label=label)
    for axis, metric, letter in [(axes[0], "False-susceptible rate", "A"), (axes[1], "Balanced accuracy", "B")]:
        axis.set_xticks(np.arange(len(scenarios)), ["Primary", "No borderline", "Recomputed", "I→S", "I→R"], rotation=30, ha="right")
        axis.set_ylim(0, 1); axis.set_ylabel(metric); axis.legend(frameon=False, fontsize=7)
        style_axis(axis); panel_label(axis, letter)
    save_figure(fig, figures, "Supplementary_Figure_S3_phenotype_sensitivity")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--figures", type=Path, default=Path("figures/main"))
    parser.add_argument("--supplementary-figures", type=Path, default=Path("figures/supplementary"))
    parser.add_argument("--source-data", type=Path, default=Path("results/source_data"))
    args = parser.parse_args()
    root = args.root.resolve()
    figures = args.figures if args.figures.is_absolute() else root / args.figures
    supplementary = args.supplementary_figures if args.supplementary_figures.is_absolute() else root / args.supplementary_figures
    source = args.source_data if args.source_data.is_absolute() else root / args.source_data
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
        "axes.labelsize": 9, "figure.dpi": 120, "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    figure1(root, figures, source)
    figure2(root, figures, source)
    figure3(root, figures, source)
    figure4(root, figures, source)
    figure5(root, figures, source)
    figure6(root, figures, source)
    figure7(root, figures, source)
    figure8(root, figures, source)
    supplementary_figures(root, supplementary, source)
    manifest = []
    figure_suffixes = {".pdf", ".png", ".tiff"}
    for directory in [figures, supplementary]:
        for path in sorted(directory.iterdir()):
            if path.is_file() and path.suffix.lower() in figure_suffixes:
                manifest.append({"file": path.relative_to(root).as_posix(), "bytes": path.stat().st_size})
    source.mkdir(parents=True, exist_ok=True)
    with (source / "main_figure_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "bytes"])
        writer.writeheader(); writer.writerows(manifest)


if __name__ == "__main__":
    main()
