#!/usr/bin/env python3
"""Render the six prespecified main figures and their source-data tables."""

from __future__ import annotations

import argparse
import csv
import textwrap
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch


COHORT_ORDER = ["HPGP_GLOBAL", "CHINA_NINGXIA_2022", "ZENODO_10369064"]
COHORT_LABELS = {
    "HPGP_GLOBAL": "HpGP\n(catalogue)",
    "CHINA_NINGXIA_2022": "Ningxia\n(primary external)",
    "ZENODO_10369064": "Raw-read cohort\n(method stress)",
}
COHORT_COLORS = {
    "HPGP_GLOBAL": "#2878B5",
    "CHINA_NINGXIA_2022": "#D95F02",
    "ZENODO_10369064": "#3A9D5D",
}
DRUG_LABELS = {"clarithromycin": "Clarithromycin", "levofloxacin": "Levofloxacin"}


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def numeric(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    result = frame.copy()
    for column in columns:
        if column in result:
            result[column] = pd.to_numeric(result[column], errors="coerce")
    return result


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.tiff", dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def write_source(frame: pd.DataFrame, source_dir: Path, name: str) -> None:
    source_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(source_dir / name, index=False)


def style_axis(axis: plt.Axes) -> None:
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(labelsize=8)


def panel_label(axis: plt.Axes, label: str) -> None:
    compact = bool(getattr(axis.figure, "_compact_panel_labels", False))
    x, y, vertical_alignment = (-0.08, 1.03, "bottom") if compact else (-0.12, 1.08, "top")
    journal_label = f"({label.lower()})"
    axis.text(x, y, journal_label, transform=axis.transAxes, fontsize=12, fontweight="bold", va=vertical_alignment)


def figure1(root: Path, figures: Path, source: Path) -> None:
    phenotypes = read_csv(root / "metadata/phenotype_manifest.csv")
    qc = read_csv(root / "results/qc/assembly_qc_with_checkm2.csv")
    callability = numeric(read_csv(root / "results/qc/panel_callability.csv"), ["n_phenotype_linked", "n_callable"])
    status = read_csv(root / "results/qc/zenodo_processing_status.csv")

    assembly_counts = qc.groupby("dataset_id", sort=False)["isolate_id"].nunique().to_dict()
    qc_pass = qc[qc["final_qc_status"] == "PASS"].groupby("dataset_id")["isolate_id"].nunique().to_dict()
    phenotype_counts = phenotypes.groupby("dataset_id")["isolate_id"].nunique().to_dict()
    stages = []
    for dataset in COHORT_ORDER:
        stages.append({
            "dataset_id": dataset,
            "phenotype_linked_isolates": int(phenotype_counts.get(dataset, 0)),
            "sequence_materialized": int(assembly_counts.get(dataset, 0)),
            "final_qc_pass": int(qc_pass.get(dataset, 0)),
            "sequence_source": "paired reads + SKESA" if dataset == "ZENODO_10369064" else "public assembly",
        })
    flow = pd.DataFrame(stages)
    write_source(flow, source, "Figure1A_study_flow.csv")

    phenotype_summary = (
        phenotypes[phenotypes["antibiotic"].isin(DRUG_LABELS)]
        .groupby(["dataset_id", "antibiotic", "susceptibility_original"], dropna=False)
        .size().rename("n").reset_index()
    )
    write_source(phenotype_summary, source, "Figure1B_phenotype_harmonization.csv")
    write_source(callability, source, "Figure1C_callability.csv")
    write_source(status, source, "Figure1D_raw_read_processing.csv")

    fig = plt.figure(figsize=(11.4, 7.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.05, 1])
    ax_flow = fig.add_subplot(grid[0, :])
    ax_pheno = fig.add_subplot(grid[1, 0])
    ax_call = fig.add_subplot(grid[1, 1])

    ax_flow.set_xlim(0, 12)
    ax_flow.set_ylim(0, 4)
    ax_flow.axis("off")
    for index, row in flow.iterrows():
        y = 3.15 - index * 1.22
        color = COHORT_COLORS[row.dataset_id]
        left = FancyBboxPatch((0.2, y - 0.38), 2.6, 0.76, boxstyle="round,pad=0.03", facecolor=color, edgecolor="none", alpha=0.92)
        middle = FancyBboxPatch((4.05, y - 0.38), 2.9, 0.76, boxstyle="round,pad=0.03", facecolor="#F2F3F5", edgecolor="#73777A")
        right = FancyBboxPatch((8.2, y - 0.38), 3.15, 0.76, boxstyle="round,pad=0.03", facecolor="#FFFFFF", edgecolor=color, linewidth=1.8)
        for patch in (left, middle, right):
            ax_flow.add_patch(patch)
        ax_flow.annotate("", xy=(4.0, y), xytext=(2.85, y), arrowprops={"arrowstyle": "->", "color": "#4D4D4D"})
        ax_flow.annotate("", xy=(8.15, y), xytext=(7.0, y), arrowprops={"arrowstyle": "->", "color": "#4D4D4D"})
        ax_flow.text(1.5, y, f"{COHORT_LABELS[row.dataset_id]}\nn={row.phenotype_linked_isolates}", ha="center", va="center", color="white", fontsize=9, fontweight="bold")
        ax_flow.text(5.5, y, f"{row.sequence_source}\nmaterialized n={row.sequence_materialized}", ha="center", va="center", fontsize=8.5)
        ax_flow.text(9.78, y, f"Frozen QC pass n={row.final_qc_pass}\nmarker callability reported", ha="center", va="center", fontsize=8.5)
    ax_flow.set_title("Study design, data lineage and frozen genome-quality gates", fontsize=12, loc="left")
    panel_label(ax_flow, "A")

    pheno_pivot = phenotype_summary.pivot_table(index=["dataset_id", "antibiotic"], columns="susceptibility_original", values="n", fill_value=0)
    pheno_pivot = pheno_pivot.reindex(pd.MultiIndex.from_product([COHORT_ORDER, DRUG_LABELS]))
    x = np.arange(len(pheno_pivot))
    bottom = np.zeros(len(pheno_pivot))
    label_colors = {"S": "#4C78A8", "I": "#F2CF5B", "R": "#D64B4B", "": "#BDBDBD"}
    for label in ["S", "I", "R", ""]:
        values = pheno_pivot[label].to_numpy() if label in pheno_pivot else np.zeros(len(x))
        if np.any(values):
            ax_pheno.bar(x, values, bottom=bottom, color=label_colors[label], label=label or "Missing")
            bottom += values
    short_cohorts = {"HPGP_GLOBAL": "HpGP", "CHINA_NINGXIA_2022": "Ningxia", "ZENODO_10369064": "Reads"}
    short_drugs = {"clarithromycin": "CLR", "levofloxacin": "LVX"}
    ax_pheno.set_xticks(x, [f"{short_cohorts[d]}\n{short_drugs[a]}" for d, a in pheno_pivot.index], rotation=0)
    ax_pheno.set_ylabel("Phenotype-linked isolates")
    ax_pheno.set_title("Published phenotype labels retained")
    ax_pheno.legend(frameon=False, ncol=3, fontsize=8)
    style_axis(ax_pheno)
    panel_label(ax_pheno, "B")

    callability["callability_percent"] = 100 * callability["n_callable"] / callability["n_phenotype_linked"]
    for drug_index, drug in enumerate(DRUG_LABELS):
        values = callability[callability["antibiotic"] == drug].set_index("dataset_id").reindex(COHORT_ORDER)
        offset = (drug_index - 0.5) * 0.34
        ax_call.bar(np.arange(3) + offset, values["callability_percent"], width=0.32, label=DRUG_LABELS[drug])
    ax_call.axhline(90, color="#777777", linestyle="--", linewidth=1)
    ax_call.set_ylim(0, 105)
    ax_call.set_ylabel("Callable predictions (%)")
    ax_call.set_xticks(np.arange(3), [COHORT_LABELS[d].splitlines()[0] for d in COHORT_ORDER])
    ax_call.set_title("Target callability is not recoded as wild type")
    ax_call.legend(frameon=False, fontsize=8)
    style_axis(ax_call)
    panel_label(ax_call, "C")
    save_figure(fig, figures, "Figure1_study_design_data_lineage")


def figure2(root: Path, figures: Path, source: Path) -> None:
    lineages = numeric(read_csv(root / "results/lineage_validation/lineage_assignments.csv"), ["PC1", "PC2"])
    groups = read_csv(root / "results/qc/near_clone_groups.csv")
    pairs = numeric(read_csv(root / "results/qc/pairwise_relatedness_candidates.csv"), [
        "mash_distance", "skani_ani_percent", "skani_align_fraction_a", "skani_align_fraction_b",
    ])
    write_source(lineages, source, "Figure2A_lineage_embedding.csv")
    write_source(groups, source, "Figure2B_near_clone_groups.csv")
    write_source(pairs, source, "Figure2C_relatedness_candidates.csv")

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True, gridspec_kw={"width_ratios": [1.35, 1, 0.9]})
    ax, ax_mix, ax_clone = axes
    for dataset in COHORT_ORDER:
        part = lineages[lineages["dataset_id"] == dataset]
        ax.scatter(part["PC1"], part["PC2"], s=18, alpha=0.72, color=COHORT_COLORS[dataset], label=COHORT_LABELS[dataset].replace("\n", " "), edgecolors="none")
    ax.set_xlabel("Development-fitted component 1")
    ax.set_ylabel("Development-fitted component 2")
    ax.set_title("Phenotype-blind core-SNP embedding")
    ax.legend(frameon=False, fontsize=7, markerscale=1.3)
    style_axis(ax)
    panel_label(ax, "A")

    mix = pd.crosstab(lineages["lineage_recomputed"], lineages["dataset_id"], normalize="columns").reindex(columns=COHORT_ORDER, fill_value=0)
    image = ax_mix.imshow(mix.to_numpy(), aspect="auto", cmap="Blues", vmin=0, vmax=max(0.01, float(mix.to_numpy().max())))
    ax_mix.set_xticks(np.arange(3), [COHORT_LABELS[d].splitlines()[0] for d in COHORT_ORDER], rotation=35, ha="right")
    ax_mix.set_yticks(np.arange(len(mix)), mix.index)
    ax_mix.set_title("Fixed-cluster composition")
    fig.colorbar(image, ax=ax_mix, label="Within-cohort proportion", shrink=0.75)
    panel_label(ax_mix, "B")

    sizes = groups.groupby("near_clone_group").size().sort_values(ascending=False)
    component_counts = Counter(sizes.to_list())
    x_values = sorted(component_counts)
    ax_clone.bar([str(value) for value in x_values], [component_counts[value] for value in x_values], color="#6B6ECF")
    cross = groups.groupby("near_clone_group")["dataset_id"].nunique()
    ax_clone.text(0.98, 0.96, f"Cross-cohort components: {(cross > 1).sum()}", transform=ax_clone.transAxes, ha="right", va="top", fontsize=8)
    ax_clone.set_xlabel("Near-clone component size")
    ax_clone.set_ylabel("Number of components")
    ax_clone.set_title("Frozen Mash/skani groups")
    style_axis(ax_clone)
    panel_label(ax_clone, "C")
    save_figure(fig, figures, "Figure2_population_structure_relatedness")


def forest_panel(axis: plt.Axes, rows: pd.DataFrame, metric: str, title: str) -> None:
    labels = []
    y = np.arange(len(rows))[::-1]
    for position, (_, row) in zip(y, rows.iterrows()):
        cohort_label = COHORT_LABELS.get(row["dataset_id"], row["dataset_id"]).replace("\n", " ")
        labels.append(f"{cohort_label} (n={int(row['n'])})")
        if not bool(row["estimable"]):
            axis.text(
                0.02, position,
                f"Not estimable ({int(row['n_resistant'])} R / {int(row['n_susceptible'])} S)",
                va="center", fontsize=7.5, color="#777777",
            )
            continue
        value = row[metric]
        low = row[f"{metric}_ci_low"]
        high = row[f"{metric}_ci_high"]
        axis.errorbar(value, position, xerr=[[value - low], [high - value]], fmt="o", color=COHORT_COLORS.get(row["dataset_id"], "#333333"), capsize=3)
    axis.axvline(0.90, color="#777777", linestyle="--", linewidth=1)
    axis.set_xlim(0, 1.02)
    axis.set_yticks(y, labels)
    axis.set_xlabel("Estimate (exact 95% CI)")
    axis.set_title(title)
    style_axis(axis)


def figure3(root: Path, figures: Path, source: Path) -> None:
    metrics = numeric(read_csv(root / "results/external_validation/frozen_panel_metrics.csv"), [
        "n", "n_resistant", "n_susceptible", "sensitivity", "sensitivity_ci_low", "sensitivity_ci_high",
        "specificity", "specificity_ci_low", "specificity_ci_high", "false_susceptible_rate",
    ])
    cohorts = metrics[(metrics["stratum_type"] == "COHORT") & metrics["dataset_id"].isin(COHORT_ORDER)].copy()
    cohorts["estimable"] = (cohorts["n_resistant"] >= 10) & (cohorts["n_susceptible"] >= 10)
    cohorts["dataset_rank"] = cohorts["dataset_id"].map({name: index for index, name in enumerate(COHORT_ORDER)})
    cohorts = cohorts.sort_values(["antibiotic", "dataset_rank"])
    write_source(cohorts.drop(columns="dataset_rank"), source, "Figure3_external_panel_performance.csv")
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.5), constrained_layout=True, sharex=True)
    for row_index, drug in enumerate(DRUG_LABELS):
        part = cohorts[cohorts["antibiotic"] == drug]
        forest_panel(axes[row_index, 0], part, "sensitivity", f"{DRUG_LABELS[drug]} sensitivity")
        forest_panel(axes[row_index, 1], part, "specificity", f"{DRUG_LABELS[drug]} specificity")
        panel_label(axes[row_index, 0], chr(ord("A") + row_index * 2))
        panel_label(axes[row_index, 1], chr(ord("B") + row_index * 2))
    save_figure(fig, figures, "Figure3_external_performance_frozen_catalogues")


def figure4(root: Path, figures: Path, source: Path) -> None:
    folds = numeric(read_csv(root / "results/lineage_validation/leakage_benchmark_folds.csv"), ["balanced_accuracy", "auroc_probability"])
    write_source(folds, source, "Figure4_validation_design_folds.csv")
    splits = ["RANDOM_ISOLATE_SPLIT", "CLONE_GROUPED_SPLIT", "LEAVE_COUNTRY_OUT", "LEAVE_LINEAGE_OUT", "HPGP_TO_EXTERNAL"]
    models = ["MUTATION_ONLY_LOGISTIC", "LINEAGE_ONLY_LOGISTIC", "MUTATION_PLUS_LINEAGE_LOGISTIC"]
    model_labels = {"MUTATION_ONLY_LOGISTIC": "Mutation", "LINEAGE_ONLY_LOGISTIC": "Lineage", "MUTATION_PLUS_LINEAGE_LOGISTIC": "Mutation + lineage"}
    colors = {models[0]: "#2878B5", models[1]: "#E39C37", models[2]: "#7A5195"}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True, sharey=True)
    rng = np.random.default_rng(20260830)
    for axis, drug in zip(axes, DRUG_LABELS):
        part = folds[folds["antibiotic"] == drug]
        for model_index, model in enumerate(models):
            for split_index, split in enumerate(splits):
                values = part[(part["model"] == model) & (part["split_type"] == split)]["balanced_accuracy"].dropna().to_numpy()
                if not len(values):
                    continue
                x = split_index + (model_index - 1) * 0.22
                jitter = rng.normal(0, 0.025, size=len(values))
                axis.scatter(x + jitter, values, s=12, alpha=0.35, color=colors[model])
                axis.errorbar(x, np.mean(values), yerr=np.std(values, ddof=1) if len(values) > 1 else 0, fmt="D", color=colors[model], capsize=3, markersize=4)
        axis.set_xticks(np.arange(len(splits)), ["Random", "Clone-\ngrouped", "Leave-\ncountry", "Leave-\nlineage", "External"], rotation=0)
        axis.set_ylim(0, 1.02)
        axis.axhline(0.5, color="#AAAAAA", linestyle=":")
        axis.set_title(DRUG_LABELS[drug])
        axis.set_ylabel("Balanced accuracy")
        style_axis(axis)
    axes[0].legend(handles=[Line2D([0], [0], marker="D", linestyle="", color=colors[m], label=model_labels[m]) for m in models], frameon=False, fontsize=8)
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    save_figure(fig, figures, "Figure4_validation_leakage_benchmark")


def figure5(root: Path, figures: Path, source: Path) -> None:
    samples = numeric(read_csv(root / "results/external_validation/sample_level_predictions.csv"), ["mic_numeric"])
    sensitivity = numeric(read_csv(root / "results/sensitivity/phenotype_sensitivity_metrics.csv"), [
        "sensitivity", "specificity", "false_susceptible_rate",
    ])
    audit = numeric(read_csv(root / "results/external_validation/ningxia_published_genotype_audit_summary.csv"), [
        "n", "n_agree", "call_concordance",
    ])
    mic = samples[(samples["dataset_id"] == "CHINA_NINGXIA_2022") & samples["mic_numeric"].notna() & (samples["mic_numeric"] > 0)].copy()
    write_source(mic, source, "Figure5A_MIC_sample_data.csv")
    write_source(sensitivity, source, "Figure5B_breakpoint_sensitivity.csv")
    write_source(audit, source, "Figure5C_Ningxia_genotype_reproducibility.csv")
    fig, axes = plt.subplots(1, 4, figsize=(15.8, 4.5), constrained_layout=True, gridspec_kw={"width_ratios": [1, 1, 1.3, 0.95]})
    rng = np.random.default_rng(20260830)
    for axis, drug, label in zip(axes[:2], DRUG_LABELS, ["A", "B"]):
        part = mic[mic["antibiotic"] == drug]
        for x, prediction in enumerate(["S", "R"]):
            values = part[part["prediction"] == prediction]
            jitter = rng.normal(0, 0.055, size=len(values))
            colors = np.where(values["phenotype"] == "R", "#D64B4B", "#4C78A8")
            axis.scatter(x + jitter, values["mic_numeric"], color=colors, s=28, alpha=0.78, edgecolor="white", linewidth=0.3)
        axis.set_yscale("log", base=2)
        axis.set_xticks([0, 1], ["Marker S", "Marker R"])
        axis.set_ylabel("MIC (mg/L; inequality retained in source)")
        plotted_n = int((part["callable"].str.lower() == "yes").sum())
        axis.set_title(f"{DRUG_LABELS[drug]} (callable n={plotted_n})")
        style_axis(axis)
        panel_label(axis, label)
    primary = sensitivity[(sensitivity["stratum_type"] == "COHORT") & sensitivity["dataset_id"].isin(COHORT_ORDER)]
    scenarios = ["PRIMARY_ORIGINAL_I_EXCLUDED", "EXCLUDE_BORDERLINE_MIC", "RECOMPUTED_I_EXCLUDED", "I_AS_S", "I_AS_R"]
    series = [
        ("CHINA_NINGXIA_2022", "levofloxacin", "Ningxia LVX", "#D95F02"),
        ("ZENODO_10369064", "clarithromycin", "Reads CLR", "#3A9D5D"),
        ("ZENODO_10369064", "levofloxacin", "Reads LVX", "#2878B5"),
    ]
    for dataset, drug, label, color in series:
        part = primary[
            (primary["dataset_id"] == dataset)
            & (primary["antibiotic"] == drug)
            & primary["scenario"].isin(scenarios)
        ].set_index("scenario").reindex(scenarios)
        axes[2].plot(np.arange(len(scenarios)), part["false_susceptible_rate"], marker="o", label=label, color=color)
    axes[2].axhline(0.10, color="#777777", linestyle="--", linewidth=1)
    axes[2].set_xticks(np.arange(len(scenarios)), ["Primary", "No borderline", "Recomputed", "I→S", "I→R"], rotation=35, ha="right")
    axes[2].set_ylim(-0.02, 0.63)
    axes[2].set_ylabel("False-susceptible rate")
    axes[2].set_title("External-cohort label sensitivity")
    axes[2].legend(frameon=False, fontsize=7)
    style_axis(axes[2])
    panel_label(axes[2], "C")

    audit_rows = audit[audit["row_type"].isin(["CALL_REPRODUCIBILITY", "INDEPENDENT_CALLER_REPRODUCIBILITY"])].copy()
    source_recall = audit_rows[
        (audit_rows["antibiotic"] == "levofloxacin")
        & (audit_rows["row_type"] == "CALL_REPRODUCIBILITY")
    ].iloc[0]
    caller_recall = audit_rows[audit_rows["row_type"] == "INDEPENDENT_CALLER_REPRODUCIBILITY"].iloc[0]
    audit_display = pd.DataFrame([
        {"label": "Published G-R\nvs deposited recall", "agreement": source_recall["call_concordance"], "n": source_recall["n"]},
        {"label": "Target BLAST\nvs Snippy", "agreement": caller_recall["call_concordance"], "n": caller_recall["n"]},
    ])
    bars = axes[3].bar(np.arange(2), audit_display["agreement"], color=["#D95F02", "#3A9D5D"], width=0.64)
    for bar, row in zip(bars, audit_display.itertuples()):
        axes[3].text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.035, f"{row.agreement:.1%}\n(n={int(row.n)})", ha="center", va="bottom", fontsize=8)
    axes[3].set_ylim(0, 1.12)
    axes[3].set_xticks(np.arange(2), audit_display["label"], fontsize=7.5)
    axes[3].set_ylabel("Call agreement")
    axes[3].set_title("Ningxia LVX reproducibility")
    style_axis(axes[3])
    panel_label(axes[3], "D")
    save_figure(fig, figures, "Figure5_MIC_false_susceptible_sensitivity")


def figure6(root: Path, figures: Path, source: Path) -> None:
    classification = read_csv(root / "results/external_validation/transportability_classification.csv")
    metrics = numeric(read_csv(root / "results/external_validation/frozen_panel_metrics.csv"), [
        "n", "n_resistant", "n_susceptible", "sensitivity", "specificity", "false_susceptible_rate",
    ])
    secondary = read_csv(root / "results/external_validation/secondary_drug_boundaries.csv")
    write_source(classification, source, "Figure6A_transportability_classification.csv")
    write_source(metrics, source, "Figure6B_cohort_lineage_boundaries.csv")
    write_source(secondary, source, "Figure6C_secondary_drug_boundaries.csv")
    eligible = metrics[metrics["stratum_type"].isin(["COHORT", "LINEAGE"])].copy()
    eligible["row_label"] = np.where(eligible["stratum_type"] == "COHORT", eligible["stratum"].map(lambda value: COHORT_LABELS.get(value, value).replace("\n", " ")), eligible["stratum"])
    eligible["estimable"] = ~(
        (eligible["stratum_type"] == "COHORT")
        & eligible["dataset_id"].isin(["CHINA_NINGXIA_2022", "ZENODO_10369064"])
        & ((eligible["n_resistant"] < 10) | (eligible["n_susceptible"] < 10))
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, max(4.5, 0.26 * len(eligible) + 1.8)), constrained_layout=True, gridspec_kw={"width_ratios": [1.55, 1]})
    metric_names = ["sensitivity", "specificity", "false_susceptible_rate"]
    rows = []
    row_labels = []
    for drug in DRUG_LABELS:
        for _, row in eligible[eligible["antibiotic"] == drug].iterrows():
            rows.append([row[name] for name in metric_names] if row["estimable"] else [np.nan, np.nan, np.nan])
            suffix = "" if row["estimable"] else f" [NE; n={int(row['n'])}]"
            row_labels.append(f"{DRUG_LABELS[drug][:3]} | {row['row_label']}{suffix}")
    matrix = np.asarray(rows, dtype=float) if rows else np.empty((0, 3))
    # Map all metrics to a common gate score: 0.5 is the frozen boundary,
    # values above it pass in the clinically favourable direction.
    colour_matrix = np.full_like(matrix, np.nan)
    if colour_matrix.size:
        colour_matrix[:, :2] = np.clip(0.5 + (matrix[:, :2] - 0.90) / 0.20, 0, 1)
        colour_matrix[:, 2] = np.clip(0.5 + (0.10 - matrix[:, 2]) / 0.20, 0, 1)
    colour_map = plt.get_cmap("RdYlGn").copy()
    colour_map.set_bad("#D9D9D9")
    image = axes[0].imshow(np.ma.masked_invalid(colour_matrix), aspect="auto", cmap=colour_map, vmin=0, vmax=1)
    axes[0].set_xticks(np.arange(3), ["Sensitivity", "Specificity", "False-susceptible\nrate"])
    axes[0].set_yticks(np.arange(len(row_labels)), row_labels, fontsize=7)
    axes[0].set_title("Prespecified cohort and lineage boundaries")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            label = "NE" if np.isnan(matrix[i, j]) else f"{matrix[i, j]:.2f}"
            axes[0].text(j, i, label, ha="center", va="center", fontsize=6.5)
    colorbar = fig.colorbar(image, ax=axes[0], shrink=0.7, label="Frozen gate score (higher is better)")
    colorbar.set_ticks([0, 0.5, 1], labels=["Fail", "Gate", "Pass"])
    panel_label(axes[0], "A")

    axes[1].axis("off")
    axes[1].set_title("Frozen transportability interpretation", loc="left")
    y = 0.95
    for _, row in classification.iterrows():
        color = "#3A9D5D" if row["transportability_label"] == "ROBUSTLY_TRANSPORTABLE" else "#D95F02"
        axes[1].text(0.02, y, DRUG_LABELS.get(row["antibiotic"], row["antibiotic"]), transform=axes[1].transAxes, fontsize=11, fontweight="bold")
        axes[1].text(0.02, y - 0.07, row["transportability_label"].replace("_", " "), transform=axes[1].transAxes, fontsize=9, color=color, fontweight="bold")
        axes[1].text(0.02, y - 0.14, textwrap.fill(row["reason"], width=52), transform=axes[1].transAxes, fontsize=8, va="top")
        y -= 0.32
    axes[1].text(0.02, y, "Other antibiotics", transform=axes[1].transAxes, fontsize=10, fontweight="bold")
    axes[1].text(0.02, y - 0.07, "Insufficient data for a reproducible frozen cross-cohort panel", transform=axes[1].transAxes, fontsize=8, wrap=True)
    axes[1].text(0.02, 0.03, "Geographic boundary: both external cohorts were Chinese; no global external-validation claim.", transform=axes[1].transAxes, fontsize=8, color="#555555", wrap=True)
    panel_label(axes[1], "B")
    save_figure(fig, figures, "Figure6_transportability_boundaries")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--figures", type=Path, default=Path("figures/main"))
    parser.add_argument("--source-data", type=Path, default=Path("results/source_data"))
    args = parser.parse_args()
    root = args.root.resolve()
    figures = args.figures if args.figures.is_absolute() else root / args.figures
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
    manifest = []
    for path in sorted(figures.glob("Figure*.*")):
        manifest.append({"file": path.relative_to(root).as_posix(), "bytes": path.stat().st_size})
    with (source / "main_figure_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "bytes"])
        writer.writeheader()
        writer.writerows(manifest)


if __name__ == "__main__":
    main()
