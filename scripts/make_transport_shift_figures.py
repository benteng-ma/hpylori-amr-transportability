#!/usr/bin/env python3
"""Render the four-domain transport-shift main figure and supplements S4-S6."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
from sklearn.metrics import auc, roc_curve

from make_manuscript_figures import (
    COHORT_COLORS,
    panel_label,
    save_figure,
    style_axis,
    write_source,
)


HPGP = "HPGP_GLOBAL"
NINGXIA = "CHINA_NINGXIA_2022"
READ_COHORT = "ZENODO_10369064"
COHORTS = [HPGP, NINGXIA, READ_COHORT]
LABELS = {HPGP: "HpGP", NINGXIA: "Ningxia", READ_COHORT: "Read cohort"}
OUTCOME_COLORS = {"TP": "#2A9D8F", "TN": "#4C78A8", "FP": "#E9C46A", "FN": "#D64B4B"}


def numeric(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    result = frame.copy()
    for column in columns:
        if column in result:
            result[column] = pd.to_numeric(result[column], errors="coerce")
    return result


def read(path: Path, numeric_columns: list[str] | None = None) -> pd.DataFrame:
    frame = pd.read_csv(path)
    return numeric(frame, numeric_columns or [])


def errorbar_point(axis, x, row, color, marker="o", label=None):
    low = row["estimate"] - row["ci_low"] if pd.notna(row["ci_low"]) else 0
    high = row["ci_high"] - row["estimate"] if pd.notna(row["ci_high"]) else 0
    axis.errorbar(
        x,
        row["estimate"],
        yerr=np.array([[low], [high]]),
        fmt=marker,
        color=color,
        markerfacecolor=color,
        markersize=6,
        capsize=3,
        linewidth=1.2,
        label=label,
    )


def figure9(root: Path, figures: Path, source: Path) -> None:
    base = root / "results/transport_shift"
    layer = read(base / "three_layer_transport_summary.csv", ["estimate"])
    discriminator = read(
        base / "cohort_discriminator_summary.csv",
        ["pooled_oof_auc", "group_bootstrap_ci_low", "group_bootstrap_ci_high", "permutation_p"],
    )
    performance = read(
        base / "levofloxacin_lineage_performance.csv",
        ["successes", "total", "estimate", "ci_low", "ci_high"],
    )
    standardized = read(
        base / "lineage_standardized_performance.csv",
        ["estimate", "ci_low", "ci_high", "numerator", "denominator"],
    )
    severity = read(
        base / "ningxia_resistant_mic_severity.csv",
        ["threshold_mg_L", "n_at_or_above", "denominator", "proportion_at_or_above"],
    )
    comparisons = read(base / "dominant_lineage_comparisons.csv", ["fisher_p"])

    write_source(layer, source, "Figure9A_three_layer_summary.csv")
    write_source(discriminator, source, "Figure9B_cohort_discriminator.csv")
    dominant = performance[
        (performance["lineage_recomputed"] == "SNP_CLUSTER_03")
        & performance["metric"].isin(["sensitivity", "specificity", "marker_negative_resistance"])
    ].copy()
    write_source(dominant, source, "Figure9CD_dominant_lineage_performance.csv")
    write_source(comparisons, source, "Figure9CD_exact_comparisons.csv")
    write_source(standardized, source, "Figure9E_lineage_standardized_performance.csv")
    write_source(severity, source, "Figure9F_resistant_MIC_severity.csv")

    fig = plt.figure(figsize=(13.5, 10.3), constrained_layout=True)
    fig._compact_panel_labels = True
    grid = fig.add_gridspec(3, 6, height_ratios=[0.82, 1.15, 1.15])
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0:2])
    ax_c = fig.add_subplot(grid[1, 2:4])
    ax_d = fig.add_subplot(grid[1, 4:6])
    ax_e = fig.add_subplot(grid[2, 0:3])
    ax_f = fig.add_subplot(grid[2, 3:6])

    # A: conceptual/data-supported transport chain.
    ax_a.set_xlim(0, 12)
    ax_a.set_ylim(0, 3.2)
    ax_a.axis("off")
    boxes = [
        (0.15, "1  Analytic availability", "Target callable?\nCLR Ningxia 6/60", "#FCE8D5"),
        (3.15, "2  Population shift", "Study identity from PCs\nAUC 0.90-0.92", "#E6EFF8"),
        (6.15, "3  Conditional transport", "Same SNP cluster; observed\nmapping differs by study", "#E6F4EA"),
        (9.15, "4  Clinical consequence", "16 resistant isolates missed\nmedian lower-bound MIC 16 mg/L", "#F8E1E1"),
    ]
    for index, (x, title, detail, color) in enumerate(boxes):
        box = FancyBboxPatch(
            (x, 0.55), 2.55, 1.85, boxstyle="round,pad=0.04,rounding_size=0.08",
            facecolor=color, edgecolor="#666666", linewidth=1.0,
        )
        ax_a.add_patch(box)
        ax_a.text(x + 1.275, 1.88, title, ha="center", va="center", fontsize=10, fontweight="bold")
        ax_a.text(x + 1.275, 1.18, detail, ha="center", va="center", fontsize=9)
        if index < len(boxes) - 1:
            ax_a.add_patch(FancyArrowPatch((x + 2.60, 1.48), (x + 2.97, 1.48), arrowstyle="-|>", mutation_scale=13, color="#555555"))
    panel_label(ax_a, "A")

    # B: adversarial validation.
    for index, row in discriminator.reset_index(drop=True).iterrows():
        y = 1 - index
        x = row["pooled_oof_auc"]
        ax_b.errorbar(
            x,
            y,
            xerr=np.array([[x - row["group_bootstrap_ci_low"]], [row["group_bootstrap_ci_high"] - x]]),
            fmt="o",
            color=COHORT_COLORS[row["external_dataset"]],
            markersize=7,
            capsize=3,
        )
        ax_b.text(x + 0.012, y, f"{x:.3f}", va="center", fontsize=8)
    ax_b.axvline(0.5, linestyle="--", color="#777777", linewidth=1)
    ax_b.set_xlim(0.45, 1.02)
    ax_b.set_yticks([1, 0], ["HpGP vs Ningxia", "HpGP vs read cohort"])
    ax_b.set_xlabel("Out-of-fold cohort-discriminator AUC")
    style_axis(ax_b)
    panel_label(ax_b, "B")

    # C: dominant-lineage sensitivity and specificity.
    cdata = dominant[dominant["metric"].isin(["sensitivity", "specificity"])]
    metric_x = {"sensitivity": 0, "specificity": 1}
    offsets = {HPGP: -0.20, NINGXIA: 0.0, READ_COHORT: 0.20}
    for dataset in COHORTS:
        for _, row in cdata[cdata["dataset_id"] == dataset].iterrows():
            x = metric_x[row["metric"]] + offsets[dataset]
            errorbar_point(ax_c, x, row, COHORT_COLORS[dataset])
            ax_c.text(x, max(0.02, row["estimate"] - 0.13), f"{int(row['successes'])}/{int(row['total'])}", ha="center", fontsize=7)
    ax_c.axhline(0.90, linestyle=":", color="#777777", linewidth=1)
    ax_c.set_ylim(-0.02, 1.08)
    ax_c.set_xlim(-0.5, 1.5)
    ax_c.set_xticks([0, 1], ["Sensitivity", "Specificity"])
    ax_c.set_ylabel("Estimate (exact 95% CI)")
    ax_c.legend(
        handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=COHORT_COLORS[x], label=LABELS[x]) for x in COHORTS],
        frameon=False, fontsize=7, loc="lower left",
    )
    style_axis(ax_c)
    panel_label(ax_c, "C")

    # D: resistance among marker-negative isolates in the dominant lineage.
    ddata = dominant[dominant["metric"] == "marker_negative_resistance"].set_index("dataset_id")
    for index, dataset in enumerate(COHORTS):
        row = ddata.loc[dataset]
        errorbar_point(ax_d, index, row, COHORT_COLORS[dataset])
        ax_d.text(index, row["estimate"] + 0.10, f"{int(row['successes'])}/{int(row['total'])}", ha="center", fontsize=8)
    p_hpgp = comparisons[(comparisons["metric"] == "marker_negative_resistance") & (comparisons["dataset_b"] == HPGP)]["fisher_p"].iloc[0]
    ax_d.text(0.98, 0.98, f"Ningxia vs HpGP\nFisher P={p_hpgp:.2g}", transform=ax_d.transAxes, ha="right", va="top", fontsize=8)
    ax_d.set_xticks(range(3), [LABELS[x] for x in COHORTS], rotation=15)
    ax_d.set_ylim(-0.02, 0.62)
    ax_d.set_ylabel("P(resistant | frozen marker-negative)")
    style_axis(ax_d)
    panel_label(ax_d, "D")

    # E: case-mix standardization.
    edata = standardized[standardized["target_dataset"] == NINGXIA].copy()
    types = ["target_observed", "HpGP_crude", "HpGP_standardized_to_target_lineages"]
    type_labels = ["Ningxia\nobserved", "HpGP\ncrude", "HpGP standardized\nto Ningxia lineages"]
    colors = [COHORT_COLORS[NINGXIA], COHORT_COLORS[HPGP], "#6AAED6"]
    width = 0.23
    for metric_index, metric in enumerate(["sensitivity", "specificity"]):
        for type_index, estimate_type in enumerate(types):
            row = edata[(edata["metric"] == metric) & (edata["estimate_type"] == estimate_type)].iloc[0]
            x = metric_index + (type_index - 1) * width
            ax_e.bar(x, row["estimate"], width=width * 0.9, color=colors[type_index], edgecolor="white")
            if pd.notna(row["ci_low"]):
                ax_e.errorbar(x, row["estimate"], yerr=[[row["estimate"] - row["ci_low"]], [row["ci_high"] - row["estimate"]]], fmt="none", color="#333333", capsize=2)
    ax_e.set_xticks([0, 1], ["Sensitivity", "Specificity"])
    ax_e.set_ylim(0, 1.08)
    ax_e.set_ylabel("Estimate")
    ax_e.legend(handles=[Patch(facecolor=colors[i], label=type_labels[i]) for i in range(3)], frameon=False, fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    style_axis(ax_e)
    panel_label(ax_e, "E")

    # F: severity spectrum among resistant isolates.
    for outcome, marker in [("TP", "o"), ("FN", "s")]:
        group = severity[severity["outcome"] == outcome].sort_values("threshold_mg_L")
        ax_f.plot(
            group["threshold_mg_L"],
            group["proportion_at_or_above"],
            marker=marker,
            linewidth=2,
            color=OUTCOME_COLORS[outcome],
            label="Detected resistant (TP)" if outcome == "TP" else "Missed resistant (FN)",
        )
        for _, row in group.iterrows():
            offset = 0.055 if outcome == "TP" else -0.075
            ax_f.text(row["threshold_mg_L"], row["proportion_at_or_above"] + offset, f"{int(row['n_at_or_above'])}/{int(row['denominator'])}", ha="center", fontsize=7)
    ax_f.set_xscale("log", base=2)
    ax_f.set_xticks([2, 4, 8, 16, 32], ["2", "4", "8", "16", "32"])
    ax_f.set_ylim(0, 1.12)
    ax_f.set_xlabel("Lower-bound levofloxacin MIC threshold (mg/L)")
    ax_f.set_ylabel("Proportion at or above threshold")
    ax_f.legend(frameon=False, fontsize=8, loc="lower left")
    style_axis(ax_f)
    panel_label(ax_f, "F")

    save_figure(fig, figures, "Figure9_three_layer_transport_shift")


def supplementary_s4(root: Path, output: Path, source: Path) -> None:
    data = read(
        root / "results/transport_shift/levofloxacin_lineage_performance.csv",
        ["successes", "total", "estimate", "ci_low", "ci_high"],
    )
    data = data[data["lineage_recomputed"] != "ALL"].copy()
    write_source(data, source, "Supplementary_Figure_S4_lineage_performance.csv")
    lineages = [f"SNP_CLUSTER_{index:02d}" for index in range(1, 9)]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), constrained_layout=True)
    metrics = [
        ("sensitivity", "Sensitivity", "A"),
        ("specificity", "Specificity", "B"),
        ("false_susceptible_rate", "False-susceptible rate", "C"),
        ("marker_negative_resistance", "P(resistant | marker-negative)", "D"),
    ]
    offsets = {HPGP: -0.18, NINGXIA: 0.0, READ_COHORT: 0.18}
    for axis, (metric, metric_label, letter) in zip(axes.ravel(), metrics):
        subset = data[(data["metric"] == metric) & (data["total"] > 0)]
        for dataset in COHORTS:
            group = subset[subset["dataset_id"] == dataset]
            for _, row in group.iterrows():
                index = lineages.index(row["lineage_recomputed"])
                x = index + offsets[dataset]
                errorbar_point(axis, x, row, COHORT_COLORS[dataset])
        axis.set_xticks(range(8), [f"C{i}" for i in range(1, 9)])
        axis.set_ylim(-0.03, 1.08)
        axis.set_ylabel(f"{metric_label} (exact 95% CI)")
        style_axis(axis)
        panel_label(axis, letter)
    axes[0, 0].legend(
        handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=COHORT_COLORS[x], label=LABELS[x]) for x in COHORTS],
        frameon=False, fontsize=8, loc="lower left",
    )
    save_figure(fig, output, "Supplementary_Figure_S4_lineage_conditional_transport")


def supplementary_s5(root: Path, output: Path, source: Path) -> None:
    base = root / "results/transport_shift"
    pred = read(base / "cohort_discriminator_oof_predictions.csv", ["cohort_binary", "oof_probability", "fold"])
    null = read(base / "cohort_discriminator_permutations.csv", ["null_mean_fold_auc", "observed_mean_fold_auc"])
    coeff = read(base / "cohort_discriminator_coefficients.csv", ["standardized_log_odds_coefficient"])
    write_source(pred, source, "Supplementary_Figure_S5_adversarial_predictions.csv")
    write_source(null, source, "Supplementary_Figure_S5_adversarial_null.csv")
    write_source(coeff, source, "Supplementary_Figure_S5_adversarial_coefficients.csv")
    comparisons = list(pred["comparison"].drop_duplicates())
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.5), constrained_layout=True)
    for axis, comparison, letter in zip(axes[0], comparisons, ["A", "B"]):
        group = pred[pred["comparison"] == comparison]
        fpr, tpr, _ = roc_curve(group["cohort_binary"], group["oof_probability"])
        score = auc(fpr, tpr)
        external = group[group["cohort_binary"] == 1]["dataset_id"].iloc[0]
        axis.plot(
            fpr,
            tpr,
            color=COHORT_COLORS[external],
            linewidth=2,
            label=f"HpGP vs {LABELS[external]}; OOF AUC={score:.3f}",
        )
        axis.plot([0, 1], [0, 1], linestyle="--", color="#777777")
        axis.set_xlabel("False-positive rate")
        axis.set_ylabel("True-positive rate")
        axis.legend(frameon=False)
        style_axis(axis)
        panel_label(axis, letter)
    ax = axes[1, 0]
    for comparison in comparisons:
        group = null[null["comparison"] == comparison]
        external = NINGXIA if "NINGXIA" in comparison else READ_COHORT
        ax.hist(group["null_mean_fold_auc"], bins=30, alpha=0.45, color=COHORT_COLORS[external], label=LABELS[external])
        ax.axvline(group["observed_mean_fold_auc"].iloc[0], color=COHORT_COLORS[external], linewidth=2)
    ax.set_xlabel("Mean fold AUC after cohort-label permutation")
    ax.set_ylabel("Permutations")
    ax.legend(frameon=False)
    style_axis(ax)
    panel_label(ax, "C")
    ax = axes[1, 1]
    matrix = coeff.pivot(index="comparison", columns="feature", values="standardized_log_odds_coefficient")
    matrix = matrix[[f"PC{i}" for i in range(1, 11)]]
    image = ax.imshow(matrix.values, aspect="auto", cmap="coolwarm", vmin=-np.abs(matrix.values).max(), vmax=np.abs(matrix.values).max())
    ax.set_xticks(range(len(matrix.columns)), matrix.columns, rotation=45)
    ax.set_yticks(range(len(matrix.index)), ["Ningxia", "Read cohort"])
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="Log-odds coefficient")
    panel_label(ax, "D")
    save_figure(fig, output, "Supplementary_Figure_S5_adversarial_cohort_shift")


def supplementary_s6(root: Path, output: Path, source: Path) -> None:
    atlas = read(
        root / "results/transport_shift/ningxia_failure_atlas.csv",
        ["atlas_order", "mic_numeric", "development_5nn_percentile"],
    )
    atlas["missense_variant_count"] = atlas["full_gyrA_missense_variants"].fillna("").map(
        lambda value: 0 if not value else len(str(value).split(";"))
    )
    write_source(atlas, source, "Supplementary_Figure_S6_Ningxia_failure_atlas.csv")
    fig = plt.figure(figsize=(13.5, 8.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=[1.1, 1])
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[1, 1])
    ax_d = fig.add_subplot(grid[1, 2])

    tracks = ["outcome", "phenotype", "prediction", "published_genotype_call", "calls_concordant"]
    mappings = {
        "outcome": {"FN": 0, "TP": 1, "FP": 2, "TN": 3},
        "phenotype": {"R": 0, "S": 1},
        "prediction": {"S": 0, "R": 1},
        "published_genotype_call": {"S": 0, "R": 1},
        "calls_concordant": {"no": 0, "yes": 1},
    }
    # Each track is rendered separately to preserve semantic colors.
    palettes = {
        "outcome": [OUTCOME_COLORS[x] for x in ["FN", "TP", "FP", "TN"]],
        "phenotype": ["#B2182B", "#2166AC"],
        "prediction": ["#D9EAF4", "#EF8A62"],
        "published_genotype_call": ["#D9EAF4", "#EF8A62"],
        "calls_concordant": ["#D64B4B", "#2A9D8F"],
    }
    ax_a.set_xlim(0, len(atlas))
    ax_a.set_ylim(0, len(tracks))
    for row_index, track in enumerate(tracks):
        values = atlas[track].fillna("missing").map(mappings[track]).fillna(-1).to_numpy()[None, :]
        cmap = ListedColormap(["#D3D3D3", *palettes[track]])
        ax_a.imshow(values + 1, aspect="auto", interpolation="nearest", extent=[0, len(atlas), row_index, row_index + 1], cmap=cmap, vmin=0, vmax=len(palettes[track]))
    for x in range(len(atlas) + 1):
        ax_a.axvline(x, color="white", linewidth=0.25)
    ax_a.set_yticks(np.arange(len(tracks)) + 0.5, ["Outcome", "Phenotype", "Frozen call", "Source genotype", "Source-recall agreement"])
    ax_a.set_xticks(np.arange(0, len(atlas), 4) + 0.5, atlas["isolate_id"].iloc[::4], rotation=90, fontsize=6)
    ax_a.set_xlabel("Ningxia isolates ordered FN → TP → FP → TN")
    panel_label(ax_a, "A")

    x = np.arange(len(atlas))
    colors = atlas["outcome"].map(OUTCOME_COLORS)
    ax_b.scatter(x, atlas["mic_numeric"], c=colors, s=22)
    ax_b.set_yscale("log", base=2)
    ax_b.set_yticks([0.008, 0.032, 0.125, 0.5, 2, 8, 32], ["0.008", "0.032", "0.125", "0.5", "2", "8", ">32"])
    ax_b.set_xlabel("Atlas order")
    ax_b.set_ylabel("Lower-bound MIC (mg/L)")
    style_axis(ax_b)
    panel_label(ax_b, "B")
    for outcome in ["FN", "TP", "FP", "TN"]:
        group = atlas[atlas["outcome"] == outcome]
        jitter = np.linspace(-0.10, 0.10, max(len(group), 1))
        ax_c.scatter(np.full(len(group), ["FN", "TP", "FP", "TN"].index(outcome)) + jitter, group["development_5nn_percentile"], color=OUTCOME_COLORS[outcome], alpha=0.85)
    ax_c.set_xticks(range(4), ["FN", "TP", "FP", "TN"])
    ax_c.set_ylim(-0.03, 1.03)
    ax_c.set_ylabel("HpGP 5-nearest-neighbour distance percentile")
    style_axis(ax_c)
    panel_label(ax_c, "C")
    groups = [atlas.loc[atlas["outcome"] == outcome, "missense_variant_count"] for outcome in ["FN", "TP", "FP", "TN"]]
    ax_d.boxplot(groups, tick_labels=["FN", "TP", "FP", "TN"], showfliers=False)
    for index, values in enumerate(groups, start=1):
        jitter = np.linspace(-0.10, 0.10, max(len(values), 1))
        ax_d.scatter(np.full(len(values), index) + jitter, values, color=OUTCOME_COLORS[["FN", "TP", "FP", "TN"][index - 1]], alpha=0.65, s=16)
    ax_d.set_ylabel("Full-gyrA missense substitutions per isolate")
    style_axis(ax_d)
    panel_label(ax_d, "D")
    fig.legend(
        handles=[Patch(facecolor=OUTCOME_COLORS[x], label=x) for x in ["FN", "TP", "FP", "TN"]],
        frameon=False, ncol=4, loc="lower center",
    )
    save_figure(fig, output, "Supplementary_Figure_S6_Ningxia_failure_atlas")


def update_manifest(root: Path, figures: Path, supplementary: Path, source: Path) -> None:
    rows = []
    figure_suffixes = {".pdf", ".png", ".tiff"}
    for directory in [figures, supplementary]:
        for path in sorted(directory.glob("*.*")):
            if path.suffix.lower() in figure_suffixes:
                rows.append({"file": path.relative_to(root).as_posix(), "bytes": path.stat().st_size})
    source.mkdir(parents=True, exist_ok=True)
    with (source / "main_figure_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "bytes"])
        writer.writeheader()
        writer.writerows(rows)


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
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "figure.dpi": 120,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    figure9(root, figures, source)
    supplementary_s4(root, supplementary, source)
    supplementary_s5(root, supplementary, source)
    supplementary_s6(root, supplementary, source)
    update_manifest(root, figures, supplementary, source)


if __name__ == "__main__":
    main()
