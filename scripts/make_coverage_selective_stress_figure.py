#!/usr/bin/env python3
"""Render Supplementary Figure S8 for coverage-aware transport stress tests."""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COHORT = {
    "HPGP_GLOBAL": "HpGP",
    "CHINA_NINGXIA_2022": "Ningxia",
    "ZENODO_10369064": "Read cohort",
}
DRUG = {"clarithromycin": "CLR", "levofloxacin": "LVX"}
COLORS = {
    "HPGP_GLOBAL": "#0072B2",
    "CHINA_NINGXIA_2022": "#D55E00",
    "ZENODO_10369064": "#009E73",
}


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.12,
        1.08,
        f"({label.lower()})",
        transform=axis.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top",
    )


def ordered_labels(frame: pd.DataFrame) -> list[str]:
    order = []
    for dataset_id in ("HPGP_GLOBAL", "CHINA_NINGXIA_2022", "ZENODO_10369064"):
        for antibiotic in ("clarithromycin", "levofloxacin"):
            if ((frame["dataset_id"] == dataset_id) & (frame["antibiotic"] == antibiotic)).any():
                order.append(f"{COHORT[dataset_id]} {DRUG[antibiotic]}")
    return order


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    analysis = root / "results/extended_analysis"
    yields = pd.read_csv(analysis / "coverage_aware_yield.csv")
    bounds = pd.read_csv(analysis / "callability_identification_bounds.csv")
    selective = pd.read_csv(analysis / "manifold_abstention_metrics.csv")
    differences = pd.read_csv(analysis / "external_performance_differences.csv")

    yields["label"] = yields.apply(
        lambda row: f"{COHORT[row['dataset_id']]} {DRUG[row['antibiotic']]}", axis=1
    )
    labels = ordered_labels(yields)
    yields["label"] = pd.Categorical(yields["label"], categories=labels, ordered=True)
    yields = yields.sort_values("label")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    figure, axes = plt.subplots(3, 2, figsize=(10.8, 12.0), constrained_layout=True)

    # A: end-to-end yield partition.
    y = np.arange(len(yields))
    linked = yields["phenotype_linked_binary_n"].to_numpy(float)
    correct = 100 * yields["correct_n"].to_numpy(float) / linked
    incorrect = 100 * yields["incorrect_n"].to_numpy(float) / linked
    unresolved = 100 * yields["unresolved_n"].to_numpy(float) / linked
    axes[0, 0].barh(y, correct, color="#009E73", label="Correct primary result")
    axes[0, 0].barh(y, incorrect, left=correct, color="#E69F00", label="Incorrect primary result")
    axes[0, 0].barh(y, unresolved, left=correct + incorrect, color="#BDBDBD", label="Unresolved")
    axes[0, 0].set_yticks(y, yields["label"].astype(str))
    axes[0, 0].invert_yaxis()
    axes[0, 0].set_xlim(0, 100)
    axes[0, 0].set_xlabel("Phenotype-linked isolates (%)")
    axes[0, 0].legend(
        frameon=False,
        fontsize=7.5,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.17),
        ncol=3,
    )
    panel_label(axes[0, 0], "A")

    # B/C: logical bounds under unresolved samples.
    for axis, metric, metric_label in (
        (axes[0, 1], "sensitivity", "Sensitivity"),
        (axes[1, 0], "specificity", "Specificity"),
    ):
        subset = bounds[bounds["metric"].eq(metric)].copy()
        subset["label"] = subset.apply(
            lambda row: f"{COHORT[row['dataset_id']]} {DRUG[row['antibiotic']]}", axis=1
        )
        subset["label"] = pd.Categorical(subset["label"], categories=labels, ordered=True)
        subset = subset.sort_values("label")
        yy = np.arange(len(subset))
        for position, row in enumerate(subset.itertuples()):
            color = COLORS[row.dataset_id]
            axis.hlines(
                position,
                100 * row.logical_lower_bound,
                100 * row.logical_upper_bound,
                color=color,
                linewidth=4,
                alpha=0.65,
            )
            if np.isfinite(row.evaluable_estimate):
                axis.plot(100 * row.evaluable_estimate, position, "o", color=color, markersize=6)
        axis.set_yticks(yy, subset["label"].astype(str))
        axis.invert_yaxis()
        axis.set_xlim(0, 102)
        axis.set_xlabel(f"{metric_label} estimate or logical bound (%)")
        axis.axvline(90, color="#555555", linestyle="--", linewidth=1)
    panel_label(axes[0, 1], "B")
    panel_label(axes[1, 0], "C")

    # D: formal external-cohort differences.
    metric_labels = {
        "sensitivity": "Sensitivity",
        "specificity": "Specificity",
        "false_susceptible_rate": "FSR",
        "balanced_accuracy": "Balanced accuracy",
    }
    differences = differences.copy()
    differences["metric_label"] = differences["metric"].map(metric_labels)
    dy = np.arange(len(differences))
    estimate = 100 * differences["absolute_difference"].to_numpy(float)
    low = 100 * differences["difference_ci_low"].to_numpy(float)
    high = 100 * differences["difference_ci_high"].to_numpy(float)
    axes[1, 1].errorbar(
        estimate,
        dy,
        xerr=np.vstack([estimate - low, high - estimate]),
        fmt="o",
        color="#7B3294",
        ecolor="#7B3294",
        capsize=3,
    )
    axes[1, 1].axvline(0, color="#333333", linewidth=1)
    axes[1, 1].set_yticks(dy, differences["metric_label"])
    axes[1, 1].invert_yaxis()
    axes[1, 1].set_xlabel("Ningxia minus read cohort (percentage points)")
    panel_label(axes[1, 1], "D")

    # E/F: phenotype-blind abstention curves for external levofloxacin.
    external = selective[
        selective["dataset_id"].isin(["CHINA_NINGXIA_2022", "ZENODO_10369064"])
        & selective["antibiotic"].eq("levofloxacin")
    ].copy()
    for dataset_id, group in external.groupby("dataset_id", sort=False):
        group = group.sort_values("development_percentile_cutoff")
        axes[2, 0].plot(
            100 * group["development_percentile_cutoff"],
            100 * group["retained_coverage"],
            marker="o",
            color=COLORS[dataset_id],
            label=COHORT[dataset_id],
        )
        axes[2, 1].plot(
            100 * group["retained_coverage"],
            100 * group["false_susceptible_rate"],
            marker="o",
            color=COLORS[dataset_id],
            label=COHORT[dataset_id],
        )
    axes[2, 0].set_xlabel("HpGP development-distance percentile cutoff")
    axes[2, 0].set_ylabel("Retained external coverage (%)")
    axes[2, 0].set_xlim(0, 102)
    axes[2, 0].set_ylim(0, 102)
    axes[2, 0].legend(frameon=False)
    panel_label(axes[2, 0], "E")

    axes[2, 1].axhline(10, color="#333333", linestyle="--", linewidth=1, label="Frozen FSR gate")
    axes[2, 1].set_xlabel("Retained external coverage (%)")
    axes[2, 1].set_ylabel("False-susceptible rate (%)")
    axes[2, 1].set_xlim(35, 102)
    axes[2, 1].set_ylim(0, 72)
    axes[2, 1].legend(frameon=False, fontsize=8)
    panel_label(axes[2, 1], "F")

    output_dir = root / "figures/supplementary"
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / "Supplementary_Figure_S8_coverage_selective_stress"
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(
        stem.with_suffix(".tiff"),
        dpi=300,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(figure)

    source_dir = root / "results/source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "coverage_aware_yield.csv",
        "callability_identification_bounds.csv",
        "manifold_abstention_metrics.csv",
        "external_performance_differences.csv",
    ):
        shutil.copy2(analysis / name, source_dir / f"Supplementary_Figure_S8_{name}")


if __name__ == "__main__":
    main()
