#!/usr/bin/env python3
"""Render Supplementary Figure S7 for the 23S callability sensitivity audit."""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = {"callable": "#0072B2", "partial": "#D55E00"}


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.12, 1.08, f"({label.lower()})", transform=axis.transAxes, fontsize=14, fontweight="bold", va="top")


def heatmap(axis: plt.Axes, grid: pd.DataFrame, task: str, label: str) -> None:
    subset = grid[grid["blast_task"] == task]
    identities = sorted(subset["minimum_identity"].unique())
    coverages = sorted(subset["minimum_query_coverage"].unique())
    matrix = np.array([
        [
            int(subset[(subset["minimum_identity"] == identity) & (subset["minimum_query_coverage"] == coverage)]["n_callable_final_qc"].iloc[0])
            for coverage in coverages
        ]
        for identity in identities
    ])
    image = axis.imshow(matrix, cmap="Blues", vmin=0, vmax=max(6, int(matrix.max())), aspect="auto")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            axis.text(column, row, str(matrix[row, column]), ha="center", va="center", fontsize=8,
                      color="white" if matrix[row, column] >= max(6, int(matrix.max())) * 0.6 else "black")
    axis.set_xticks(range(len(coverages)), [f"{value:.2f}" for value in coverages], rotation=45, ha="right")
    axis.set_yticks(range(len(identities)), [f"{value:.2f}" for value in identities])
    axis.set_xlabel("Minimum aligned-query coverage")
    axis.set_ylabel(f"Minimum identity ({task}; final-QC n=57)")
    panel_label(axis, label)
    return image


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    analysis = root / "results/extended_analysis"
    grid = pd.read_csv(analysis / "23s_callability_sensitivity_grid.csv")
    samples = pd.read_csv(analysis / "23s_callability_sensitivity_samples.csv")
    megablast = samples[samples["blast_task"] == "megablast"].copy()

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    figure, axes = plt.subplots(2, 2, figsize=(10.5, 8.2), constrained_layout=True)
    image = heatmap(axes[0, 0], grid, "megablast", "A")
    heatmap(axes[0, 1], grid, "blastn", "B")
    colorbar = figure.colorbar(image, ax=axes[0, :], shrink=0.75, pad=0.02)
    colorbar.set_label("Callable assemblies")

    megablast["coverage_bin"] = megablast["max_coverage"].round(3)
    megablast["identity_bin"] = megablast["max_identity"].round(3)
    grouped = megablast.groupby(["coverage_bin", "identity_bin", "recovery_class"], dropna=False).size().reset_index(name="n")
    color_map = {
        "CALLABLE_MARKER_SPANNING": COLORS["callable"],
        "PARTIAL_23S_NO_MARKER_SPAN": COLORS["partial"],
        "RESCUED_ONLY_UNDER_RELAXED_THRESHOLD": "#009E73",
        "NO_23S_LIKE_HIT": "#999999",
    }
    axes[1, 0].scatter(grouped["coverage_bin"], grouped["identity_bin"],
                       c=grouped["recovery_class"].map(color_map), s=35 + grouped["n"] * 12, alpha=0.85,
                       edgecolor="white", linewidth=0.7)
    for row in grouped.itertuples():
        if row.n > 1:
            axes[1, 0].text(row.coverage_bin, row.identity_bin, str(row.n), ha="center", va="center",
                            fontsize=7, color="white", fontweight="bold")
    axes[1, 0].axhline(0.90, color="#333333", linestyle="--", linewidth=1, label="Frozen identity")
    axes[1, 0].axvline(0.05, color="#666666", linestyle=":", linewidth=1, label="Frozen coverage")
    axes[1, 0].set_xlabel("Maximum aligned-query coverage")
    axes[1, 0].set_ylabel("Maximum hit identity")
    axes[1, 0].set_xlim(-0.02, 1.05)
    axes[1, 0].set_ylim(max(0.65, float(megablast["max_identity"].min()) - 0.03), 1.01)
    axes[1, 0].legend(frameon=False, fontsize=8, loc="lower right")
    panel_label(axes[1, 0], "C")

    categories = ["CALLABLE_MARKER_SPANNING", "PARTIAL_23S_NO_MARKER_SPAN",
                  "RESCUED_ONLY_UNDER_RELAXED_THRESHOLD", "NO_23S_LIKE_HIT"]
    labels = ["Marker-spanning", "Partial; no marker span", "Relaxed-only rescue", "No 23S-like hit"]
    counts = [int((megablast["recovery_class"] == category).sum()) for category in categories]
    bars = axes[1, 1].barh(labels, counts, color=[COLORS["callable"], COLORS["partial"], "#009E73", "#999999"])
    for bar, count in zip(bars, counts):
        axes[1, 1].text(count + 0.7, bar.get_y() + bar.get_height() / 2, str(count), va="center")
    axes[1, 1].set_xlim(0, max(counts) * 1.15)
    axes[1, 1].set_xlabel("Ningxia assemblies")
    axes[1, 1].invert_yaxis()
    panel_label(axes[1, 1], "D")

    output_dir = root / "figures/supplementary"
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / "Supplementary_Figure_S7_23S_callability_sensitivity"
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(stem.with_suffix(".tiff"), dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(figure)

    source_dir = root / "results/source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "23s_callability_sensitivity_grid.csv",
        "23s_callability_sensitivity_samples.csv",
        "23s_callability_sensitivity_summary.csv",
    ):
        shutil.copy2(analysis / name, source_dir / f"Supplementary_Figure_S7_{name.removeprefix('23s_callability_sensitivity_')}")


if __name__ == "__main__":
    main()
