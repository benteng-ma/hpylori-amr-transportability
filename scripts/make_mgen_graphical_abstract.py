#!/usr/bin/env python3
"""Create the optional 3:1 graphical abstract for Microbial Genomics."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


def read_row(path: Path, **filters: str) -> pd.Series:
    frame = pd.read_csv(path)
    for column, value in filters.items():
        frame = frame[frame[column].astype(str) == value]
    if len(frame) != 1:
        raise ValueError(f"Expected one row in {path} for {filters}; observed {len(frame)}")
    return frame.iloc[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, default=Path("submission/microbial_genomics/assets"))
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output_dir if args.output_dir.is_absolute() else root / args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    callability = read_row(
        root / "results/qc/panel_callability.csv",
        dataset_id="CHINA_NINGXIA_2022",
        antibiotic="clarithromycin",
    )
    discriminators = pd.read_csv(root / "results/transport_shift/cohort_discriminator_summary.csv")
    performance = read_row(
        root / "results/external_validation/frozen_panel_metrics.csv",
        stratum_type="COHORT",
        dataset_id="CHINA_NINGXIA_2022",
        antibiotic="levofloxacin",
    )
    severity = read_row(
        root / "results/transport_shift/ningxia_resistant_mic_severity.csv",
        outcome="FN",
        threshold_mg_L="8",
    )

    auc_low = discriminators["pooled_oof_auc"].min()
    auc_high = discriminators["pooled_oof_auc"].max()
    sensitivity = 100 * float(performance["sensitivity"])
    false_susceptible = 100 * float(performance["false_susceptible_rate"])

    fig, axis = plt.subplots(figsize=(15, 5), constrained_layout=True)
    axis.set_xlim(0, 15)
    axis.set_ylim(0, 5)
    axis.axis("off")

    axis.text(
        7.5,
        4.62,
        r"When $\it{H.\ pylori}$ resistance genotyping travels, every link must hold",
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold",
        color="#172B4D",
    )

    boxes = [
        (
            0.25,
            "1  ANALYTIC AVAILABILITY",
            "Is the target callable?",
            f"23S: {int(callability['n_callable'])}/{int(callability['n_phenotype_linked'])}\nNingxia assemblies",
            "#FDEBD0",
            "#D95F02",
        ),
        (
            4.02,
            "2  POPULATION SHIFT",
            "Are cohorts exchangeable?",
            f"Study identity remained learnable\nPC-only AUC {auc_low:.3f}-{auc_high:.3f}",
            "#EAF2F8",
            "#2878B5",
        ),
        (
            7.79,
            "3  CONDITIONAL TRANSPORT",
            "Does genotype map to phenotype?",
            f"Ningxia levofloxacin\nsensitivity {sensitivity:.1f}% | false-S {false_susceptible:.1f}%",
            "#E8F5E9",
            "#3A9D5D",
        ),
        (
            11.56,
            "4  CLINICAL CONSEQUENCE",
            "How severe are the misses?",
            f"{int(severity['n_at_or_above'])}/{int(severity['denominator'])} missed resistant isolates\nMIC lower bound at least 8 mg/L",
            "#FBE9E7",
            "#C74343",
        ),
    ]

    for index, (x, heading, question, result, fill, edge) in enumerate(boxes):
        patch = FancyBboxPatch(
            (x, 1.35),
            3.15,
            2.55,
            boxstyle="round,pad=0.04,rounding_size=0.10",
            facecolor=fill,
            edgecolor=edge,
            linewidth=2.0,
        )
        axis.add_patch(patch)
        axis.text(x + 1.575, 3.55, heading, ha="center", va="center", fontsize=10.2, fontweight="bold", color=edge)
        axis.text(x + 1.575, 2.96, question, ha="center", va="center", fontsize=11.2, fontweight="bold", color="#172B4D")
        axis.text(x + 1.575, 2.12, result, ha="center", va="center", fontsize=10.5, color="#172B4D", linespacing=1.3)
        if index < len(boxes) - 1:
            axis.add_patch(
                FancyArrowPatch(
                    (x + 3.20, 2.62),
                    (boxes[index + 1][0] - 0.05, 2.62),
                    arrowstyle="-|>",
                    mutation_scale=17,
                    linewidth=1.6,
                    color="#5B6770",
                )
            )

    axis.text(
        7.5,
        0.63,
        "Audit the chain before clinical use: uncallable is not susceptible, and internal accuracy is not transportability.",
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
        color="#172B4D",
    )

    stem = output / "MGEN_Graphical_Abstract"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


if __name__ == "__main__":
    main()
