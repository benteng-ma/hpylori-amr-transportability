#!/usr/bin/env python3
"""Reproduce QC-conditioned Tables S2/S6 and Supplementary Figure S1.

Run with the supplied 1,052-row qc_callability_samples.csv.
Only the final-genome-QC-pass rows enter assembly-quality associations.
Target callability in the paired-read cohort is the primary raw-read call;
assembly-quality metrics describe the corresponding support assembly.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu
import sklearn
from sklearn.metrics import roc_auc_score

COHORTS = ["CHINA_NINGXIA_2022", "HPGP_GLOBAL", "ZENODO_10369064"]
DRUGS = ["clarithromycin", "levofloxacin"]
METRICS = ["contigs", "log10_n50_bp", "completeness_percent", "contamination_percent"]
LABELS = ["Contig count", "log10 N50 (bp)", "Completeness (%)", "Contamination (%)"]


def read_samples(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == '.xlsx':
        rows = pd.read_excel(path, sheet_name='S06_QC_Associations', engine='openpyxl')
        data = rows.loc[rows['source_file'].eq('results/extended_analysis/qc_callability_samples.csv')].copy()
        data = data.drop(columns=['source_file']).dropna(axis=1, how='all')
    else:
        data = pd.read_csv(path, float_precision="round_trip")
    data[METRICS + ['callable_binary']] = data[METRICS + ['callable_binary']].apply(pd.to_numeric)
    assert len(data) == 1052
    assert not data.duplicated(["dataset_id", "isolate_id", "antibiotic"]).any()
    assert set(data["final_qc_status"]) == {"PASS", "FAIL"}
    assert set(data["callable"]) == {"yes", "no"}
    assert (data["callable_binary"] == data["callable"].eq("yes").astype(int)).all()
    assert not data.loc[data["final_qc_status"].eq("FAIL"), "callable_binary"].any()
    assert not data[METRICS].isna().any().any()
    return data


def summaries(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    accounting = []
    associations = []
    for cohort in COHORTS:
        for drug in DRUGS:
            all_records = data.loc[data["dataset_id"].eq(cohort) & data["antibiotic"].eq(drug)]
            eligible = all_records.loc[all_records["final_qc_status"].eq("PASS")]
            n_linked = len(all_records)
            n_pass = len(eligible)
            n_callable = int(eligible["callable_binary"].sum())
            accounting.append({
                "dataset_id": cohort, "antibiotic": drug, "n_linked": n_linked,
                "n_qc_pass": n_pass, "n_qc_fail": n_linked - n_pass,
                "n_callable_after_qc": n_callable, "n_uncallable_after_qc": n_pass - n_callable,
                "conditional_callability": n_callable / n_pass,
                "end_to_end_callability": n_callable / n_linked,
            })
            for metric in METRICS:
                yes = eligible.loc[eligible["callable_binary"].eq(1), metric].to_numpy(dtype=float)
                no = eligible.loc[eligible["callable_binary"].eq(0), metric].to_numpy(dtype=float)
                row = {
                    "dataset_id": cohort, "antibiotic": drug, "metric": metric,
                    "n_callable": len(yes), "n_uncallable": len(no),
                    "callable_median": float(np.median(yes)) if len(yes) else np.nan,
                    "uncallable_median": float(np.median(no)) if len(no) else np.nan,
                    "mann_whitney_u": np.nan, "mann_whitney_p": np.nan,
                    "cliffs_delta_callable_minus_uncallable": np.nan,
                    "auc_metric_predicting_callable": np.nan,
                    "analysis_scope": "final_genome_QC_PASS",
                }
                # Same minimum group size and two-sided auto method as the prior analysis.
                if len(yes) >= 3 and len(no) >= 3:
                    test = mannwhitneyu(yes, no, alternative="two-sided", method="auto")
                    row["mann_whitney_u"] = float(test.statistic)
                    row["mann_whitney_p"] = float(test.pvalue)
                    row["cliffs_delta_callable_minus_uncallable"] = float(
                        2 * test.statistic / (len(yes) * len(no)) - 1
                    )
                    row["auc_metric_predicting_callable"] = float(
                        roc_auc_score(eligible["callable_binary"], eligible[metric])
                    )
                associations.append(row)
    return pd.DataFrame(accounting), pd.DataFrame(associations)


def draw_figure(associations: pd.DataFrame, output: Path) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.labelsize": 10,
        "axes.titlesize": 11, "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(2, 1, figsize=(7.6, 5.8), sharex=True, constrained_layout=True)
    panels = [
        ("CHINA_NINGXIA_2022", "clarithromycin", "#D95F02",
         "(a) Ningxia clarithromycin", "6 callable / 51 uncallable"),
        ("ZENODO_10369064", "levofloxacin", "#3A9D5D",
         "(b) Paired-read levofloxacin", "42 callable / 4 uncallable"),
    ]
    for axis, (cohort, drug, color, title, sizes) in zip(axes, panels):
        values = associations.loc[
            associations["dataset_id"].eq(cohort) & associations["antibiotic"].eq(drug)
        ].set_index("metric").loc[METRICS, "cliffs_delta_callable_minus_uncallable"]
        axis.barh(np.arange(4), values, height=0.58, color=color)
        axis.set_yticks(np.arange(4), LABELS)
        axis.invert_yaxis()
        axis.axvline(0, color="#444444", linewidth=0.85)
        axis.set_xlim(-1, 1)
        axis.set_xticks([-1, -0.5, 0, 0.5, 1])
        axis.set_title(f"{title}\n{sizes}", loc="left", pad=8)
        axis.spines[["top", "right"]].set_visible(False)
        axis.spines[["left", "bottom"]].set_color("#777777")
        axis.tick_params(axis="y", length=0, pad=7)
        axis.grid(axis="x", color="#DDDDDD", linewidth=0.5, zorder=0)
        axis.set_axisbelow(True)
    axes[-1].set_xlabel("Cliff's delta (callable minus uncallable)")
    fig.savefig(output / "Supplementary_Figure_S1_QC_conditioned.pdf")
    fig.savefig(output / "Supplementary_Figure_S1_QC_conditioned.png", dpi=400)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path("MGEN_Supplementary_Data.xlsx"))
    parser.add_argument("--output", type=Path, default=Path("qc_reanalysis_output"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data = read_samples(args.input)
    table_s2, table_s6 = summaries(data)
    expected = [(60,57,3,6,51), (60,57,3,57,0), (414,414,0,414,0),
                (414,414,0,414,0), (52,46,6,46,0), (52,46,6,42,4)]
    assert list(table_s2[["n_linked","n_qc_pass","n_qc_fail","n_callable_after_qc","n_uncallable_after_qc"]].itertuples(index=False,name=None)) == expected
    table_s2.to_csv(args.output / "Table_S2_QC_accounting.csv", index=False)
    table_s6.to_csv(args.output / "Table_S6_QC_conditioned.csv", index=False)
    data.loc[data["final_qc_status"].eq("PASS")].to_csv(
        args.output / "QC_passing_analysis_records.csv", index=False
    )
    table_s6.dropna(subset=["cliffs_delta_callable_minus_uncallable"]).to_csv(
        args.output / "Figure_S1_source_data.csv", index=False
    )
    draw_figure(table_s6, args.output)
    environment = {
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
        "scipy": scipy.__version__, "scikit_learn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
        "mann_whitney": {"alternative": "two-sided", "method": "auto", "use_continuity": True,
                        "minimum_each_group": 3},
        "effect_direction": "callable minus uncallable",
        "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "input_rows": len(data), "QC_pass_sample_drug_rows": int(data["final_qc_status"].eq("PASS").sum()),
        "scope": "secondary assembly-quality comparison within final-genome-QC-passing records",
        "randomness": "none",
    }
    (args.output / "analysis_environment.json").write_text(json.dumps(environment, indent=2), encoding="utf-8")
    print(table_s2.to_string(index=False))
    print(table_s6.dropna(subset=["mann_whitney_p"]).to_string(index=False))


if __name__ == "__main__":
    main()

