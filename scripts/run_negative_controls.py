#!/usr/bin/env python3
"""Run the prespecified phenotype-blind negative controls and robustness checks."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedShuffleSplit

from benchmark_frozen_panels import (
    assembly_predictions, point_metrics, read_csv, summarize, write_csv,
)
from cluster_lineages import allele_feature_matrix
from run_leakage_benchmark import (
    MARKERS, eligible_group, evaluate, feature_table, probability_metrics,
)


SEED = 20260830
REPEATS = 100
EXTERNAL_COHORTS = ("CHINA_NINGXIA_2022", "ZENODO_10369064")


def select_prevalence_matched(
    candidate_prevalence: np.ndarray,
    targets: list[float],
    rng: np.random.Generator,
    tolerance: float = 0.005,
) -> list[int]:
    """Select distinct candidates from the closest phenotype-blind frequency band."""
    available = np.ones(len(candidate_prevalence), dtype=bool)
    selected: list[int] = []
    for target in targets:
        differences = np.abs(candidate_prevalence - target)
        minimum = float(differences[available].min())
        pool = np.flatnonzero(available & (differences <= minimum + tolerance))
        chosen = int(rng.choice(pool))
        selected.append(chosen)
        available[chosen] = False
    return selected


def clone_thinned_metrics(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    eligible = [row for row in rows if row["analysis_status"] == "PRIMARY"]
    retained: list[dict[str, str]] = []
    groups: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in eligible:
        unit = row["near_clone_group"] or f"{row['dataset_id']}::{row['isolate_id']}"
        groups[(row["dataset_id"], row["antibiotic"], unit)].append(row)
    for values in groups.values():
        retained.append(min(values, key=lambda row: (row["dataset_id"], row["isolate_id"])))
    by_cohort: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in retained:
        by_cohort[(row["dataset_id"], row["antibiotic"])].append(row)
    return [
        summarize(values, {
            "scenario": "LEXICOGRAPHIC_FIRST_PER_NEAR_CLONE_COMPONENT",
            "stratum_type": "COHORT", "stratum": dataset,
            "dataset_id": dataset, "antibiotic": antibiotic,
        })
        for (dataset, antibiotic), values in sorted(by_cohort.items())
    ]


def raw_assembly_concordance(root: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    raw = {
        (row["dataset_id"], row["isolate_id"], row["antibiotic"]): row
        for row in read_csv(root / "results/panels/zenodo_read_marker_predictions.csv")
    }
    support = {
        (row["dataset_id"], row["isolate_id"], row["antibiotic"]): row
        for row in assembly_predictions(root, include_zenodo=True)
        if row["dataset_id"] == "ZENODO_10369064"
    }
    sample_rows: list[dict[str, object]] = []
    for key, read_row in sorted(raw.items()):
        assembly_row = support.get(key, {})
        raw_call = read_row["prediction"]
        assembly_call = assembly_row.get("prediction", "MISSING")
        both = raw_call in {"S", "R"} and assembly_call in {"S", "R"}
        sample_rows.append({
            "dataset_id": key[0], "isolate_id": key[1], "antibiotic": key[2],
            "raw_read_prediction": raw_call, "support_assembly_prediction": assembly_call,
            "raw_read_callable": "yes" if raw_call in {"S", "R"} else "no",
            "support_assembly_callable": "yes" if assembly_call in {"S", "R"} else "no",
            "both_callable": "yes" if both else "no",
            "concordant_when_both_callable": "yes" if both and raw_call == assembly_call else "no" if both else "",
            "interpretation": "raw-read call is authoritative; assembly is support-only",
        })
    summary: list[dict[str, object]] = []
    for antibiotic in sorted({row["antibiotic"] for row in sample_rows}):
        values = [row for row in sample_rows if row["antibiotic"] == antibiotic]
        both = [row for row in values if row["both_callable"] == "yes"]
        agreements = sum(row["concordant_when_both_callable"] == "yes" for row in both)
        summary.append({
            "dataset_id": "ZENODO_10369064", "antibiotic": antibiotic,
            "n": len(values),
            "n_raw_callable": sum(row["raw_read_callable"] == "yes" for row in values),
            "n_support_assembly_callable": sum(row["support_assembly_callable"] == "yes" for row in values),
            "n_both_callable": len(both), "n_concordant": agreements,
            "n_discordant": len(both) - agreements,
            "concordance_when_both_callable": agreements / len(both) if both else np.nan,
        })
    return sample_rows, summary


def label_permutation_controls(table: pd.DataFrame) -> list[dict[str, object]]:
    rng = np.random.default_rng(SEED)
    output: list[dict[str, object]] = []
    for drug in ("clarithromycin", "levofloxacin"):
        development = table[
            (table["dataset_id"] == "HPGP_GLOBAL") & (table["antibiotic"] == drug)
        ].reset_index(drop=True)
        splitter = StratifiedShuffleSplit(n_splits=REPEATS, test_size=0.20, random_state=SEED)
        for repeat, (train_index, test_index) in enumerate(
            splitter.split(development, development["outcome"]), start=1
        ):
            train = development.iloc[train_index].copy()
            train["outcome"] = rng.permutation(train["outcome"].to_numpy())
            result = evaluate(
                train, development.iloc[test_index], drug, "MUTATION_ONLY_LOGISTIC",
                "RANDOM_ISOLATE_LABEL_PERMUTATION", str(repeat),
            )
            output.append({"negative_control": "TRAINING_LABEL_PERMUTATION", **result})

            external_train = development.copy()
            external_train["outcome"] = rng.permutation(external_train["outcome"].to_numpy())
            for external in EXTERNAL_COHORTS:
                test = table[(table["dataset_id"] == external) & (table["antibiotic"] == drug)]
                if not eligible_group(test):
                    continue
                result = evaluate(
                    external_train, test, drug, "MUTATION_ONLY_LOGISTIC",
                    "HPGP_TO_EXTERNAL_LABEL_PERMUTATION", f"{external}__{repeat:03d}",
                )
                output.append({"negative_control": "TRAINING_LABEL_PERMUTATION", **result})
    return output


def matrix_evaluation(
    train_x, train_y: np.ndarray, test_x, test_y: np.ndarray,
    drug: str, split_type: str, fold: str,
) -> dict[str, object]:
    model = LogisticRegression(
        C=1.0, solver="liblinear", random_state=SEED, max_iter=1000,
    )
    model.fit(train_x, train_y)
    probabilities = model.predict_proba(test_x)[:, 1]
    predictions = (probabilities >= 0.5).astype(int)
    return {
        "negative_control": "SIZE_PREVALENCE_MATCHED_RANDOM_CORE_SNP_PANEL",
        "antibiotic": drug, "model": "RANDOM_CORE_SNP_LOGISTIC",
        "split_type": split_type, "fold": fold,
        "n_train": len(train_y), "n_test": len(test_y),
        **point_metrics(test_y.tolist(), predictions.tolist()),
        **probability_metrics(test_y.tolist(), probabilities.tolist()),
    }


def random_snp_controls(
    root: Path, table: pd.DataFrame,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    names, matrix, features = allele_feature_matrix(
        root / "data/phylogeny/core/core.aln", training_dataset="HPGP_GLOBAL"
    )
    row_by_name = {name: index for index, name in enumerate(names)}
    hp_rows = np.asarray([name.startswith("HPGP_GLOBAL__") for name in names], dtype=bool)
    candidate_prevalence = np.asarray(matrix[hp_rows].mean(axis=0)).ravel()
    rng = np.random.default_rng(SEED)
    output: list[dict[str, object]] = []
    trace: list[dict[str, object]] = []

    def alignment_indices(frame: pd.DataFrame) -> np.ndarray:
        identifiers = [f"{row.dataset_id}__{row.isolate_id}" for row in frame.itertuples()]
        missing = [identifier for identifier in identifiers if identifier not in row_by_name]
        if missing:
            raise KeyError(f"{len(missing)} phenotype rows absent from core alignment; first={missing[0]}")
        return np.asarray([row_by_name[identifier] for identifier in identifiers], dtype=int)

    for drug in ("clarithromycin", "levofloxacin"):
        development = table[
            (table["dataset_id"] == "HPGP_GLOBAL") & (table["antibiotic"] == drug)
        ].reset_index(drop=True)
        development_rows = alignment_indices(development)
        target_prevalence = [float(development[column].mean()) for column in MARKERS[drug]]
        splitter = StratifiedShuffleSplit(n_splits=REPEATS, test_size=0.20, random_state=SEED)
        splits = list(splitter.split(development, development["outcome"]))
        for repeat in range(1, REPEATS + 1):
            selected = select_prevalence_matched(candidate_prevalence, target_prevalence, rng)
            for marker, target, index in zip(MARKERS[drug], target_prevalence, selected):
                position, allele = features[index]
                trace.append({
                    "antibiotic": drug, "repeat": repeat, "matched_panel_marker": marker,
                    "core_alignment_position_1_based": position + 1,
                    "alternate_allele": chr(allele),
                    "target_marker_prevalence": target,
                    "selected_snp_prevalence": float(candidate_prevalence[index]),
                    "absolute_prevalence_difference": abs(float(candidate_prevalence[index]) - target),
                })
            train_index, test_index = splits[repeat - 1]
            selected_matrix = matrix[:, selected]
            output.append(matrix_evaluation(
                selected_matrix[development_rows[train_index]], development.iloc[train_index]["outcome"].to_numpy(),
                selected_matrix[development_rows[test_index]], development.iloc[test_index]["outcome"].to_numpy(),
                drug, "RANDOM_ISOLATE_RANDOM_SNP_PANEL", str(repeat),
            ))
            for external in EXTERNAL_COHORTS:
                test = table[(table["dataset_id"] == external) & (table["antibiotic"] == drug)].reset_index(drop=True)
                if not eligible_group(test):
                    continue
                test_rows = alignment_indices(test)
                output.append(matrix_evaluation(
                    selected_matrix[development_rows], development["outcome"].to_numpy(),
                    selected_matrix[test_rows], test["outcome"].to_numpy(),
                    drug, "HPGP_TO_EXTERNAL_RANDOM_SNP_PANEL", f"{external}__{repeat:03d}",
                ))
    return output, trace


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.root.resolve()
    output = root / "results/negative_controls"
    output.mkdir(parents=True, exist_ok=True)
    sample_rows = read_csv(root / "results/external_validation/sample_level_predictions.csv")
    write_csv(output / "clone_thinned_panel_metrics.csv", clone_thinned_metrics(sample_rows))
    concordance_samples, concordance_summary = raw_assembly_concordance(root)
    write_csv(output / "raw_assembly_concordance_samples.csv", concordance_samples)
    write_csv(output / "raw_assembly_concordance_summary.csv", concordance_summary)
    table = feature_table(root)
    write_csv(output / "label_permutation_metrics.csv", label_permutation_controls(table))
    random_metrics, random_trace = random_snp_controls(root, table)
    write_csv(output / "random_snp_panel_metrics.csv", random_metrics)
    write_csv(output / "random_snp_panel_trace.csv", random_trace)
    status = [
        {"control_id": 1, "control": "random_vs_leave_study_out", "status": "COMPUTED_ELSEWHERE", "output": "leakage_benchmark_summary.csv", "reason": "random HpGP splits versus HpGP-to-external tests"},
        {"control_id": 2, "control": "random_vs_leave_lineage_out", "status": "COMPUTED_ELSEWHERE", "output": "leakage_benchmark_summary.csv", "reason": "same frozen-feature models"},
        {"control_id": 3, "control": "exclude_near_clones", "status": "COMPUTED", "output": "clone_thinned_panel_metrics.csv", "reason": "lexicographic first isolate per connected component"},
        {"control_id": 4, "control": "one_isolate_per_patient", "status": "NOT_ESTIMABLE", "output": "", "reason": "public patient identifiers unavailable"},
        {"control_id": 5, "control": "raw_read_supported_only", "status": "COMPUTED_ELSEWHERE", "output": "phenotype_sensitivity_metrics.csv", "reason": "Zenodo raw-read predictions authoritative"},
        {"control_id": 6, "control": "assembly_supported_only", "status": "COMPUTED_ELSEWHERE", "output": "phenotype_sensitivity_metrics.csv", "reason": "HpGP and Ningxia assembly predictions; Zenodo support assembly separate"},
        {"control_id": 7, "control": "alternative_breakpoint_interpretation", "status": "COMPUTED_ELSEWHERE", "output": "phenotype_sensitivity_metrics.csv", "reason": "original and recomputed Ningxia labels with I sensitivities"},
        {"control_id": 8, "control": "exclude_borderline_MIC", "status": "COMPUTED_ELSEWHERE", "output": "phenotype_sensitivity_metrics.csv", "reason": "prespecified borderline flag"},
        {"control_id": 9, "control": "AST_method_strata", "status": "COMPUTED_ELSEWHERE", "output": "phenotype_sensitivity_metrics.csv", "reason": "method-specific descriptive estimates"},
        {"control_id": 10, "control": "unified_laboratory_AST", "status": "COMPUTED_ELSEWHERE", "output": "frozen_panel_metrics.csv", "reason": "cohort-specific estimates preserve laboratory/study separation"},
        {"control_id": 11, "control": "high_quality_genomes_only", "status": "COMPUTED_ELSEWHERE", "output": "phenotype_sensitivity_metrics.csv", "reason": "frozen complete genome-QC pass"},
        {"control_id": 12, "control": "country_strata", "status": "COMPUTED_ELSEWHERE", "output": "frozen_panel_metrics.csv", "reason": "country-specific descriptive estimates"},
        {"control_id": 13, "control": "lineage_strata", "status": "COMPUTED_ELSEWHERE", "output": "frozen_panel_metrics.csv", "reason": "predeclared estimability gate plus descriptive sensitivity table"},
        {"control_id": 14, "control": "lineage_only_model", "status": "COMPUTED_ELSEWHERE", "output": "leakage_benchmark_summary.csv", "reason": "prespecified negative-control model"},
        {"control_id": 15, "control": "random_core_SNP_panels_100", "status": "COMPUTED", "output": "random_snp_panel_metrics.csv", "reason": "HpGP sequence frequency matching without phenotype selection"},
        {"control_id": 16, "control": "training_label_permutation_100", "status": "COMPUTED", "output": "label_permutation_metrics.csv", "reason": "HpGP training labels only"},
    ]
    write_csv(output / "negative_control_status.csv", status)


if __name__ == "__main__":
    main()
