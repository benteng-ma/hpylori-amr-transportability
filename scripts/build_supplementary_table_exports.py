#!/usr/bin/env python3
"""Consolidate the machine-readable result inventory into Tables S1-S20."""

from __future__ import annotations

import csv
from pathlib import Path


TABLES = {
    "S01_Crosswalk": ["results/external_validation/isolate_sequence_phenotype_crosswalk.csv"],
    "S02_QC_Callability": ["results/qc/assembly_qc_with_checkm2.csv", "results/qc/panel_callability.csv"],
    "S03_Lineage_Relatedness": ["results/qc/near_clone_groups.csv", "results/qc/pairwise_relatedness_candidates.csv", "results/lineage_validation/lineage_assignments.csv"],
    "S04_Predictions": ["results/external_validation/sample_level_predictions.csv"],
    "S05_Performance": ["results/external_validation/frozen_panel_metrics.csv"],
    "S06_Validation_Design": ["results/lineage_validation/leakage_benchmark_folds.csv", "results/lineage_validation/leakage_benchmark_summary.csv"],
    "S07_Phenotype_Sensitivity": ["results/sensitivity/phenotype_sensitivity_metrics.csv"],
    "S08_Negative_Controls": [
        "results/negative_controls/negative_control_status.csv",
        "results/negative_controls/clone_thinned_panel_metrics.csv",
        "results/negative_controls/raw_assembly_concordance_summary.csv",
        "results/negative_controls/random_snp_panel_metrics.csv",
        "results/negative_controls/label_permutation_metrics.csv",
    ],
    "S09_MIC_Summary": ["results/mic/ningxia_mic_lower_bound_summary.csv"],
    "S10_Genotype_Audit": ["results/external_validation/ningxia_published_genotype_audit_samples.csv", "results/external_validation/ningxia_published_genotype_audit_summary.csv"],
    "S11_Classification": ["results/external_validation/transportability_classification.csv"],
    "S12_Figure_Manifest": ["results/source_data/main_figure_manifest.csv"],
    "S13_Marker_Architecture": [
        "results/extended_analysis/mutation_spectrum_samples.csv",
        "results/extended_analysis/mutation_spectrum_summary.csv",
        "results/extended_analysis/marker_combination_summary.csv",
        "results/extended_analysis/marker_phenotype_prevalence_shift.csv",
        "results/extended_analysis/marker_prevalence_by_lineage.csv",
    ],
    "S14_QC_Associations": ["results/extended_analysis/qc_callability_samples.csv", "results/extended_analysis/qc_callability_associations.csv"],
    "S15_Population_Shift": [
        "results/extended_analysis/development_manifold_distance.csv",
        "results/extended_analysis/development_manifold_summary.csv",
        "results/extended_analysis/external_error_manifold_samples.csv",
        "results/extended_analysis/external_error_manifold_associations.csv",
    ],
    "S16_Residual_Mechanisms": [
        "results/extended_analysis/ningxia_lvx_error_audit_enriched.csv",
        "results/extended_analysis/ningxia_lvx_false_susceptible_decomposition.csv",
        "results/extended_analysis/ningxia_gyrA_missense_samples.csv",
        "results/extended_analysis/ningxia_gyrA_missense_associations.csv",
        "results/extended_analysis/exploratory_N87Y_rescue.csv",
    ],
    "S17_Stability_Scenarios": [
        "results/extended_analysis/external_bootstrap_distributions.csv",
        "results/extended_analysis/ningxia_lvx_leave_one_out_influence.csv",
        "results/extended_analysis/predictive_values_by_assumed_prevalence.csv",
        "results/extended_analysis/transport_gate_robustness_grid.csv",
    ],
    "S18_Transport_Shift": [
        "results/transport_shift/three_layer_transport_summary.csv",
        "results/transport_shift/levofloxacin_lineage_performance.csv",
        "results/transport_shift/dominant_lineage_comparisons.csv",
        "results/transport_shift/lineage_standardization_weights.csv",
        "results/transport_shift/lineage_standardized_performance.csv",
        "results/transport_shift/cohort_discriminator_summary.csv",
        "results/transport_shift/cohort_discriminator_coefficients.csv",
        "results/transport_shift/cohort_discriminator_permutations.csv",
        "results/transport_shift/ningxia_resistant_mic_summary.csv",
        "results/transport_shift/ningxia_failure_atlas.csv",
    ],
    "S19_23S_Sensitivity": [
        "results/extended_analysis/23s_callability_sensitivity_summary.csv",
        "results/extended_analysis/23s_callability_sensitivity_grid.csv",
        "results/extended_analysis/23s_callability_sensitivity_samples.csv",
    ],
    "S20_Coverage_Selective": [
        "results/extended_analysis/coverage_aware_yield.csv",
        "results/extended_analysis/callability_identification_bounds.csv",
        "results/extended_analysis/manifold_abstention_metrics.csv",
        "results/extended_analysis/external_performance_differences.csv",
        "results/extended_analysis/coverage_selective_summary.csv",
    ],
}


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    output_dir = root / "results/supplementary_tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for table, relative_paths in TABLES.items():
        fields = ["source_file"]
        source_rows: list[tuple[str, list[dict[str, str]]]] = []
        for relative in relative_paths:
            path = root / relative
            source_fields, rows = read_csv(path)
            for field in source_fields:
                if field not in fields:
                    fields.append(field)
            source_rows.append((relative, rows))
            manifest_rows.append({"table": table, "source_file": relative, "n_rows": len(rows)})
        output = output_dir / f"{table}.csv"
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for relative, rows in source_rows:
                for row in rows:
                    writer.writerow({"source_file": relative, **row})
    with (output_dir / "supplementary_table_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["table", "source_file", "n_rows"])
        writer.writeheader()
        writer.writerows(manifest_rows)


if __name__ == "__main__":
    main()
