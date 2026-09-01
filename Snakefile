from pathlib import Path
import json

DECISION = Path("reports/phase0_decision.json")

def phase0_authorized():
    if not DECISION.exists():
        return False
    payload = json.loads(DECISION.read_text(encoding="utf-8"))
    return payload.get("decision", "").startswith(("GO_", "CONDITIONAL_GO_"))

rule all:
    input:
        "reports/phase0_decision.json"

rule phase0_guard:
    input:
        DECISION
    output:
        touch("results/logs/phase0_guard.ok")
    shell:
        "/usr/bin/python3 scripts/check_phase0_decision.py {input} && touch {output}"


rule phase2:
    input:
        "results/external_validation/transportability_classification.json",
        "results/lineage_validation/leakage_benchmark_summary.csv",
        "results/negative_controls/negative_control_status.csv",
        "results/sensitivity/phenotype_sensitivity_metrics.csv",
        "results/mic/ningxia_mic_lower_bound_summary.csv",


rule acquire_phase2_assemblies:
    input:
        guard="results/logs/phase0_guard.ok",
        isolates="metadata/isolate_manifest.csv",
        index="metadata/phase0/ncbi_assembly_index.csv",
    output:
        "metadata/phase2/acquisition_assemblies.csv",
    log:
        "results/logs/phase2_acquire_assemblies.snakemake.log",
    conda:
        "environment.yml"
    shell:
        "python scripts/acquire_phase2.py --scope assemblies --workers 6 --output {output} > {log} 2>&1"


rule acquire_phase2_reads:
    input:
        guard="results/logs/phase0_guard.ok",
        isolates="metadata/isolate_manifest.csv",
        record="metadata/raw/zenodo_10369064_2026-08-30.json",
    output:
        "metadata/phase2/acquisition_reads.csv",
    log:
        "results/logs/phase2_acquire_reads.snakemake.log",
    conda:
        "environment.yml"
    shell:
        "python scripts/acquire_phase2.py --scope reads --workers 8 --output {output} > {log} 2>&1"


rule process_zenodo_reads:
    input:
        "metadata/phase2/acquisition_reads.csv",
    output:
        status="results/qc/zenodo_processing_status.csv",
        predictions="results/panels/zenodo_read_marker_predictions.csv",
    log:
        "results/logs/phase2_process_zenodo_reads.log",
    conda:
        "environment.yml"
    shell:
        "python scripts/process_zenodo_reads.py --workers 2 --threads-per-sample 3 > {log} 2>&1"


rule process_phase2_assemblies:
    input:
        assemblies="metadata/phase2/acquisition_assemblies.csv",
        zenodo="results/qc/zenodo_processing_status.csv",
    output:
        "results/qc/assembly_qc.csv",
    log:
        "results/logs/phase2_assembly_qc.log",
    conda:
        "environment.yml"
    shell:
        "python scripts/process_phase2_assemblies.py --workers 8 > {log} 2>&1"


rule call_phase2_assembly_markers:
    input:
        "results/qc/assembly_qc.csv",
    output:
        "results/panels/assembly_marker_calls.csv",
    log:
        "results/logs/phase2_assembly_markers.log",
    conda:
        "environment.yml"
    shell:
        "python scripts/call_phase2_assembly_markers.py --workers 6 > {log} 2>&1"


rule checkm2_database:
    output:
        "data/raw/checkm2_db/complete.json",
    log:
        "results/logs/checkm2_database.log",
    shell:
        "python scripts/acquire_checkm2_database.py > {log} 2>&1"


rule phase2_checkm2:
    input:
        qc="results/qc/assembly_qc.csv",
        database="data/raw/checkm2_db/complete.json",
    output:
        "results/qc/assembly_qc_with_checkm2.csv",
    log:
        "results/logs/checkm2_wrapper.log",
    conda:
        "workflow/envs/checkm2.yaml"
    shell:
        "ln -sfn \"$CONDA_PREFIX\" /tmp/hpylori-checkm2-env && /tmp/hpylori-checkm2-env/bin/python scripts/run_checkm2.py --database data/raw/checkm2_db --threads 6 > {log} 2>&1"


rule full_relatedness:
    input:
        "results/qc/assembly_qc_with_checkm2.csv",
    output:
        pairs="results/qc/pairwise_relatedness_candidates.csv",
        groups="results/qc/near_clone_groups.csv",
    log:
        "results/logs/full_relatedness.log",
    conda:
        "environment.yml"
    shell:
        "python scripts/build_relatedness.py --threads 6 > {log} 2>&1"


rule core_snp:
    input:
        "results/qc/assembly_qc_with_checkm2.csv",
    output:
        alignment="data/phylogeny/core/core.aln",
        status="results/qc/core_snp_processing_status.csv",
    log:
        "results/logs/core_snp_wrapper.log",
    conda:
        "workflow/envs/snippy.yaml"
    shell:
        "python scripts/prepare_snippy_runtime.py --source-prefix \"$CONDA_PREFIX\" --runtime /tmp/hpylori-snippy-runtime && python scripts/run_core_snp.py --snippy /tmp/hpylori-snippy-runtime/bin/snippy --snippy-core /tmp/hpylori-snippy-runtime/bin/snippy-core --tool-prefix \"$CONDA_PREFIX\" --workers 4 > {log} 2>&1"


rule cluster_lineages:
    input:
        "data/phylogeny/core/core.aln",
    output:
        "results/lineage_validation/lineage_assignments.csv",
    conda:
        "environment.yml"
    shell:
        "python scripts/cluster_lineages.py"


rule benchmark_frozen_panels:
    input:
        assembly_calls="results/panels/assembly_marker_calls.csv",
        read_calls="results/panels/zenodo_read_marker_predictions.csv",
        qc="results/qc/assembly_qc_with_checkm2.csv",
        clones="results/qc/near_clone_groups.csv",
        lineages="results/lineage_validation/lineage_assignments.csv",
    output:
        samples="results/external_validation/sample_level_predictions.csv",
        metrics="results/external_validation/frozen_panel_metrics.csv",
    conda:
        "environment.yml"
    shell:
        "python scripts/benchmark_frozen_panels.py"


rule sensitivity_analyses:
    input:
        "results/external_validation/sample_level_predictions.csv",
    output:
        phenotype="results/sensitivity/phenotype_sensitivity_metrics.csv",
        mic="results/mic/ningxia_mic_lower_bound_summary.csv",
        callability="results/qc/panel_callability.csv",
        boundaries="results/external_validation/secondary_drug_boundaries.csv",
    conda:
        "environment.yml"
    shell:
        "python scripts/run_sensitivity_analyses.py"


rule classify_transportability:
    input:
        "results/external_validation/frozen_panel_metrics.csv",
    output:
        csv="results/external_validation/transportability_classification.csv",
        json="results/external_validation/transportability_classification.json",
    conda:
        "environment.yml"
    shell:
        "python scripts/classify_transportability.py"


rule leakage_benchmark:
    input:
        samples="results/external_validation/sample_level_predictions.csv",
        calls="results/panels/assembly_marker_calls.csv",
        lineages="results/lineage_validation/lineage_assignments.csv",
        clones="results/qc/near_clone_groups.csv",
    output:
        folds="results/lineage_validation/leakage_benchmark_folds.csv",
        summary="results/lineage_validation/leakage_benchmark_summary.csv",
    conda:
        "environment.yml"
    shell:
        "python scripts/run_leakage_benchmark.py"


rule negative_controls:
    input:
        samples="results/external_validation/sample_level_predictions.csv",
        alignment="data/phylogeny/core/core.aln",
        calls="results/panels/assembly_marker_calls.csv",
        lineages="results/lineage_validation/lineage_assignments.csv",
        clones="results/qc/near_clone_groups.csv",
    output:
        status="results/negative_controls/negative_control_status.csv",
        clone="results/negative_controls/clone_thinned_panel_metrics.csv",
        concordance="results/negative_controls/raw_assembly_concordance_summary.csv",
        permutation="results/negative_controls/label_permutation_metrics.csv",
        random_panels="results/negative_controls/random_snp_panel_metrics.csv",
        random_trace="results/negative_controls/random_snp_panel_trace.csv",
    conda:
        "environment.yml"
    shell:
        "python scripts/run_negative_controls.py"


rule transport_shift_analyses:
    input:
        samples="results/external_validation/sample_level_predictions.csv",
        lineages="results/lineage_validation/lineage_assignments.csv",
        manifold="results/extended_analysis/development_manifold_distance.csv",
        audit="results/extended_analysis/ningxia_lvx_error_audit_enriched.csv",
        callability="results/qc/panel_callability.csv",
    output:
        summary="results/transport_shift/three_layer_transport_summary.csv",
        performance="results/transport_shift/levofloxacin_lineage_performance.csv",
        comparisons="results/transport_shift/dominant_lineage_comparisons.csv",
        standardization="results/transport_shift/lineage_standardized_performance.csv",
        discriminator="results/transport_shift/cohort_discriminator_summary.csv",
        discriminator_null="results/transport_shift/cohort_discriminator_permutations.csv",
        mic="results/transport_shift/ningxia_resistant_mic_severity.csv",
        atlas="results/transport_shift/ningxia_failure_atlas.csv",
    conda:
        "environment.yml"
    shell:
        "python scripts/run_transport_shift_analyses.py --root . --permutations 1000 --bootstrap-repeats 2000"


rule transport_shift_figures:
    input:
        "results/transport_shift/three_layer_transport_summary.csv",
        "results/transport_shift/cohort_discriminator_summary.csv",
        "results/transport_shift/levofloxacin_lineage_performance.csv",
        "results/transport_shift/ningxia_failure_atlas.csv",
    output:
        main="figures/main/Figure9_three_layer_transport_shift.png",
        s4="figures/supplementary/Supplementary_Figure_S4_lineage_conditional_transport.png",
        s5="figures/supplementary/Supplementary_Figure_S5_adversarial_cohort_shift.png",
        s6="figures/supplementary/Supplementary_Figure_S6_Ningxia_failure_atlas.png",
    conda:
        "environment.yml"
    shell:
        "python scripts/make_transport_shift_figures.py --root ."
