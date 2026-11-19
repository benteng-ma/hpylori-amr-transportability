from pathlib import Path
import math
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_extended_analyses import bh_adjust, compute_metrics, frozen_mutations, parse_missense, qc_callability


def test_parse_missense_normalizes_three_letter_amino_acids() -> None:
    assert parse_missense("missense_variant c.259A>T p.Asn87Tyr") == "N87Y"
    assert parse_missense("missense_variant p.Ter100Gln") == "*100Q"
    assert parse_missense("synonymous_variant p.Asn87Asn") == "N87N"
    assert parse_missense("upstream_gene_variant") is None


def test_frozen_mutations_never_promotes_off_panel_variant() -> None:
    assert frozen_mutations("levofloxacin", "N87K;N87Y;D91N") == ["N87K", "D91N"]
    assert frozen_mutations("clarithromycin", "A2142:0.20;A2143:0.19") == ["A2142 resistant allele"]
    assert frozen_mutations("clarithromycin", "A2142:0.19;A2143:0.20") == ["A2143G"]


def test_compute_metrics_recovers_false_susceptible_rate() -> None:
    frame = pd.DataFrame(
        {
            "phenotype": ["R", "R", "R", "S", "S", "S"],
            "prediction": ["R", "S", "S", "S", "S", "R"],
        }
    )
    metrics = compute_metrics(frame)
    assert (metrics["tp"], metrics["tn"], metrics["fp"], metrics["fn"]) == (1, 2, 1, 2)
    assert math.isclose(metrics["sensitivity"], 1 / 3)
    assert math.isclose(metrics["specificity"], 2 / 3)
    assert math.isclose(metrics["false_susceptible_rate"], 2 / 3)
    assert math.isclose(metrics["balanced_accuracy"], 0.5)


def test_bh_adjust_preserves_order_and_monotonicity() -> None:
    adjusted = bh_adjust([0.01, 0.04, 0.03])
    assert adjusted == [0.03, 0.04, 0.04]


def test_qc_callability_excludes_genome_qc_failures_from_associations(tmp_path) -> None:
    """QC failures stay in sample exports and do not become target-failure controls."""
    isolate_ids = [f"isolate_{i}" for i in range(7)]
    predictions = pd.DataFrame(
        {
            "dataset_id": ["test_cohort"] * 7,
            "isolate_id": isolate_ids,
            "antibiotic": ["clarithromycin"] * 7,
            "callable": ["yes"] * 3 + ["no"] * 4,
        }
    )
    quality = pd.DataFrame(
        {
            "dataset_id": ["test_cohort"] * 7,
            "isolate_id": isolate_ids,
            "final_qc_status": ["PASS"] * 6 + ["FAIL"],
            "contigs": [10, 11, 12, 20, 21, 22, 99999],
            "assembly_size_bp": [1600000] * 7,
            "n50_bp": [100000] * 7,
            "longest_contig_bp": [200000] * 7,
            "gc_percent": [39] * 7,
            "ambiguous_bases": [0] * 7,
            "completeness_percent": [99] * 6 + [50],
            "contamination_percent": [0.1] * 6 + [20],
        }
    )

    qc_callability(predictions, quality, tmp_path)

    samples = pd.read_csv(tmp_path / "qc_callability_samples.csv")
    associations = pd.read_csv(tmp_path / "qc_callability_associations.csv")
    assert len(samples) == 7
    assert samples["final_qc_status"].eq("FAIL").sum() == 1
    assert associations["n_callable"].eq(3).all()
    assert associations["n_uncallable"].eq(3).all()
    contigs = associations.loc[associations["metric"].eq("contigs")].iloc[0]
    assert contigs["callable_median"] == 11
    assert contigs["uncallable_median"] == 21
    assert contigs["cliffs_delta_callable_minus_uncallable"] == -1
    assert contigs["auc_metric_predicting_callable"] == 0
