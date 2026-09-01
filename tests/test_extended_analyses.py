from pathlib import Path
import math
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_extended_analyses import bh_adjust, compute_metrics, frozen_mutations, parse_missense


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
