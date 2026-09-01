from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "results/extended_analysis"


def test_end_to_end_yields_and_clarithromycin_bounds() -> None:
    yields = pd.read_csv(ANALYSIS / "coverage_aware_yield.csv")
    ningxia = yields[yields["dataset_id"].eq("CHINA_NINGXIA_2022")].set_index("antibiotic")
    assert int(ningxia.loc["clarithromycin", "correct_n"]) == 6
    assert int(ningxia.loc["clarithromycin", "phenotype_linked_binary_n"]) == 60
    assert ningxia.loc["clarithromycin", "actionable_correct_yield"] == pytest.approx(0.10)
    assert int(ningxia.loc["levofloxacin", "correct_n"]) == 34
    assert ningxia.loc["levofloxacin", "actionable_correct_yield"] == pytest.approx(34 / 60)

    bounds = pd.read_csv(ANALYSIS / "callability_identification_bounds.csv")
    clarithromycin = bounds[
        bounds["dataset_id"].eq("CHINA_NINGXIA_2022")
        & bounds["antibiotic"].eq("clarithromycin")
    ].set_index("metric")
    assert clarithromycin.loc["sensitivity", "logical_lower_bound"] == pytest.approx(2 / 31)
    assert clarithromycin.loc["sensitivity", "logical_upper_bound"] == pytest.approx(1.0)
    assert clarithromycin.loc["specificity", "logical_lower_bound"] == pytest.approx(4 / 29)
    assert clarithromycin.loc["specificity", "logical_upper_bound"] == pytest.approx(1.0)


def test_manifold_abstention_does_not_rescue_ningxia_levofloxacin() -> None:
    table = pd.read_csv(ANALYSIS / "manifold_abstention_metrics.csv")
    eligible = table[
        table["dataset_id"].eq("CHINA_NINGXIA_2022")
        & table["antibiotic"].eq("levofloxacin")
        & table["n_resistant"].ge(10)
        & table["n_susceptible"].ge(10)
    ]
    assert not eligible.empty
    assert eligible["passes_frozen_safety_gate"].eq("no").all()
    assert eligible["false_susceptible_rate"].min() == pytest.approx(0.50)
    best = eligible.loc[eligible["false_susceptible_rate"].idxmin()]
    assert best["retained_coverage"] == pytest.approx(50 / 57)


def test_external_levofloxacin_difference_is_formally_quantified() -> None:
    differences = pd.read_csv(ANALYSIS / "external_performance_differences.csv").set_index("metric")
    sensitivity = differences.loc["sensitivity"]
    assert sensitivity["absolute_difference"] == pytest.approx(12 / 28 - 14 / 16)
    assert sensitivity["difference_ci_high"] < 0
    assert sensitivity["fisher_exact_p"] == pytest.approx(0.004614271335482203)
    specificity = differences.loc["specificity"]
    assert specificity["difference_ci_low"] < 0 < specificity["difference_ci_high"]
    balanced = differences.loc["balanced_accuracy"]
    assert balanced["absolute_difference"] == pytest.approx(0.5935960591 - 0.9157608696, abs=1e-8)
    assert balanced["difference_ci_high"] < 0


def test_frozen_transportability_labels_are_unchanged() -> None:
    labels = pd.read_csv(ROOT / "results/external_validation/transportability_classification.csv")
    observed = dict(zip(labels["antibiotic"], labels["transportability_label"]))
    assert observed == {
        "clarithromycin": "INSUFFICIENT_DATA",
        "levofloxacin": "HIGH_FALSE_SUSCEPTIBLE_RISK",
    }
