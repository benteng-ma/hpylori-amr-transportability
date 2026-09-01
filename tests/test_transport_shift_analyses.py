import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "transport_shift", ROOT / "scripts" / "run_transport_shift_analyses.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_exact_interval_has_defined_boundaries():
    assert MODULE.exact_interval(0, 10)[0] == 0.0
    assert MODULE.exact_interval(10, 10)[1] == 1.0


def test_lineage_performance_keeps_false_susceptible_denominator_resistant_only():
    frame = pd.DataFrame(
        {
            "dataset_id": ["A"] * 4,
            "lineage_recomputed": ["L1"] * 4,
            "phenotype": ["R", "R", "S", "S"],
            "prediction": ["R", "S", "R", "S"],
        }
    )
    result = MODULE.lineage_performance(MODULE.add_binary_columns(frame))
    row = result[
        (result["lineage_recomputed"] == "L1")
        & (result["metric"] == "false_susceptible_rate")
    ].iloc[0]
    assert row["successes"] == 1
    assert row["total"] == 2
    assert row["estimate"] == 0.5


def test_standardization_uses_target_phenotype_specific_lineage_weights():
    source = pd.DataFrame(
        {
            "phenotype": ["R", "R", "R", "R"],
            "prediction": ["R", "R", "R", "S"],
            "lineage_recomputed": ["L1", "L1", "L2", "L2"],
        }
    )
    target = pd.DataFrame(
        {
            "phenotype": ["R", "R", "R", "R"],
            "prediction": ["R", "R", "R", "R"],
            "lineage_recomputed": ["L1", "L1", "L1", "L2"],
        }
    )
    estimate, pieces = MODULE.standardized_rate(source, target, "sensitivity")
    assert estimate == 0.875
    assert sum(item["target_weight"] for item in pieces) == 1.0


def test_primary_filter_excludes_uncallable_and_nonprimary_rows():
    frame = pd.DataFrame(
        {
            "antibiotic": ["levofloxacin"] * 3,
            "analysis_status": ["PRIMARY", "EXCLUDED_AMBIGUOUS_OR_UNCALLABLE", "PRIMARY"],
            "phenotype": ["R", "R", "I"],
            "prediction": ["S", "UNCALLABLE", "S"],
        }
    )
    result = MODULE.primary_levofloxacin(frame)
    assert len(result) == 1
    assert result.iloc[0]["outcome"] == "FN"
