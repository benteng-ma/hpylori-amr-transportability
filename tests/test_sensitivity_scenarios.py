from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_sensitivity_analyses import scenario_label


def test_recomputed_intermediate_scenarios_do_not_modify_original_label() -> None:
    row = {"phenotype": "S", "phenotype_recomputed": "I", "borderline_mic": "yes"}
    assert scenario_label(row, "PRIMARY_ORIGINAL_I_EXCLUDED") == "S"
    assert scenario_label(row, "RECOMPUTED_I_EXCLUDED") == ""
    assert scenario_label(row, "RECOMPUTED_I_AS_S") == "S"
    assert scenario_label(row, "RECOMPUTED_I_AS_R") == "R"
    assert scenario_label(row, "EXCLUDE_BORDERLINE_MIC") == ""
