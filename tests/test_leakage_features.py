from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_leakage_benchmark import raw_features


def test_raw_clarithromycin_features_use_frozen_mixed_threshold():
    features = raw_features({"antibiotic": "clarithromycin", "marker_summary": "A2142:0.1999;A2143:0.2000"})
    assert features == {"CLR_A2142_ANY": 0, "CLR_A2143G": 1}


def test_raw_levofloxacin_features_are_exact_catalogue_changes():
    features = raw_features({"antibiotic": "levofloxacin", "marker_summary": "N87K;D91D"})
    assert features["LVX_N87K"] == 1
    assert features["LVX_D91G"] == 0
