import numpy as np
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_negative_controls import select_prevalence_matched


def test_prevalence_matching_is_distinct_and_near_targets() -> None:
    prevalence = np.asarray([0.05, 0.10, 0.11, 0.49, 0.50, 0.51])
    selected = select_prevalence_matched(prevalence, [0.10, 0.50], np.random.default_rng(7))
    assert len(selected) == len(set(selected)) == 2
    assert abs(prevalence[selected[0]] - 0.10) <= 0.005
    assert abs(prevalence[selected[1]] - 0.50) <= 0.005
