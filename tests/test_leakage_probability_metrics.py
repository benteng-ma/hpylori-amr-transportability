import math

from run_leakage_benchmark import probability_metrics


def test_probability_metrics_use_probabilities_not_hard_predictions() -> None:
    metrics = probability_metrics([0, 0, 1, 1], [0.1, 0.4, 0.6, 0.9])
    assert metrics["auroc_probability"] == 1.0
    assert metrics["auprc_probability"] == 1.0
    assert math.isclose(metrics["brier_probability"], 0.085)
    assert math.isfinite(metrics["calibration_intercept"])
    assert math.isfinite(metrics["calibration_slope"])
