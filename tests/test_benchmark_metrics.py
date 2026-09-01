from pathlib import Path
import csv
import math
import sys

from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from benchmark_frozen_panels import exact_interval, point_metrics, raw_read_predictions


def test_toy_metrics_recover_expected_confusion():
    metrics = point_metrics([1, 1, 1, 0, 0, 0, 0], [1, 1, 0, 0, 0, 1, 0])
    assert (metrics["tp"], metrics["tn"], metrics["fp"], metrics["fn"]) == (2, 3, 1, 1)
    assert math.isclose(metrics["sensitivity"], 2 / 3)
    assert math.isclose(metrics["specificity"], 3 / 4)


def test_closed_form_binary_scores_match_sklearn():
    cases = [
        ([1, 1, 1, 0, 0, 0, 0], [1, 1, 0, 0, 0, 1, 0]),
        ([1, 1, 0, 0], [0, 0, 0, 0]),
        ([1, 1, 0, 0], [1, 1, 1, 1]),
        ([1, 0, 1, 0, 1, 0], [1, 0, 1, 0, 1, 0]),
    ]
    for truth, prediction in cases:
        metrics = point_metrics(truth, prediction)
        assert math.isclose(metrics["mcc"], matthews_corrcoef(truth, prediction))
        assert math.isclose(metrics["auroc_binary_score"], roc_auc_score(truth, prediction))
        assert math.isclose(metrics["auprc_binary_score"], average_precision_score(truth, prediction))
        assert math.isclose(
            metrics["brier_binary_score"],
            sum(t != p for t, p in zip(truth, prediction)) / len(truth),
        )


def test_exact_interval_has_boundary_behavior():
    assert exact_interval(0, 10)[0] == 0
    assert exact_interval(10, 10)[1] == 1


def test_raw_predictions_respect_final_sample_qc(tmp_path):
    panel = tmp_path / "results/panels/zenodo_read_marker_predictions.csv"
    panel.parent.mkdir(parents=True)
    fields = ["dataset_id", "isolate_id", "antibiotic", "prediction", "sequence_support", "marker_summary", "callable"]
    with panel.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([
            {"dataset_id": "ZENODO_10369064", "isolate_id": "PASS1", "antibiotic": "clarithromycin", "prediction": "R", "sequence_support": "RAW_READ_SUPPORTED", "marker_summary": "A2143:1.0", "callable": "yes"},
            {"dataset_id": "ZENODO_10369064", "isolate_id": "FAIL1", "antibiotic": "clarithromycin", "prediction": "R", "sequence_support": "RAW_READ_SUPPORTED", "marker_summary": "A2143:1.0", "callable": "yes"},
        ])
    qc = tmp_path / "results/qc/assembly_qc_with_checkm2.csv"
    qc.parent.mkdir(parents=True)
    with qc.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["dataset_id", "isolate_id", "final_qc_status"])
        writer.writeheader()
        writer.writerows([
            {"dataset_id": "ZENODO_10369064", "isolate_id": "PASS1", "final_qc_status": "PASS"},
            {"dataset_id": "ZENODO_10369064", "isolate_id": "FAIL1", "final_qc_status": "FAIL"},
        ])
    rows = {row["isolate_id"]: row for row in raw_read_predictions(tmp_path)}
    assert rows["PASS1"]["prediction"] == "R"
    assert rows["PASS1"]["basic_qc_status"] == "PASS"
    assert rows["FAIL1"]["prediction"] == "UNCALLABLE"
    assert rows["FAIL1"]["callable"] == "no"
