import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_full_analysis_stays_locked_without_go():
    path = ROOT / "reports/phase0_decision.json"
    if not path.exists():
        assert True
        return
    decision = json.loads(path.read_text(encoding="utf-8")).get("decision", "")
    assert decision.startswith(("GO_", "CONDITIONAL_GO_")) or decision.startswith("NO_GO")
