from scripts.run_23s_callability_sensitivity import qualifies, recovery_class


def audit_hit(task: str, identity: float, coverage: float, spanning: bool) -> dict[str, object]:
    return {
        "blast_task": task,
        "identity": identity,
        "coverage": coverage,
        "marker_spanning": spanning,
    }


def test_qualifies_requires_marker_span() -> None:
    assert not qualifies(audit_hit("megablast", 1.0, 1.0, False), 0.90, 0.05)


def test_recovery_class_distinguishes_relaxed_threshold() -> None:
    rows = [audit_hit("blastn", 0.85, 0.10, True)]
    assert recovery_class(rows) == "RESCUED_ONLY_UNDER_RELAXED_THRESHOLD"


def test_recovery_class_distinguishes_partial_hit() -> None:
    rows = [audit_hit("megablast", 0.95, 0.50, False)]
    assert recovery_class(rows) == "PARTIAL_23S_NO_MARKER_SPAN"
