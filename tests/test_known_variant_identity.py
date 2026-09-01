from scripts.call_known_variants import best_nonredundant_hits


def hit(identity: float, length: int = 100) -> dict[str, str]:
    return {
        "qseqid": "q", "sseqid": "s", "pident": str(identity), "length": str(length),
        "qstart": "1", "qend": str(length), "sstart": "1", "send": str(length),
        "qseq": "A" * length, "sseq": "A" * length,
    }


def test_identity_threshold_is_applied_before_nonredundancy() -> None:
    retained = best_nonredundant_hits([hit(89.99), hit(90.0)], 0.50, 100, 0.90)
    assert [float(row["pident"]) for row in retained] == [90.0]


def test_identity_threshold_can_be_relaxed_for_post_freeze_sensitivity() -> None:
    retained = best_nonredundant_hits([hit(80.0)], 0.50, 100, 0.80)
    assert len(retained) == 1
