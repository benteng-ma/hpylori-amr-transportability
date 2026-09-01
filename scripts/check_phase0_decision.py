#!/usr/bin/env python3
"""Fail closed unless the machine-readable Phase 0 decision authorizes work."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_phase0_decision.py reports/phase0_decision.json")
    path = Path(sys.argv[1])
    if not path.exists():
        raise SystemExit("Phase 0 decision file is missing")
    decision = str(json.loads(path.read_text(encoding="utf-8")).get("decision", ""))
    if not decision.startswith(("GO_", "CONDITIONAL_GO_")):
        raise SystemExit(f"Phase 0 does not authorize analysis: {decision or 'empty decision'}")


if __name__ == "__main__":
    main()
