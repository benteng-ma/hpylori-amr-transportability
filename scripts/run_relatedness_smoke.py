#!/usr/bin/env python3
"""Run pairwise Mash/skani checks on a small Phase 0 assembly panel."""

from __future__ import annotations

import argparse
import csv
import itertools
import subprocess
from pathlib import Path


def cohort(path: Path) -> str:
    name = path.name.upper()
    if name.startswith("HPGP_"):
        return "HPGP_GLOBAL"
    if name.startswith("CHINA_"):
        return "CHINA_NINGXIA_2022"
    if name.startswith("HP"):
        return "LINQU_2025"
    return "UNKNOWN"


def run(command: list[str]) -> str:
    return subprocess.run(command, check=True, capture_output=True, text=True).stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mash", default="mash")
    parser.add_argument("--skani", default="skani")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("assemblies", nargs="+", type=Path)
    args = parser.parse_args()

    rows = []
    for first, second in itertools.combinations(args.assemblies, 2):
        mash_fields = run([args.mash, "dist", str(first), str(second)]).split("\t")
        skani_lines = run([args.skani, "dist", "-q", str(first), "-r", str(second)]).splitlines()
        skani_fields = skani_lines[-1].split("\t")
        ani = float(skani_fields[2])
        align_ref = float(skani_fields[3]) / 100
        align_query = float(skani_fields[4]) / 100
        mash_distance = float(mash_fields[2])
        near_clone = ani >= 99.9 and min(align_ref, align_query) >= 0.90 and mash_distance <= 0.001
        rows.append({
            "assembly_a": first.name,
            "cohort_a": cohort(first),
            "assembly_b": second.name,
            "cohort_b": cohort(second),
            "mash_distance": mash_distance,
            "mash_p_value": mash_fields[3],
            "mash_shared_hashes": mash_fields[4],
            "skani_ani_percent": ani,
            "skani_align_fraction_a": align_query,
            "skani_align_fraction_b": align_ref,
            "near_clone_phase0": "yes" if near_clone else "no",
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
