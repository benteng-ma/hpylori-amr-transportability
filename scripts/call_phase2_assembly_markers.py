#!/usr/bin/env python3
"""Call the frozen 23S/gyrA panel across all materialized assemblies."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import yaml

from call_known_variants import call_one, extract_reference_features


def worker(
    assembly: Path,
    root: Path,
    rrna: Path,
    gyra: Path,
    blastn: str,
    minimum_identity: float,
) -> list[dict[str, object]]:
    dataset_id = assembly.parent.name
    isolate_id = assembly.stem
    rows = call_one(assembly, rrna, gyra, blastn, minimum_identity)
    for row in rows:
        row["dataset_id"] = dataset_id
        row["isolate_id"] = isolate_id
        row["assembly_path"] = assembly.relative_to(root).as_posix()
        row["sequence_support"] = "ASSEMBLY_ONLY"
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--output", type=Path, default=Path("results/panels/assembly_marker_calls.csv"))
    args = parser.parse_args()
    root = args.root.resolve()
    with (root / "config/qc.yaml").open(encoding="utf-8") as handle:
        qc_config = yaml.safe_load(handle)
    target_config = qc_config["target_gene"]
    minimum_identity = float(target_config["minimum_identity"])
    if not 0 <= minimum_identity <= 1:
        raise SystemExit("target_gene.minimum_identity must be a fraction from 0 to 1")
    if not bool(target_config["require_target_spanning_alignment"]):
        raise SystemExit("the frozen panel requires target-spanning alignments")
    rrna, gyra = extract_reference_features(
        root / "data/raw/reference/U27270.1.gb",
        root / "data/interim/reference/ncbi_dataset/data/GCF_000008525.1/genomic.gbff",
        root / "data/interim/reference/markers",
    )
    assemblies = sorted((root / "data/assemblies").glob("*/*.fna"))
    if not assemblies:
        raise SystemExit("no materialized assemblies")
    rows: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(worker, path, root, rrna, gyra, args.blastn, minimum_identity): path
            for path in assemblies
        }
        for index, future in enumerate(as_completed(futures), start=1):
            path = futures[future]
            calls = future.result()
            rows.extend(calls)
            print(f"[{index}/{len(assemblies)}] PASS {path.parent.name} {path.stem}", flush=True)
    rows.sort(key=lambda row: (str(row["dataset_id"]), str(row["isolate_id"]), str(row["gene"]), str(row["copy"]), int(row["position"])))
    fieldnames = [
        "dataset_id", "isolate_id", "assembly_path", "sequence_support", "assembly", "gene", "copy",
        "coordinate_system", "position", "subject_sequence", "subject_position", "reference", "observed",
        "change", "known_resistance_marker", "hit_identity_percent", "hit_coverage", "status",
    ]
    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
