#!/usr/bin/env python3
"""Small assembly QC utility for Phase 0 smoke tests."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def fasta_lengths(path: Path) -> list[int]:
    lengths: list[int] = []
    current = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    return lengths


def n50(lengths: list[int]) -> int:
    half = sum(lengths) / 2
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= half:
            return length
    return 0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("fastas", nargs="+", type=Path)
    args = parser.parse_args()
    rows = []
    for path in args.fastas:
        lengths = fasta_lengths(path)
        rows.append({
            "path": path.as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "contigs": len(lengths),
            "assembly_size_bp": sum(lengths),
            "n50_bp": n50(lengths),
            "longest_contig_bp": max(lengths, default=0),
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()

