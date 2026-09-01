#!/usr/bin/env python3
"""Materialize and audit the frozen NCBI assembly set."""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def n50(lengths: list[int]) -> int:
    threshold = sum(lengths) / 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= threshold:
            return length
    return 0


def fasta_metrics(path: Path) -> dict[str, object]:
    lengths: list[int] = []
    current = 0
    gc = 0
    atgc = 0
    ambiguous = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
                continue
            sequence = line.strip().upper()
            current += len(sequence)
            gc += sequence.count("G") + sequence.count("C")
            atgc += sum(sequence.count(base) for base in "ATGC")
            ambiguous += len(sequence) - sum(sequence.count(base) for base in "ATGCN") + sequence.count("N")
    if current:
        lengths.append(current)
    total = sum(lengths)
    return {
        "contigs": len(lengths), "assembly_size_bp": total, "n50_bp": n50(lengths),
        "longest_contig_bp": max(lengths, default=0),
        "gc_percent": round(100 * gc / atgc, 5) if atgc else "",
        "ambiguous_bases": ambiguous,
    }


def materialize(row: dict[str, str], root: Path) -> dict[str, object]:
    source = root / row["local_path"]
    output = root / "data/assemblies" / row["dataset_id"] / f"{row['isolate_id']}.fna"
    output.parent.mkdir(parents=True, exist_ok=True)
    if not output.exists() or output.stat().st_mtime_ns < source.stat().st_mtime_ns:
        temporary = output.with_suffix(".fna.part")
        with gzip.open(source, "rb") as input_handle, temporary.open("wb") as output_handle:
            shutil.copyfileobj(input_handle, output_handle, length=4 * 1024 * 1024)
        temporary.replace(output)
    metrics = fasta_metrics(output)
    size_pass = 1_400_000 <= int(metrics["assembly_size_bp"]) <= 1_900_000
    contig_pass = int(metrics["contigs"]) <= 300
    return {
        "dataset_id": row["dataset_id"], "isolate_id": row["isolate_id"], "accession": row["accession"],
        "assembly_path": output.relative_to(root).as_posix(), **metrics,
        "size_gate": "PASS" if size_pass else "FAIL", "contig_gate": "PASS" if contig_pass else "FAIL",
        "completeness_percent": "NOT_COMPUTED", "contamination_percent": "NOT_COMPUTED",
        "species_confirmation": "NCBI_H_PYLORI_RECORD",
        "basic_qc_status": "PASS" if size_pass and contig_pass else "FAIL",
    }


def audit_existing(path: Path, root: Path) -> dict[str, object]:
    metrics = fasta_metrics(path)
    size_pass = 1_400_000 <= int(metrics["assembly_size_bp"]) <= 1_900_000
    contig_pass = int(metrics["contigs"]) <= 300
    return {
        "dataset_id": path.parent.name, "isolate_id": path.stem, "accession": "10.5281/zenodo.10369064",
        "assembly_path": path.relative_to(root).as_posix(), **metrics,
        "size_gate": "PASS" if size_pass else "FAIL", "contig_gate": "PASS" if contig_pass else "FAIL",
        "completeness_percent": "NOT_COMPUTED", "contamination_percent": "NOT_COMPUTED",
        "species_confirmation": "PENDING_GENOME_QC",
        "basic_qc_status": "PASS" if size_pass and contig_pass else "FAIL",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--acquisition-manifest", type=Path, default=Path("metadata/phase2/acquisition_assemblies.csv"))
    parser.add_argument("--output", type=Path, default=Path("results/qc/assembly_qc.csv"))
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = args.acquisition_manifest if args.acquisition_manifest.is_absolute() else root / args.acquisition_manifest
    rows = [row for row in read_csv(manifest) if row["validation_status"] == "PASS" and row["sequence_type"] == "assembly"]
    results: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(materialize, row, root): row for row in rows}
        for index, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)
            print(f"[{index}/{len(rows)}] {result['basic_qc_status']} {result['dataset_id']} {result['isolate_id']}", flush=True)
    for path in sorted((root / "data/assemblies/ZENODO_10369064").glob("*.fna")):
        result = audit_existing(path, root)
        results.append(result)
        print(f"[ZENODO] {result['basic_qc_status']} {result['isolate_id']}", flush=True)
    results.sort(key=lambda row: (str(row["dataset_id"]), str(row["isolate_id"])))
    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    main()
