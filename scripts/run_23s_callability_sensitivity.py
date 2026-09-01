#!/usr/bin/env python3
"""Post-freeze sensitivity audit of assembly-based 23S target recovery.

The frozen primary caller is not changed by this analysis. The grid varies the
BLAST task, identity threshold and aligned-query coverage while retaining the
biologically essential requirement that both resistance-marker bases are
present in the alignment.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import SeqIO

try:
    from call_known_variants import blast, extract_reference_features, query_to_subject_map
except ModuleNotFoundError:  # imported as scripts.run_23s_callability_sensitivity in tests
    from scripts.call_known_variants import blast, extract_reference_features, query_to_subject_map


TASKS = ("megablast", "blastn")
IDENTITIES = (0.70, 0.80, 0.90, 0.95, 0.98)
COVERAGES = (0.00, 0.01, 0.05, 0.25, 0.50, 0.90)
MARKER_POSITIONS = (2143, 2144)
FROZEN_TASK = "megablast"
FROZEN_IDENTITY = 0.90
FROZEN_COVERAGE = 0.05


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def hit_summary(hit: dict[str, str], query_length: int) -> dict[str, object]:
    mapping = query_to_subject_map(hit)
    covered = sorted(mapping)
    marker_spanning = all(position in mapping for position in MARKER_POSITIONS)
    if marker_spanning:
        nearest_distance = 0
    elif covered:
        nearest_distance = min(abs(position - marker) for position in covered for marker in MARKER_POSITIONS)
    else:
        nearest_distance = query_length
    return {
        "subject_sequence": hit["sseqid"],
        "query_start": int(hit["qstart"]),
        "query_end": int(hit["qend"]),
        "subject_start": int(hit["sstart"]),
        "subject_end": int(hit["send"]),
        "identity": float(hit["pident"]) / 100,
        "coverage": int(hit["length"]) / query_length,
        "marker_spanning": marker_spanning,
        "nearest_marker_distance_bp": nearest_distance,
    }


def qualifies(row: dict[str, object], identity: float, coverage: float) -> bool:
    return bool(row["marker_spanning"]) and float(row["identity"]) >= identity and float(row["coverage"]) >= coverage


def recovery_class(rows: list[dict[str, object]]) -> str:
    frozen = any(qualifies(row, FROZEN_IDENTITY, FROZEN_COVERAGE) for row in rows if row["blast_task"] == FROZEN_TASK)
    if frozen:
        return "CALLABLE_MARKER_SPANNING"
    relaxed = any(
        qualifies(row, min(IDENTITIES), min(COVERAGES))
        for row in rows
    )
    if relaxed:
        return "RESCUED_ONLY_UNDER_RELAXED_THRESHOLD"
    plausible = any(float(row["identity"]) >= min(IDENTITIES) and float(row["coverage"]) >= 0.01 for row in rows)
    if plausible:
        return "PARTIAL_23S_NO_MARKER_SPAN"
    return "NO_23S_LIKE_HIT"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--dataset", default="CHINA_NINGXIA_2022")
    args = parser.parse_args()
    root = args.root.resolve()
    rrna, _ = extract_reference_features(
        root / "data/raw/reference/U27270.1.gb",
        root / "data/interim/reference/ncbi_dataset/data/GCF_000008525.1/genomic.gbff",
        root / "data/interim/reference/markers",
    )
    query_length = len(SeqIO.read(rrna, "fasta").seq)
    assemblies = sorted((root / "data/assemblies" / args.dataset).glob("*.fna"))
    if not assemblies:
        raise SystemExit(f"no assemblies found for {args.dataset}")

    qc = {
        (row["dataset_id"], row["isolate_id"]): row
        for row in read_csv(root / "results/qc/assembly_qc_with_checkm2.csv")
    }
    primary_calls = read_csv(root / "results/panels/assembly_marker_calls.csv")
    primary_callable = {
        (row["dataset_id"], row["isolate_id"])
        for row in primary_calls
        if row["gene"] == "23S_rRNA" and row["status"] == "PASS"
    }

    all_hits: list[dict[str, object]] = []
    by_sample: dict[str, list[dict[str, object]]] = {path.stem: [] for path in assemblies}
    for index, assembly in enumerate(assemblies, start=1):
        for task in TASKS:
            raw_hits = blast(rrna, assembly, args.blastn, task=task)
            for rank, raw_hit in enumerate(raw_hits, start=1):
                row = {
                    "dataset_id": args.dataset,
                    "isolate_id": assembly.stem,
                    "blast_task": task,
                    "hit_rank": rank,
                    **hit_summary(raw_hit, query_length),
                }
                all_hits.append(row)
                by_sample[assembly.stem].append(row)
        print(f"[{index}/{len(assemblies)}] {assembly.stem}", flush=True)

    sample_rows: list[dict[str, object]] = []
    for isolate_id, hits in by_sample.items():
        key = (args.dataset, isolate_id)
        final_qc = qc.get(key, {}).get("final_qc_status", "MISSING")
        for task in TASKS:
            task_hits = [row for row in hits if row["blast_task"] == task]
            spanning_hits = [row for row in task_hits if row["marker_spanning"]]
            sample_rows.append({
                "dataset_id": args.dataset,
                "isolate_id": isolate_id,
                "final_qc_status": final_qc,
                "blast_task": task,
                "n_hits": len(task_hits),
                "n_marker_spanning_hits": len(spanning_hits),
                "max_identity": max((float(row["identity"]) for row in task_hits), default=""),
                "max_coverage": max((float(row["coverage"]) for row in task_hits), default=""),
                "max_spanning_identity": max((float(row["identity"]) for row in spanning_hits), default=""),
                "max_spanning_coverage": max((float(row["coverage"]) for row in spanning_hits), default=""),
                "nearest_marker_distance_bp": min((int(row["nearest_marker_distance_bp"]) for row in task_hits), default=""),
                "frozen_callable": "yes" if task == FROZEN_TASK and any(qualifies(row, FROZEN_IDENTITY, FROZEN_COVERAGE) for row in task_hits) else "no",
                "primary_caller_callable": "yes" if key in primary_callable else "no",
                "recovery_class": recovery_class(hits),
            })

    grid_rows: list[dict[str, object]] = []
    for task in TASKS:
        for identity in IDENTITIES:
            for coverage in COVERAGES:
                callable_all = 0
                callable_final_qc = 0
                total_final_qc = 0
                for isolate_id, hits in by_sample.items():
                    key = (args.dataset, isolate_id)
                    task_hits = [row for row in hits if row["blast_task"] == task]
                    is_callable = any(qualifies(row, identity, coverage) for row in task_hits)
                    callable_all += int(is_callable)
                    if qc.get(key, {}).get("final_qc_status") == "PASS":
                        total_final_qc += 1
                        callable_final_qc += int(is_callable)
                grid_rows.append({
                    "blast_task": task,
                    "minimum_identity": identity,
                    "minimum_query_coverage": coverage,
                    "require_both_marker_bases": "yes",
                    "n_assemblies": len(assemblies),
                    "n_callable_all": callable_all,
                    "n_final_qc": total_final_qc,
                    "n_callable_final_qc": callable_final_qc,
                    "callable_fraction_final_qc": callable_final_qc / total_final_qc if total_final_qc else "",
                    "is_frozen_configuration": "yes" if (task, identity, coverage) == (FROZEN_TASK, FROZEN_IDENTITY, FROZEN_COVERAGE) else "no",
                })

    frozen_samples = [row for row in sample_rows if row["blast_task"] == FROZEN_TASK]
    mismatches = [
        row for row in frozen_samples
        if row["frozen_callable"] != row["primary_caller_callable"]
    ]
    if mismatches:
        raise AssertionError(f"sensitivity audit disagrees with primary caller for {len(mismatches)} samples")
    frozen_final_qc = sum(
        row["final_qc_status"] == "PASS" and row["frozen_callable"] == "yes"
        for row in frozen_samples
    )
    if frozen_final_qc != 6:
        raise AssertionError(f"expected six frozen-callable final-QC Ningxia assemblies, observed {frozen_final_qc}")

    class_counts: dict[str, int] = {}
    for row in frozen_samples:
        label = str(row["recovery_class"])
        class_counts[label] = class_counts.get(label, 0) + 1
    summary_rows: list[dict[str, object]] = [
        {"section": "configuration", "metric": "frozen_blast_task", "value": FROZEN_TASK, "denominator": "", "notes": "primary caller"},
        {"section": "configuration", "metric": "frozen_minimum_identity", "value": FROZEN_IDENTITY, "denominator": "", "notes": "fraction"},
        {"section": "configuration", "metric": "frozen_minimum_query_coverage", "value": FROZEN_COVERAGE, "denominator": "", "notes": "fraction"},
        {"section": "result", "metric": "frozen_callable_final_qc", "value": frozen_final_qc, "denominator": sum(row["final_qc_status"] == "PASS" for row in frozen_samples), "notes": "unchanged primary callability"},
    ]
    summary_rows.extend(
        {"section": "recovery_class", "metric": label, "value": count, "denominator": len(frozen_samples), "notes": "all Ningxia assemblies"}
        for label, count in sorted(class_counts.items())
    )

    output_dir = root / "results/extended_analysis"
    write_csv(output_dir / "23s_callability_sensitivity_hits.csv", all_hits)
    write_csv(output_dir / "23s_callability_sensitivity_samples.csv", sample_rows)
    write_csv(output_dir / "23s_callability_sensitivity_grid.csv", grid_rows)
    write_csv(output_dir / "23s_callability_sensitivity_summary.csv", summary_rows)


if __name__ == "__main__":
    main()
