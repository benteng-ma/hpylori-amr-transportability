#!/usr/bin/env python3
"""Build an accession-to-analysis crosswalk for every isolate in the benchmark."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def prefixed(source: dict[str, str], prefix: str, fields: list[str]) -> dict[str, str]:
    return {f"{prefix}_{field}": source.get(field, "") for field in fields}


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    isolates = read_csv(root / "metadata/isolate_manifest.csv")
    assemblies = {
        (row["dataset_id"], row["isolate_id"]): row
        for row in read_csv(root / "metadata/phase2/acquisition_assemblies.csv")
    }
    reads = {
        (row["dataset_id"], row["isolate_id"]): row
        for row in read_csv(root / "metadata/phase2/acquisition_reads.csv")
    }
    qc = {
        (row["dataset_id"], row["isolate_id"]): row
        for row in read_csv(root / "results/qc/assembly_qc_with_checkm2.csv")
    }
    lineages = {
        (row["dataset_id"], row["isolate_id"]): row
        for row in read_csv(root / "results/lineage_validation/lineage_assignments.csv")
    }
    predictions = {
        (row["dataset_id"], row["isolate_id"], row["antibiotic"]): row
        for row in read_csv(root / "results/external_validation/sample_level_predictions.csv")
    }

    acquisition_fields = [
        "sequence_type", "accession", "source_url", "local_path", "bytes", "sha256", "md5",
        "validation_status", "acquisition_status", "accessed_at",
    ]
    phenotype_fields = [
        "phenotype", "phenotype_quality", "phenotype_recomputed", "mic_raw", "mic_numeric",
        "mic_operator", "borderline_mic", "ast_method", "medium", "breakpoint_standard",
        "breakpoint_version", "prediction", "sequence_support", "marker_summary", "callable",
        "analysis_status", "correct",
    ]
    rows: list[dict[str, object]] = []
    for isolate in isolates:
        key = (isolate["dataset_id"], isolate["isolate_id"])
        if key not in qc:
            continue
        assembly = assemblies.get(key, {})
        read = reads.get(key, {})
        quality = qc.get(key, {})
        lineage = lineages.get(key, {})
        clarithromycin = predictions.get((*key, "clarithromycin"), {})
        levofloxacin = predictions.get((*key, "levofloxacin"), {})
        rows.append({
            "dataset_id": key[0],
            "isolate_id": key[1],
            "patient_id": isolate.get("patient_id", ""),
            "sample_id": isolate.get("sample_id", ""),
            "assembly_id_manifest": isolate.get("assembly_id", ""),
            "country": isolate.get("country", ""),
            "site": isolate.get("site", ""),
            "collection_year": isolate.get("collection_year", ""),
            "clinical_diagnosis": isolate.get("clinical_diagnosis", ""),
            "primary_or_post_treatment": isolate.get("primary_or_post_treatment", ""),
            "included": isolate.get("included", ""),
            "exclusion_reason": isolate.get("exclusion_reason", ""),
            **prefixed(assembly, "assembly_acquisition", acquisition_fields),
            **prefixed(read, "read_acquisition", acquisition_fields),
            "analysis_assembly_path": quality.get("assembly_path", ""),
            "contigs": quality.get("contigs", ""),
            "assembly_size_bp": quality.get("assembly_size_bp", ""),
            "n50_bp": quality.get("n50_bp", ""),
            "gc_percent": quality.get("gc_percent", ""),
            "completeness_percent": quality.get("completeness_percent", ""),
            "contamination_percent": quality.get("contamination_percent", ""),
            "species_confirmation": quality.get("species_confirmation", ""),
            "final_qc_status": quality.get("final_qc_status", ""),
            "lineage_recomputed": lineage.get("lineage_recomputed", isolate.get("lineage_recomputed", "")),
            "lineage_method": lineage.get("clustering_method", ""),
            "PC1": lineage.get("PC1", ""),
            "PC2": lineage.get("PC2", ""),
            "near_clone_group": clarithromycin.get("near_clone_group", levofloxacin.get("near_clone_group", isolate.get("near_clone_group", ""))),
            **prefixed(clarithromycin, "clarithromycin", phenotype_fields),
            **prefixed(levofloxacin, "levofloxacin", phenotype_fields),
        })

    if len(rows) != 526 or len({(row["dataset_id"], row["isolate_id"]) for row in rows}) != 526:
        raise AssertionError("crosswalk must contain exactly one row for each of 526 isolates")
    write_csv(root / "results/external_validation/isolate_sequence_phenotype_crosswalk.csv", rows)

    summary: list[dict[str, object]] = []
    for dataset_id in sorted({str(row["dataset_id"]) for row in rows}):
        group = [row for row in rows if row["dataset_id"] == dataset_id]
        sequence_types = Counter(
            "raw_reads" if row["read_acquisition_accession"] else "public_assembly"
            for row in group
        )
        summary.append({
            "dataset_id": dataset_id,
            "n_isolates": len(group),
            "n_public_assemblies": sequence_types["public_assembly"],
            "n_raw_read_packages": sequence_types["raw_reads"],
            "n_final_qc_pass": sum(row["final_qc_status"] == "PASS" for row in group),
            "n_clarithromycin_phenotypes": sum(row["clarithromycin_phenotype"] in {"R", "S"} for row in group),
            "n_clarithromycin_callable": sum(row["clarithromycin_callable"] == "yes" for row in group),
            "n_levofloxacin_phenotypes": sum(row["levofloxacin_phenotype"] in {"R", "S"} for row in group),
            "n_levofloxacin_callable": sum(row["levofloxacin_callable"] == "yes" for row in group),
            "n_recomputed_lineages": len({row["lineage_recomputed"] for row in group if row["lineage_recomputed"]}),
        })
    write_csv(root / "results/external_validation/isolate_sequence_phenotype_crosswalk_summary.csv", summary)


if __name__ == "__main__":
    main()
