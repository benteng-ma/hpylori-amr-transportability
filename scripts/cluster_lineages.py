#!/usr/bin/env python3
"""Phenotype-blind core-SNP clustering for leave-lineage-out validation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from Bio import SeqIO
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD


def allele_feature_matrix(
    alignment: Path, training_dataset: str | None = None
) -> tuple[list[str], sparse.csr_matrix, list[tuple[int, int]]]:
    records = list(SeqIO.parse(alignment, "fasta"))
    if len(records) < 2:
        raise ValueError("core alignment has fewer than two records")
    reference_records = [record for record in records if record.id == "Reference"]
    if len(reference_records) != 1:
        raise ValueError(f"expected exactly one Reference record, found {len(reference_records)}")
    reference_record = reference_records[0]
    reference = np.frombuffer(str(reference_record.seq).upper().encode("ascii"), dtype=np.uint8)
    names: list[str] = []
    sample_changes: list[list[tuple[int, int]]] = []
    canonical = np.frombuffer(b"ACGT", dtype=np.uint8)
    for record in records:
        if record is reference_record:
            continue
        sequence = np.frombuffer(str(record.seq).upper().encode("ascii"), dtype=np.uint8)
        if len(sequence) != len(reference):
            raise ValueError("unaligned FASTA lengths")
        valid = np.isin(sequence, canonical) & np.isin(reference, canonical)
        indices = np.flatnonzero(valid & (sequence != reference))
        sample_changes.append([(int(position), int(sequence[position])) for position in indices])
        names.append(record.id)
    if training_dataset is None:
        training_rows = np.ones(len(names), dtype=bool)
    else:
        training_rows = np.asarray([name.startswith(f"{training_dataset}__") for name in names], dtype=bool)
        if not training_rows.any():
            raise ValueError(f"no alignment records for training dataset {training_dataset}")
    allele_counts: dict[tuple[int, int], int] = {}
    for include, changes in zip(training_rows, sample_changes):
        if include:
            for key in changes:
                allele_counts[key] = allele_counts.get(key, 0) + 1
    n_training = int(training_rows.sum())
    allele_features = sorted(key for key, count in allele_counts.items() if 0 < count < n_training)
    feature_index = {key: index for index, key in enumerate(allele_features)}
    row_indices: list[int] = []
    column_indices: list[int] = []
    for row_index, changes in enumerate(sample_changes):
        for key in changes:
            if key in feature_index:
                row_indices.append(row_index)
                column_indices.append(feature_index[key])
    data = np.ones(len(row_indices), dtype=np.float32)
    matrix = sparse.csr_matrix((data, (row_indices, column_indices)), shape=(len(names), len(allele_features)))
    return names, matrix, allele_features


def variant_matrix(alignment: Path, training_dataset: str | None = None) -> tuple[list[str], sparse.csr_matrix, int]:
    names, matrix, allele_features = allele_feature_matrix(alignment, training_dataset)
    variable_sites = len({position for position, _base in allele_features})
    return names, matrix, variable_sites


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--alignment", type=Path, default=Path("data/phylogeny/core/core.aln"))
    parser.add_argument("--clusters", type=int, default=8)
    parser.add_argument("--components", type=int, default=20)
    parser.add_argument("--training-dataset", default="HPGP_GLOBAL")
    parser.add_argument("--seed", type=int, default=20260830)
    parser.add_argument("--output", type=Path, default=Path("results/lineage_validation/lineage_assignments.csv"))
    args = parser.parse_args()
    root = args.root.resolve()
    alignment = args.alignment if args.alignment.is_absolute() else root / args.alignment
    names, matrix, n_variants = variant_matrix(alignment, args.training_dataset)
    training_rows = np.asarray([name.startswith(f"{args.training_dataset}__") for name in names], dtype=bool)
    n_components = min(args.components, int(training_rows.sum()) - 1, matrix.shape[1] - 1)
    reducer = TruncatedSVD(n_components=n_components, random_state=args.seed)
    reducer.fit(matrix[training_rows])
    embedding = reducer.transform(matrix)
    clusterer = KMeans(n_clusters=args.clusters, random_state=args.seed, n_init=50)
    clusterer.fit(embedding[training_rows])
    raw_labels = clusterer.predict(embedding)
    ordering = {
        label: index + 1
        for index, label in enumerate(sorted(set(clusterer.labels_), key=lambda value: float(embedding[training_rows][clusterer.labels_ == value, 0].mean())))
    }
    rows = []
    for name, label, coordinates in zip(names, raw_labels, embedding):
        if "__" not in name:
            raise ValueError(f"sample name lacks dataset delimiter: {name}")
        dataset, isolate = name.split("__", 1)
        row = {
            "dataset_id": dataset, "isolate_id": isolate,
            "lineage_recomputed": f"SNP_CLUSTER_{ordering[int(label)]:02d}",
            "clustering_method": f"core_SNP_TruncatedSVD20_KMeans8_fit_{args.training_dataset}",
            "core_variable_sites": n_variants,
        }
        for index, value in enumerate(coordinates[:10], start=1):
            row[f"PC{index}"] = value
        rows.append(row)
    rows.sort(key=lambda row: (str(row["lineage_recomputed"]), str(row["dataset_id"]), str(row["isolate_id"])))
    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["dataset_id", "isolate_id", "lineage_recomputed", "clustering_method", "core_variable_sites", *[f"PC{i}" for i in range(1, 11)]]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
