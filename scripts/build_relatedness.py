#!/usr/bin/env python3
"""Full-cohort Mash/skani relatedness and frozen near-clone grouping."""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def run(command: list[str], output: Path | None = None) -> str:
    if output is None:
        return subprocess.run(command, check=True, capture_output=True, text=True).stdout
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        subprocess.run(command, check=True, stdout=handle, text=True)
    return ""


def sample(path_text: str, root: Path) -> tuple[str, str, str]:
    path = Path(path_text)
    try:
        # Mash echoes the absolute paths written by this script. Avoid a
        # filesystem-resolving stat call for every member of every pair;
        # 517 genomes otherwise trigger more than half a million slow calls
        # through the Windows/WSL mount for the same lexical operation.
        relative = path.relative_to(root)
    except ValueError:
        try:
            relative = path.resolve().relative_to(root)
        except ValueError:
            relative = path
    parts = relative.parts
    dataset = parts[-2]
    isolate = Path(parts[-1]).stem
    return dataset, isolate, relative.as_posix()


def skani_one(pair: dict[str, object], skani: str) -> dict[str, object]:
    completed = subprocess.run(
        [skani, "dist", "-q", str(pair["path_a"]), "-r", str(pair["path_b"])],
        check=True, capture_output=True, text=True,
    )
    fields = completed.stdout.strip().splitlines()[-1].split("\t")
    pair["skani_ani_percent"] = float(fields[2])
    pair["skani_align_fraction_a"] = float(fields[4]) / 100
    pair["skani_align_fraction_b"] = float(fields[3]) / 100
    pair["near_clone"] = "yes" if (
        float(pair["skani_ani_percent"]) >= 99.9
        and float(pair["skani_align_fraction_a"]) >= 0.90
        and float(pair["skani_align_fraction_b"]) >= 0.90
        and float(pair["mash_distance"]) <= 0.001
    ) else "no"
    return pair


class UnionFind:
    def __init__(self, values: list[tuple[str, str]]):
        self.parent = {value: value for value in values}

    def find(self, value: tuple[str, str]) -> tuple[str, str]:
        if self.parent[value] != value:
            self.parent[value] = self.find(self.parent[value])
        return self.parent[value]

    def union(self, left: tuple[str, str], right: tuple[str, str]) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--mash", default="mash")
    parser.add_argument("--skani", default="skani")
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--candidate-distance", type=float, default=0.005)
    parser.add_argument("--pair-output", type=Path, default=Path("results/qc/pairwise_relatedness_candidates.csv"))
    parser.add_argument("--group-output", type=Path, default=Path("results/qc/near_clone_groups.csv"))
    args = parser.parse_args()
    root = args.root.resolve()
    final_qc = root / "results/qc/assembly_qc_with_checkm2.csv"
    basic_qc = root / "results/qc/assembly_qc.csv"
    qc_path = final_qc if final_qc.exists() else basic_qc
    assemblies = sorted(
        root / row["assembly_path"] for row in read_csv(qc_path)
        if row.get("final_qc_status", row["basic_qc_status"]) == "PASS"
    )
    if len(assemblies) < 2:
        raise SystemExit("fewer than two assemblies")
    working = root / "data/phylogeny/mash"
    working.mkdir(parents=True, exist_ok=True)
    list_path = working / "assemblies.txt"
    list_path.write_text("".join(f"{path.resolve()}\n" for path in assemblies), encoding="utf-8")
    prefix = working / "all_assemblies"
    run([args.mash, "sketch", "-p", str(args.threads), "-l", str(list_path), "-o", str(prefix)])
    distances = working / "all_pairwise_mash.tsv"
    run([args.mash, "dist", "-p", str(args.threads), str(prefix) + ".msh", str(prefix) + ".msh"], distances)

    candidates: list[dict[str, object]] = []
    seen: set[tuple[str, str, str, str]] = set()
    with distances.open(encoding="utf-8") as handle:
        for line in handle:
            left, right, distance, p_value, shared = line.rstrip("\n").split("\t")
            dataset_a, isolate_a, relative_a = sample(left, root)
            dataset_b, isolate_b, relative_b = sample(right, root)
            if (dataset_a, isolate_a) == (dataset_b, isolate_b):
                continue
            ordered = sorted(((dataset_a, isolate_a, relative_a), (dataset_b, isolate_b, relative_b)))
            key = (ordered[0][0], ordered[0][1], ordered[1][0], ordered[1][1])
            if key in seen:
                continue
            seen.add(key)
            value = float(distance)
            if value > args.candidate_distance:
                continue
            candidates.append({
                "dataset_a": ordered[0][0], "isolate_a": ordered[0][1], "path_a": str(root / ordered[0][2]),
                "dataset_b": ordered[1][0], "isolate_b": ordered[1][1], "path_b": str(root / ordered[1][2]),
                "mash_distance": value, "mash_p_value": p_value, "mash_shared_hashes": shared,
            })
    completed: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(skani_one, pair, args.skani): pair for pair in candidates}
        for index, future in enumerate(as_completed(futures), start=1):
            completed.append(future.result())
            print(f"[{index}/{len(candidates)}] skani candidate", flush=True)
    completed.sort(key=lambda row: (str(row["dataset_a"]), str(row["isolate_a"]), str(row["dataset_b"]), str(row["isolate_b"])))
    pair_output = args.pair_output if args.pair_output.is_absolute() else root / args.pair_output
    pair_output.parent.mkdir(parents=True, exist_ok=True)
    pair_fields = [
        "dataset_a", "isolate_a", "dataset_b", "isolate_b", "mash_distance", "mash_p_value", "mash_shared_hashes",
        "skani_ani_percent", "skani_align_fraction_a", "skani_align_fraction_b", "near_clone",
    ]
    with pair_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=pair_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(completed)

    values = [(path.parent.name, path.stem) for path in assemblies]
    union_find = UnionFind(values)
    for row in completed:
        if row["near_clone"] == "yes":
            union_find.union((str(row["dataset_a"]), str(row["isolate_a"])), (str(row["dataset_b"]), str(row["isolate_b"])))
    members: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    for value in values:
        members[union_find.find(value)].append(value)
    ordered_groups = sorted(members.values(), key=lambda group: min(group))
    group_rows = []
    for index, group in enumerate(ordered_groups, start=1):
        group_id = f"NCG{index:04d}"
        for dataset, isolate in sorted(group):
            group_rows.append({"dataset_id": dataset, "isolate_id": isolate, "near_clone_group": group_id, "group_size": len(group)})
    group_output = args.group_output if args.group_output.is_absolute() else root / args.group_output
    group_output.parent.mkdir(parents=True, exist_ok=True)
    with group_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["dataset_id", "isolate_id", "near_clone_group", "group_size"])
        writer.writeheader()
        writer.writerows(group_rows)


if __name__ == "__main__":
    main()
