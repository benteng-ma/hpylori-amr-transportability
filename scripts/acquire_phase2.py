#!/usr/bin/env python3
"""Acquire the frozen Phase 2 public sequence set with resumable checks.

Every task is constructed from the committed isolate and repository manifests;
no downstream performance result is used to discover or select samples.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import requests


@dataclass(frozen=True)
class Task:
    dataset_id: str
    isolate_id: str
    sequence_type: str
    accession: str
    source_url: str
    local_path: str
    expected_bytes: int | None = None
    expected_checksum_type: str = ""
    expected_checksum: str = ""


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def digest(path: Path, algorithm: str) -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def validate_fasta_gzip(path: Path) -> None:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        if not handle.readline().startswith(">"):
            raise ValueError(f"not a gzip FASTA: {path}")


def validate_zip(path: Path) -> None:
    import zipfile

    with zipfile.ZipFile(path) as archive:
        bad = archive.testzip()
        if bad is not None:
            raise ValueError(f"corrupt ZIP member {bad}: {path}")
        names = [name for name in archive.namelist() if not name.endswith("/")]
        if len(names) < 2:
            raise ValueError(f"paired-read archive has fewer than two files: {path}")


def ncbi_tasks(root: Path) -> list[Task]:
    isolates = [
        row for row in read_csv(root / "metadata/isolate_manifest.csv")
        if row["included"] == "yes" and row["dataset_id"] in {"HPGP_GLOBAL", "CHINA_NINGXIA_2022"}
    ]
    index = read_csv(root / "metadata/phase0/ncbi_assembly_index.csv")
    by_accession: dict[tuple[str, str], dict[str, str]] = {}
    for assembly_row in index:
        for field in ("assembly_accession", "genbank_accession", "refseq_accession"):
            if assembly_row.get(field):
                by_accession[(assembly_row["dataset_id"], assembly_row[field])] = assembly_row
    tasks: list[Task] = []
    for row in isolates:
        key = (row["dataset_id"], row["assembly_id"])
        assembly = by_accession.get(key)
        if assembly is None:
            raise ValueError(f"no exact assembly index row for {key}")
        base = assembly["ftp_genbank"] or assembly["ftp_refseq"]
        if not base:
            raise ValueError(f"no public FTP path for {key}")
        base = base.replace("ftp://", "https://")
        filename = base.rsplit("/", 1)[-1] + "_genomic.fna.gz"
        tasks.append(Task(
            dataset_id=row["dataset_id"], isolate_id=row["isolate_id"], sequence_type="assembly",
            accession=row["assembly_id"], source_url=f"{base}/{filename}",
            local_path=(Path("data/raw/assemblies") / row["dataset_id"] / f"{row['isolate_id']}.fna.gz").as_posix(),
        ))
    return tasks


def zenodo_tasks(root: Path) -> list[Task]:
    isolates = {
        row["isolate_id"] for row in read_csv(root / "metadata/isolate_manifest.csv")
        if row["included"] == "yes" and row["dataset_id"] == "ZENODO_10369064"
    }
    record = json.loads((root / "metadata/raw/zenodo_10369064_2026-08-30.json").read_text(encoding="utf-8"))
    tasks: list[Task] = []
    for item in record["files"]:
        key = item["key"]
        isolate = Path(key).stem
        if isolate not in isolates or not key.lower().endswith(".zip"):
            continue
        checksum_type, checksum = item.get("checksum", ":").split(":", 1)
        tasks.append(Task(
            dataset_id="ZENODO_10369064", isolate_id=isolate, sequence_type="paired_reads",
            accession="10.5281/zenodo.10369064", source_url=item["links"].get("content") or item["links"]["self"],
            local_path=(Path("data/raw/reads/zenodo_10369064") / key).as_posix(),
            expected_bytes=int(item["size"]), expected_checksum_type=checksum_type, expected_checksum=checksum,
        ))
    missing = isolates - {task.isolate_id for task in tasks}
    if missing:
        raise ValueError(f"Zenodo file missing for isolates: {sorted(missing)}")
    return tasks


def download(task: Task, root: Path, retries: int, timeout: int) -> dict[str, object]:
    destination = root / task.local_path
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")

    def valid_existing() -> bool:
        if not destination.exists():
            return False
        if task.expected_bytes is not None and destination.stat().st_size != task.expected_bytes:
            return False
        if task.expected_checksum:
            return digest(destination, task.expected_checksum_type) == task.expected_checksum.lower()
        return True

    status = "already_present" if valid_existing() else "downloaded"
    if status == "downloaded":
        for attempt in range(1, retries + 1):
            try:
                offset = partial.stat().st_size if partial.exists() else 0
                headers = {"Range": f"bytes={offset}-"} if offset else {}
                with requests.get(task.source_url, stream=True, timeout=(30, timeout), headers=headers) as response:
                    if offset and response.status_code != 206:
                        partial.unlink(missing_ok=True)
                        offset = 0
                    response.raise_for_status()
                    mode = "ab" if offset and response.status_code == 206 else "wb"
                    with partial.open(mode) as handle:
                        for chunk in response.iter_content(chunk_size=4 * 1024 * 1024):
                            if chunk:
                                handle.write(chunk)
                partial.replace(destination)
                if task.expected_bytes is not None and destination.stat().st_size != task.expected_bytes:
                    raise ValueError(f"size mismatch {destination.stat().st_size} != {task.expected_bytes}")
                if task.expected_checksum:
                    observed = digest(destination, task.expected_checksum_type)
                    if observed != task.expected_checksum.lower():
                        raise ValueError(f"{task.expected_checksum_type} mismatch")
                break
            except Exception:
                if attempt == retries:
                    raise
                time.sleep(min(30, 2**attempt))

    if task.sequence_type == "assembly":
        validate_fasta_gzip(destination)
    else:
        validate_zip(destination)
    return {
        **asdict(task), "bytes": destination.stat().st_size,
        "sha256": digest(destination, "sha256"), "md5": digest(destination, "md5"),
        "validation_status": "PASS", "acquisition_status": status,
        "accessed_at": datetime.now(timezone.utc).isoformat(timespec="seconds"), "error": "",
    }


def write_manifest(path: Path, rows: Iterable[dict[str, object]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "dataset_id", "isolate_id", "sequence_type", "accession", "source_url", "local_path",
        "expected_bytes", "expected_checksum_type", "expected_checksum", "bytes", "sha256", "md5",
        "validation_status", "acquisition_status", "accessed_at", "error",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--scope", choices=["assemblies", "reads", "all"], default="all")
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--retries", type=int, default=5)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--output", type=Path, default=Path("metadata/phase2/acquisition_manifest.csv"))
    args = parser.parse_args()
    root = args.root.resolve()
    tasks = []
    if args.scope in {"assemblies", "all"}:
        tasks.extend(ncbi_tasks(root))
    if args.scope in {"reads", "all"}:
        tasks.extend(zenodo_tasks(root))

    rows: list[dict[str, object]] = []
    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(download, task, root, args.retries, args.timeout): task for task in tasks}
        for index, future in enumerate(as_completed(futures), start=1):
            task = futures[future]
            try:
                rows.append(future.result())
                print(f"[{index}/{len(tasks)}] PASS {task.dataset_id} {task.isolate_id}", flush=True)
            except Exception as exc:
                failures.append(f"{task.dataset_id}:{task.isolate_id}: {exc}")
                rows.append({
                    **asdict(task), "bytes": "", "sha256": "", "md5": "", "validation_status": "FAIL",
                    "acquisition_status": "failed", "accessed_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    "error": str(exc),
                })
                print(f"[{index}/{len(tasks)}] FAIL {task.dataset_id} {task.isolate_id}: {exc}", flush=True)
    rows.sort(key=lambda row: (str(row["dataset_id"]), str(row["isolate_id"])))
    output = args.output if args.output.is_absolute() else root / args.output
    write_manifest(output, rows)
    if failures:
        raise SystemExit("\n".join(failures))


if __name__ == "__main__":
    main()
