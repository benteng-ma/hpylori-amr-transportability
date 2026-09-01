#!/usr/bin/env python3
"""Native-curl acquisition path for large Zenodo archives.

This is transport-equivalent to ``acquire_phase2.py`` but avoids the lower
throughput observed for long Python streams on the Windows host.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from acquire_phase2 import Task, digest, validate_zip, zenodo_tasks


def acquire(task: Task, root: Path, curl: str) -> dict[str, object]:
    destination = root / task.local_path
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    valid = destination.exists() and destination.stat().st_size == task.expected_bytes
    if valid and task.expected_checksum:
        valid = digest(destination, task.expected_checksum_type) == task.expected_checksum.lower()
    status = "already_present" if valid else "downloaded"
    if not valid:
        command = [
            curl, "--location", "--fail", "--silent", "--show-error", "--retry", "8", "--retry-all-errors",
            "--connect-timeout", "30", "--continue-at", "-", "--output", str(partial), task.source_url,
        ]
        completed = subprocess.run(command, capture_output=True, text=True)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr.strip() or f"curl exit {completed.returncode}")
        partial.replace(destination)
    if destination.stat().st_size != task.expected_bytes:
        raise ValueError(f"size mismatch {destination.stat().st_size} != {task.expected_bytes}")
    observed = digest(destination, task.expected_checksum_type)
    if observed != task.expected_checksum.lower():
        raise ValueError(f"{task.expected_checksum_type} mismatch")
    validate_zip(destination)
    return {
        **asdict(task), "bytes": destination.stat().st_size, "sha256": digest(destination, "sha256"),
        "md5": digest(destination, "md5"), "validation_status": "PASS", "acquisition_status": status,
        "accessed_at": datetime.now(timezone.utc).isoformat(timespec="seconds"), "error": "",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--curl", default="curl.exe")
    parser.add_argument("--output", type=Path, default=Path("metadata/phase2/acquisition_reads.csv"))
    args = parser.parse_args()
    root = args.root.resolve()
    tasks = zenodo_tasks(root)
    rows: list[dict[str, object]] = []
    failures = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(acquire, task, root, args.curl): task for task in tasks}
        for index, future in enumerate(as_completed(futures), start=1):
            task = futures[future]
            try:
                rows.append(future.result())
                print(f"[{index}/{len(tasks)}] PASS {task.isolate_id}", flush=True)
            except Exception as exc:
                failures.append(f"{task.isolate_id}: {exc}")
                rows.append({
                    **asdict(task), "bytes": "", "sha256": "", "md5": "", "validation_status": "FAIL",
                    "acquisition_status": "failed", "accessed_at": datetime.now(timezone.utc).isoformat(timespec="seconds"), "error": str(exc),
                })
                print(f"[{index}/{len(tasks)}] FAIL {task.isolate_id}: {exc}", flush=True)
    rows.sort(key=lambda row: str(row["isolate_id"]))
    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "dataset_id", "isolate_id", "sequence_type", "accession", "source_url", "local_path",
        "expected_bytes", "expected_checksum_type", "expected_checksum", "bytes", "sha256", "md5",
        "validation_status", "acquisition_status", "accessed_at", "error",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    if failures:
        raise SystemExit("\n".join(failures))


if __name__ == "__main__":
    main()
