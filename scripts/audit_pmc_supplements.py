#!/usr/bin/env python3
"""Verify downloaded PMC AWS supplement files against the published MD5 values."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import parse_qs, urlparse


def file_hash(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-dir", type=Path, required=True)
    parser.add_argument("--supplement-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows: list[dict[str, str | int]] = []
    for metadata_path in sorted(args.metadata_dir.glob("PMC*.json")):
        if metadata_path.name.endswith(".meta.json"):
            continue
        payload = json.loads(metadata_path.read_text(encoding="utf-8"))
        for media_url in payload.get("media_urls", []):
            parsed = urlparse(media_url)
            relative = parsed.path.lstrip("/")
            if not relative.lower().endswith((".pdf", ".docx", ".xlsx", ".csv", ".tsv")):
                continue
            local_path = args.supplement_root / relative
            expected_md5 = parse_qs(parsed.query).get("md5", [""])[0]
            exists = local_path.exists()
            observed_md5 = file_hash(local_path, "md5") if exists else ""
            rows.append(
                {
                    "pmcid_version": relative.split("/", 1)[0],
                    "file_name": Path(relative).name,
                    "source_s3_url": media_url,
                    "local_path": local_path.as_posix(),
                    "bytes": local_path.stat().st_size if exists else 0,
                    "expected_md5": expected_md5,
                    "observed_md5": observed_md5,
                    "md5_status": "PASS" if exists and expected_md5 == observed_md5 else "NOT_DOWNLOADED" if not exists else "FAIL",
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "pmcid_version",
        "file_name",
        "source_s3_url",
        "local_path",
        "bytes",
        "expected_md5",
        "observed_md5",
        "md5_status",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
