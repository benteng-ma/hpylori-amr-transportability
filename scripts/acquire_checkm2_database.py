#!/usr/bin/env python3
"""Acquire and verify the frozen CheckM2 v3 reference database."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tarfile
from datetime import datetime, timezone
from pathlib import Path


RECORD_ID = "14897628"
ARCHIVE_URL = "https://zenodo.org/api/records/14897628/files/checkm2_database.tar.gz/content"
EXPECTED_BYTES = 1_735_095_710
EXPECTED_MD5 = "07c10655620843b517d0df0c160d911f"


def digest(path: Path, algorithm: str) -> str:
    checksum = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def validate_archive(path: Path) -> None:
    if not path.is_file() or path.stat().st_size != EXPECTED_BYTES:
        observed = path.stat().st_size if path.exists() else "missing"
        raise ValueError(f"CheckM2 archive size mismatch: {observed} != {EXPECTED_BYTES}")
    if digest(path, "md5") != EXPECTED_MD5:
        raise ValueError("CheckM2 archive MD5 mismatch")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--curl", default="curl.exe" if shutil.which("curl.exe") else "curl")
    parser.add_argument("--output-dir", type=Path, default=Path("data/raw/checkm2_db"))
    args = parser.parse_args()
    root = args.root.resolve()
    output_dir = args.output_dir if args.output_dir.is_absolute() else root / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    archive = output_dir / "checkm2_database.tar.gz"
    partial = archive.with_suffix(archive.suffix + ".part")
    if archive.exists():
        validate_archive(archive)
    else:
        completed = subprocess.run([
            args.curl, "--location", "--fail", "--show-error", "--retry", "8", "--retry-all-errors",
            "--connect-timeout", "30", "--continue-at", "-", "--output", str(partial), ARCHIVE_URL,
        ])
        if completed.returncode != 0:
            raise SystemExit(f"curl failed with exit {completed.returncode}")
        validate_archive(partial)
        partial.replace(archive)

    database = output_dir / "CheckM2_database" / "uniref100.KO.1.dmnd"
    if not database.is_file():
        with tarfile.open(archive, "r:gz") as handle:
            handle.extractall(output_dir, filter="data")
    if not database.is_file() or database.stat().st_size == 0:
        raise FileNotFoundError(f"extracted CheckM2 database not found: {database}")
    payload = {
        "record_id": RECORD_ID,
        "doi": f"10.5281/zenodo.{RECORD_ID}",
        "archive_url": ARCHIVE_URL,
        "archive_path": archive.relative_to(root).as_posix(),
        "archive_bytes": archive.stat().st_size,
        "archive_md5": EXPECTED_MD5,
        "archive_sha256": digest(archive, "sha256"),
        "database_path": database.relative_to(root).as_posix(),
        "database_bytes": database.stat().st_size,
        "verified_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    (output_dir / "complete.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
