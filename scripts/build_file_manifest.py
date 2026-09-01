#!/usr/bin/env python3
"""Hash auditable project artifacts while excluding environments and bulk sequence data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


INCLUDE_ROOT_FILES = {
    ".gitignore", "CITATION.cff", "DECISION_LOG.md", "LICENSE", "README.md", "STATUS.md",
    "Snakefile", "_quarto.yml", "_targets.R", "environment.yml", "renv.lock",
}
INCLUDE_DIRS = {"config", "literature", "metadata", "reports", "scripts", "tests", "workflow"}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def role(relative: Path) -> str:
    if relative.parts[0] == "literature":
        return "immutable public literature input" if "raw" in relative.parts or "supplements" in relative.parts else "literature audit"
    if relative.parts[0] == "metadata":
        return "machine-readable manifest or derived audit table"
    if relative.parts[0] == "reports":
        return "human-readable audit report"
    if relative.parts[0] in {"scripts", "tests", "workflow", "config"}:
        return "reproducible analysis control"
    if relative.parts[0] == "results":
        return "Phase 0 smoke result"
    return "repository control"


def source_url(path: Path) -> str:
    meta_path = path.with_suffix(path.suffix + ".meta.json")
    if meta_path.exists():
        try:
            return str(json.loads(meta_path.read_text(encoding="utf-8")).get("source_url", ""))
        except (json.JSONDecodeError, OSError):
            return ""
    return ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=Path("metadata/file_manifest.csv"))
    args = parser.parse_args()
    root = args.root.resolve()
    output = (root / args.output).resolve() if not args.output.is_absolute() else args.output

    files = [root / name for name in INCLUDE_ROOT_FILES if (root / name).is_file()]
    for directory in sorted(INCLUDE_DIRS):
        files.extend(path for path in (root / directory).rglob("*") if path.is_file())
    for directory in (root / "results" / "smoke", root / "results" / "logs"):
        if directory.exists():
            files.extend(path for path in directory.rglob("*") if path.is_file())
    files = [path for path in files if path.resolve() != output]

    rows = []
    for path in sorted(files):
        relative = path.relative_to(root)
        rows.append({
            "relative_path": relative.as_posix(), "file_role": role(relative), "source_url": source_url(path),
            "accessed_at": "2026-08-30", "bytes": path.stat().st_size, "sha256": digest(path),
            "immutable": "yes" if role(relative).startswith("immutable") else "no",
            "notes": "generated before Phase 1 full performance inspection" if relative.parts[0] in {"metadata", "reports", "results"} else "",
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "file_role", "source_url", "accessed_at", "bytes", "sha256", "immutable", "notes"])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
