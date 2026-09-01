#!/usr/bin/env python3
"""Run CheckM2 genome completeness/contamination and merge its QC gates."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def resolve_database(path: Path) -> Path:
    """Resolve either a CheckM2 .dmnd file or its download parent directory."""
    candidate = path.resolve()
    if candidate.is_file():
        return candidate
    nested = candidate / "CheckM2_database" / "uniref100.KO.1.dmnd"
    if nested.is_file():
        return nested
    direct = candidate / "uniref100.KO.1.dmnd"
    if direct.is_file():
        return direct
    raise FileNotFoundError(f"CheckM2 DIAMOND database not found under {candidate}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=6)
    args = parser.parse_args()
    root = args.root.resolve()
    qc_rows = read_csv(root / "results/qc/assembly_qc.csv")
    output_dir = root / "data/processed/checkm2"
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    database = resolve_database(args.database)
    # CheckM2 1.1.0 constructs Prodigal and DIAMOND commands as unquoted
    # shell strings. Stage only execution paths under the system temp
    # directory so a project path containing spaces cannot split arguments.
    with tempfile.TemporaryDirectory(prefix="hpylori_checkm2_") as runtime_text:
        runtime = Path(runtime_text)
        input_dir = runtime / "input"
        runtime_output = runtime / "output"
        input_dir.mkdir()
        for row in qc_rows:
            source = root / row["assembly_path"]
            target = input_dir / f"{row['dataset_id']}__{row['isolate_id']}.fna"
            shutil.copy2(source, target)
        runtime_database = runtime / "database.dmnd"
        runtime_database.symlink_to(database)
        command = [
            sys.executable, "-m", "checkm2.main", "predict", "--threads", str(args.threads),
            "--input", str(input_dir), "--output-directory", str(runtime_output),
            "--database_path", str(runtime_database), "--extension", "fna", "--lowmem", "--force",
        ]
        completed = subprocess.run(command, capture_output=True, text=True)
        if completed.returncode == 0:
            shutil.copytree(runtime_output, output_dir, dirs_exist_ok=True)
    log = root / "results/logs/checkm2.log"
    log.write_text(completed.stdout + "\n--- STDERR ---\n" + completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise SystemExit(f"CheckM2 failed with exit {completed.returncode}")
    with (output_dir / "quality_report.tsv").open(newline="", encoding="utf-8-sig") as handle:
        quality = list(csv.DictReader(handle, delimiter="\t"))
    by_name = {row["Name"]: row for row in quality}
    merged = []
    for row in qc_rows:
        name = f"{row['dataset_id']}__{row['isolate_id']}"
        checkm = by_name.get(name)
        if checkm is None:
            raise ValueError(f"no CheckM2 row for {name}")
        completeness = float(checkm["Completeness"])
        contamination = float(checkm["Contamination"])
        row["completeness_percent"] = completeness
        row["contamination_percent"] = contamination
        row["completeness_gate"] = "PASS" if completeness >= 90 else "FAIL"
        row["contamination_gate"] = "PASS" if contamination <= 10 else "FAIL"
        row["final_qc_status"] = "PASS" if row["basic_qc_status"] == "PASS" and completeness >= 90 and contamination <= 10 else "FAIL"
        merged.append(row)
    output = root / "results/qc/assembly_qc_with_checkm2.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(merged[0].keys()))
        writer.writeheader()
        writer.writerows(merged)


if __name__ == "__main__":
    main()
