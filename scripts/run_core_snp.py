#!/usr/bin/env python3
"""Run phenotype-blind Snippy core-SNP processing on all QC-passing assemblies."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def call_one(
    row: dict[str, str], root: Path, snippy: str, reference: Path,
    tool_prefix: Path | None,
) -> dict[str, str]:
    sample = f"{row['dataset_id']}__{row['isolate_id']}"
    output = root / "data/phylogeny/snippy" / sample
    expected = output / "snps.vcf"
    aligned = output / "snps.aligned.fa"
    if expected.exists() and expected.stat().st_size > 0 and aligned.exists() and aligned.stat().st_size > 0:
        return {"dataset_id": row["dataset_id"], "isolate_id": row["isolate_id"], "sample_name": sample, "status": "already_complete", "output_dir": output.relative_to(root).as_posix(), "error": ""}
    environment = os.environ.copy()
    if tool_prefix is not None:
        environment["PATH"] = os.pathsep.join([
            str(Path(snippy).resolve().parent), str(tool_prefix / "bin"),
            environment.get("PATH", ""),
        ])
    internal_log = ""
    with tempfile.TemporaryDirectory(prefix=f"hpylori_snippy_{sample}_") as temporary:
        staged_output = Path(temporary) / "output"
        command = [
            snippy, "--outdir", str(staged_output), "--ctgs", str(root / row["assembly_path"]), "--ref", str(reference),
            "--cpus", "1", "--mincov", "10", "--minqual", "20", "--mapqual", "20", "--basequal", "20",
            "--fbopt", "-i -X -u --haplotype-length 0", "--rgid", sample, "--force",
        ]
        completed = subprocess.run(command, capture_output=True, text=True, env=environment)
        staged_log = staged_output / "snps.log"
        if staged_log.exists():
            internal_log = staged_log.read_text(encoding="utf-8", errors="replace")
        staged_expected = staged_output / "snps.vcf"
        staged_aligned = staged_output / "snps.aligned.fa"
        if completed.returncode == 0 and staged_expected.exists() and staged_aligned.exists():
            if output.exists():
                shutil.rmtree(output)
            output.mkdir(parents=True)
            retained = {
                "snps.aligned.fa", "snps.bam", "snps.bam.bai", "snps.bed",
                "snps.filt.vcf", "snps.log", "snps.raw.vcf", "snps.tab",
                "snps.txt", "snps.vcf",
            }
            for name in retained:
                source = staged_output / name
                if source.exists():
                    shutil.copy2(source, output / name)
    log = root / "results/logs/snippy" / f"{sample}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(
        completed.stdout + "\n--- STDERR ---\n" + completed.stderr
        + "\n--- SNIPPY INTERNAL LOG ---\n" + internal_log,
        encoding="utf-8",
    )
    if completed.returncode != 0 or not expected.exists() or not aligned.exists():
        return {"dataset_id": row["dataset_id"], "isolate_id": row["isolate_id"], "sample_name": sample, "status": "FAIL", "output_dir": output.relative_to(root).as_posix(), "error": f"snippy exit {completed.returncode}"}
    return {"dataset_id": row["dataset_id"], "isolate_id": row["isolate_id"], "sample_name": sample, "status": "PASS", "output_dir": output.relative_to(root).as_posix(), "error": ""}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--snippy", default="snippy")
    parser.add_argument("--snippy-core", default="snippy-core")
    parser.add_argument("--tool-prefix", type=Path, help="source conda prefix when using a no-space Snippy runtime")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--calls-only", action="store_true", help="materialize per-isolate calls without building a provisional core alignment")
    parser.add_argument("--reference", type=Path, default=Path("data/interim/reference/ncbi_dataset/data/GCF_000008525.1/genomic.gbff"))
    args = parser.parse_args()
    root = args.root.resolve()
    tool_prefix = args.tool_prefix.resolve() if args.tool_prefix else None
    reference = args.reference if args.reference.is_absolute() else root / args.reference
    final_qc_path = root / "results/qc/assembly_qc_with_checkm2.csv"
    qc_path = final_qc_path if final_qc_path.exists() else root / "results/qc/assembly_qc.csv"
    qc = [row for row in read_csv(qc_path) if row.get("final_qc_status", row["basic_qc_status"]) == "PASS"]
    statuses: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(call_one, row, root, args.snippy, reference, tool_prefix): row
            for row in qc
        }
        for index, future in enumerate(as_completed(futures), start=1):
            status = future.result()
            statuses.append(status)
            print(f"[{index}/{len(qc)}] {status['status']} {status['sample_name']}", flush=True)
    statuses.sort(key=lambda row: row["sample_name"])
    status_path = root / "results/qc/core_snp_processing_status.csv"
    status_path.parent.mkdir(parents=True, exist_ok=True)
    with status_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(statuses[0].keys()))
        writer.writeheader()
        writer.writerows(statuses)
    failures = [status for status in statuses if status["status"] == "FAIL"]
    if failures:
        raise SystemExit(f"{len(failures)} Snippy samples failed")
    if args.calls_only:
        return
    directories = [str(root / status["output_dir"]) for status in statuses]
    prefix = root / "data/phylogeny/core/core"
    prefix.parent.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    if tool_prefix is not None:
        environment["PATH"] = os.pathsep.join([
            str(Path(args.snippy_core).resolve().parent), str(tool_prefix / "bin"),
            environment.get("PATH", ""),
        ])
    # snippy-core 4.6.0 builds its final snp-sites command as an unquoted
    # shell string. Stage the prefix in a guaranteed no-space directory so
    # the canonical project path cannot be split at the final extraction
    # step. Publish no partial core artifacts until the staged run succeeds.
    with tempfile.TemporaryDirectory(prefix="hpylori_snippy_core_") as temporary:
        staged_prefix = Path(temporary) / "core"
        command = [args.snippy_core, "--ref", str(reference), "--prefix", str(staged_prefix), *directories]
        completed = subprocess.run(command, capture_output=True, text=True, env=environment)
        (root / "results/logs/snippy_core.log").write_text(
            completed.stdout + "\n--- STDERR ---\n" + completed.stderr,
            encoding="utf-8",
        )
        staged_alignment = Path(str(staged_prefix) + ".aln")
        if completed.returncode != 0 or not staged_alignment.exists() or staged_alignment.stat().st_size == 0:
            raise SystemExit(f"snippy-core failed with exit {completed.returncode}")
        for staged in Path(temporary).glob("core.*"):
            shutil.copy2(staged, prefix.parent / staged.name)
    if not Path(str(prefix) + ".aln").exists() or Path(str(prefix) + ".aln").stat().st_size == 0:
        raise SystemExit("staged snippy-core outputs were not published")


if __name__ == "__main__":
    main()
