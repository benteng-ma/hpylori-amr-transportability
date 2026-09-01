#!/usr/bin/env python3
"""Map one public read pair and audit coverage at canonical AMR loci."""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


MARKERS = [
    ("23S_copy1_A2142", "NC_000915.1", 447394),
    ("23S_copy1_A2143", "NC_000915.1", 447395),
    ("23S_copy2_A2142", "NC_000915.1", 1474748),
    ("23S_copy2_A2143", "NC_000915.1", 1474747),
    ("gyrA_N87_codon1", "NC_000915.1", 752770),
    ("gyrA_A88_codon1", "NC_000915.1", 752773),
    ("gyrA_D91_codon1", "NC_000915.1", 752782),
]


def run(command: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, text=True, **kwargs)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--bam-dir", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--bwa", default="bwa")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--fastp", default="fastp")
    args = parser.parse_args()

    args.bam_dir.mkdir(parents=True, exist_ok=True)
    args.result_dir.mkdir(parents=True, exist_ok=True)
    bam = args.bam_dir / "S27.sorted.bam"

    with tempfile.TemporaryDirectory(prefix="hpylori_fastp_smoke_") as temp_dir:
        run([
            args.fastp, "--in1", str(args.read1), "--in2", str(args.read2),
            "--out1", str(Path(temp_dir) / "discard_R1.fastq.gz"),
            "--out2", str(Path(temp_dir) / "discard_R2.fastq.gz"),
            "--disable_adapter_trimming", "--disable_quality_filtering", "--disable_length_filtering",
            "--json", str(args.result_dir / "S27_fastp.json"), "--html", str(args.result_dir / "S27_fastp.html"),
            "--thread", "4",
        ], capture_output=True)

    if not Path(str(args.reference) + ".bwt").exists():
        run([args.bwa, "index", str(args.reference)], capture_output=True)
    bwa_process = subprocess.Popen(
        [args.bwa, "mem", "-t", "4", str(args.reference), str(args.read1), str(args.read2)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False,
    )
    sort_process = subprocess.run(
        [args.samtools, "sort", "-@", "2", "-o", str(bam), "-"],
        stdin=bwa_process.stdout, capture_output=True,
    )
    assert bwa_process.stdout is not None
    bwa_process.stdout.close()
    bwa_stderr = bwa_process.stderr.read().decode(errors="replace") if bwa_process.stderr else ""
    bwa_returncode = bwa_process.wait()
    if bwa_returncode != 0 or sort_process.returncode != 0:
        raise SystemExit(f"mapping failed: bwa={bwa_returncode}; sort={sort_process.returncode}; {bwa_stderr}")
    run([args.samtools, "index", str(bam)], capture_output=True)

    flagstat = run([args.samtools, "flagstat", str(bam)], capture_output=True).stdout
    (args.result_dir / "S27_flagstat.txt").write_text(flagstat, encoding="utf-8")

    rows = []
    for marker, sequence, position in MARKERS:
        output = run(
            [args.samtools, "mpileup", "-aa", "-f", str(args.reference), "-r", f"{sequence}:{position}-{position}", str(bam)],
            capture_output=True,
        ).stdout.strip()
        fields = output.split("\t") if output else []
        depth = int(fields[3]) if len(fields) > 3 else 0
        rows.append({
            "sample": "S27", "marker": marker, "reference_sequence": sequence, "reference_position": position,
            "reference_base": fields[2] if len(fields) > 2 else "", "depth": depth,
            "pileup_bases": fields[4] if len(fields) > 4 else "",
            "status": "PASS" if depth >= 10 else "LOW_COVERAGE" if depth > 0 else "NO_COVERAGE",
        })
    with (args.result_dir / "S27_marker_pileup.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
