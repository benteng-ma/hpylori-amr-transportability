#!/usr/bin/env python3
"""QC, assemble, and call frozen AMR markers from Zenodo paired reads."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import tempfile
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Lock

from Bio import SeqIO
from Bio.Seq import Seq


DEPTH_MIN = 10
BASE_QUALITY_MIN = 20
MAPPING_QUALITY_MIN = 20
MAJOR_ALLELE_MIN = 0.80
MIXED_ALLELE_MIN = 0.20
ASSEMBLY_LOCK = Lock()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def run(command: list[str], **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(command, check=True, **kwargs)


def extract_pair(archive_path: Path, destination: Path) -> tuple[Path, Path]:
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path) as archive:
        members = [name for name in archive.namelist() if not name.endswith("/") and name.lower().endswith((".fastq.gz", ".fq.gz"))]
        read1 = [name for name in members if re.search(r"(?:_|\.)R?1(?:_|\.|$)", Path(name).name, re.I)]
        read2 = [name for name in members if re.search(r"(?:_|\.)R?2(?:_|\.|$)", Path(name).name, re.I)]
        if len(read1) != 1 or len(read2) != 1:
            raise ValueError(f"cannot identify one paired read set in {archive_path}: {members}")
        outputs = []
        for index, member in enumerate((read1[0], read2[0]), start=1):
            output = destination / f"raw_R{index}.fastq.gz"
            with archive.open(member) as source, output.open("wb") as target:
                shutil.copyfileobj(source, target, length=4 * 1024 * 1024)
            outputs.append(output)
    return outputs[0], outputs[1]


def parse_pileup_bases(text: str, reference: str) -> dict[str, int]:
    counts = {base: 0 for base in "ACGT"}
    index = 0
    while index < len(text):
        token = text[index]
        if token == "^":
            index += 2
            continue
        if token == "$":
            index += 1
            continue
        if token in "+-":
            match = re.match(r"[+-](\d+)", text[index:])
            if match is None:
                index += 1
                continue
            index += len(match.group(0)) + int(match.group(1))
            continue
        if token in ".,":
            counts[reference.upper()] += 1
        elif token.upper() in counts:
            counts[token.upper()] += 1
        index += 1
    return counts


def pileup_one(samtools: str, reference: Path, bam: Path, sequence: str, position: int) -> dict[str, object]:
    completed = run([
        samtools, "mpileup", "-aa", "-q", str(MAPPING_QUALITY_MIN), "-Q", str(BASE_QUALITY_MIN),
        "-f", str(reference), "-r", f"{sequence}:{position}-{position}", str(bam),
    ], capture_output=True, text=True)
    fields = completed.stdout.strip().split("\t") if completed.stdout.strip() else []
    reference_base = fields[2].upper() if len(fields) > 2 else ""
    counts = parse_pileup_bases(fields[4], reference_base) if len(fields) > 4 and reference_base in "ACGT" else {base: 0 for base in "ACGT"}
    depth = sum(counts.values())
    major_base = max(counts, key=counts.get) if depth else ""
    major_fraction = counts[major_base] / depth if depth else 0.0
    return {
        "reference_base": reference_base, "depth": depth, **{f"count_{base}": counts[base] for base in "ACGT"},
        "major_base": major_base, "major_fraction": round(major_fraction, 6),
        "callability": "PASS" if depth >= DEPTH_MIN and major_fraction >= MAJOR_ALLELE_MIN else "MIXED" if depth >= DEPTH_MIN else "LOW_COVERAGE",
    }


def map_reads(bwa: str, samtools: str, reference: Path, read1: Path, read2: Path, bam: Path, threads: int, log: Path) -> None:
    bam.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    # BWA emits progress records throughout long, high-depth alignments. Writing
    # stderr directly prevents the fixed-size PIPE buffer from filling while
    # samtools is still consuming stdout.
    with log.open("wb") as log_handle:
        bwa_process = subprocess.Popen(
            [bwa, "mem", "-t", str(threads), str(reference), str(read1), str(read2)],
            stdout=subprocess.PIPE, stderr=log_handle,
        )
        assert bwa_process.stdout is not None
        sort_result = subprocess.run([samtools, "sort", "-@", "2", "-o", str(bam), "-"], stdin=bwa_process.stdout, capture_output=True)
        bwa_process.stdout.close()
        bwa_code = bwa_process.wait()
        log_handle.write(b"\n--- samtools sort ---\n" + sort_result.stderr)
    if bwa_code != 0 or sort_result.returncode != 0:
        raise RuntimeError(f"mapping failure bwa={bwa_code} samtools={sort_result.returncode}")
    run([samtools, "index", str(bam)], capture_output=True)


def gyr_positions(genbank: Path) -> tuple[str, int, dict[int, list[tuple[int, str]]]]:
    record = SeqIO.read(genbank, "genbank")
    feature = next(feature for feature in record.features if feature.type == "CDS" and "gyrA" in feature.qualifiers.get("gene", []))
    coding = str(feature.extract(record).seq).upper()
    result: dict[int, list[tuple[int, str]]] = {}
    start = int(feature.location.start)
    end = int(feature.location.end)
    strand = int(feature.location.strand or 1)
    for residue in (87, 88, 91):
        values = []
        for coding_index in range((residue - 1) * 3, residue * 3):
            genomic = start + coding_index + 1 if strand == 1 else end - coding_index
            values.append((genomic, coding[coding_index]))
        result[residue] = values
    return record.id, strand, result


def summarize_calls(sample_id: str, copy_calls: list[dict[str, object]], pooled_calls: list[dict[str, object]], gyr_calls: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    pooled_by_label = {str(row["label"]): row for row in pooled_calls}
    resistant_fractions = []
    for label, resistant in (("A2142", {"G", "C"}), ("A2143", {"G"})):
        call = pooled_by_label[label]
        depth = int(call["depth"])
        resistant_count = sum(int(call[f"count_{base}"]) for base in resistant)
        resistant_fraction = resistant_count / depth if depth else 0.0
        resistant_fractions.append(resistant_fraction)
    pooled_callable = all(int(row["depth"]) >= DEPTH_MIN for row in pooled_calls)
    if pooled_callable and max(resistant_fractions) >= MIXED_ALLELE_MIN:
        clr_prediction = "R"
    elif pooled_callable and all(str(row["callability"]) == "PASS" for row in pooled_calls):
        clr_prediction = "S"
    else:
        clr_prediction = "UNCALLABLE"
    rows.append({
        "dataset_id": "ZENODO_10369064", "isolate_id": sample_id, "antibiotic": "clarithromycin",
        "prediction": clr_prediction, "sequence_support": "RAW_READ_POOLED_23S_WITH_COPY_AUDIT",
        "marker_summary": ";".join(f"{label}:{fraction:.4f}" for label, fraction in zip(("A2142", "A2143"), resistant_fractions)),
        "callable": "yes" if clr_prediction != "UNCALLABLE" else "no",
    })

    gyr_by_residue: dict[int, list[dict[str, object]]] = {}
    for call in gyr_calls:
        gyr_by_residue.setdefault(int(call["residue"]), []).append(call)
    changes = []
    all_callable = True
    known = {87: {"K", "I"}, 88: {"V", "P"}, 91: {"G", "N", "Y"}}
    for residue in (87, 88, 91):
        codon_calls = sorted(gyr_by_residue[residue], key=lambda row: int(row["coding_base_index"]))
        if len(codon_calls) != 3 or any(str(call["callability"]) != "PASS" for call in codon_calls):
            all_callable = False
            continue
        codon = "".join(str(call["coding_major_base"]) for call in codon_calls)
        amino_acid = str(Seq(codon).translate())
        reference_aa = {87: "N", 88: "A", 91: "D"}[residue]
        changes.append(f"{reference_aa}{residue}{amino_acid}")
    if all_callable:
        lvx_prediction = "R" if any(change[-1] in known[int(re.search(r"\d+", change).group())] for change in changes) else "S"
    else:
        lvx_prediction = "UNCALLABLE"
    rows.append({
        "dataset_id": "ZENODO_10369064", "isolate_id": sample_id, "antibiotic": "levofloxacin",
        "prediction": lvx_prediction, "sequence_support": "RAW_READ_SUPPORTED", "marker_summary": ";".join(changes),
        "callable": "yes" if lvx_prediction != "UNCALLABLE" else "no",
    })
    return rows


def process_sample(row: dict[str, str], root: Path, args: argparse.Namespace, reference: Path, rrna_reference: Path, gyr_sequence: str, gyr_strand: int, gyr_map: dict[int, list[tuple[int, str]]]) -> dict[str, object]:
    sample_id = row["isolate_id"]
    per_sample = root / "data/variants/zenodo_10369064" / sample_id
    done = per_sample / "complete.json"
    if done.exists() and (per_sample / "predictions.csv").exists():
        payload = json.loads(done.read_text(encoding="utf-8"))
        payload["processing_status"] = "already_complete"
        return payload
    per_sample.mkdir(parents=True, exist_ok=True)
    # High-depth archives create multi-gigabyte sorted BAM intermediates.
    # Keep execution-only files on the native WSL temporary filesystem and
    # write every auditable result back to the canonical project directory.
    with tempfile.TemporaryDirectory(prefix=f"hpylori_zenodo_{sample_id}_") as temporary_text:
        temporary = Path(temporary_text)
        raw1, raw2 = extract_pair(root / row["local_path"], temporary)
        clean1, clean2 = temporary / "clean_R1.fastq.gz", temporary / "clean_R2.fastq.gz"
        fastp_json = root / "results/qc/zenodo_fastp" / f"{sample_id}.json"
        fastp_html = root / "results/qc/zenodo_fastp" / f"{sample_id}.html"
        fastp_json.parent.mkdir(parents=True, exist_ok=True)
        run([
            args.fastp, "--in1", str(raw1), "--in2", str(raw2), "--out1", str(clean1), "--out2", str(clean2),
            "--detect_adapter_for_pe", "--cut_tail", "--cut_window_size", "4", "--cut_mean_quality", "20",
            "--length_required", "50", "--correction", "--thread", str(args.threads_per_sample),
            "--json", str(fastp_json), "--html", str(fastp_html),
        ], capture_output=True)

        whole_bam = temporary / "whole.sorted.bam"
        pooled_bam = temporary / "rrna.sorted.bam"
        map_reads(args.bwa, args.samtools, reference, clean1, clean2, whole_bam, args.threads_per_sample, root / "results/logs/zenodo" / f"{sample_id}.whole_bwa.log")
        map_reads(args.bwa, args.samtools, rrna_reference, clean1, clean2, pooled_bam, args.threads_per_sample, root / "results/logs/zenodo" / f"{sample_id}.23S_bwa.log")
        flagstat = run([args.samtools, "flagstat", str(whole_bam)], capture_output=True, text=True).stdout
        (per_sample / "flagstat.txt").write_text(flagstat, encoding="utf-8")

        copy_specs = [
            ("copy1", "A2142", 447394, {"G", "C"}), ("copy1", "A2143", 447395, {"G"}),
            ("copy2", "A2142", 1474748, {"C", "G"}), ("copy2", "A2143", 1474747, {"C"}),
        ]
        copy_calls = []
        for copy_id, label, position, resistant in copy_specs:
            call = pileup_one(args.samtools, reference, whole_bam, "NC_000915.1", position)
            call.update({"sample_id": sample_id, "assay": "23S_COPY_SPECIFIC", "copy": copy_id, "label": label, "position": position, "resistant_bases": "".join(sorted(resistant))})
            copy_calls.append(call)
        pooled_calls = []
        for label, position, resistant in (("A2142", 2143, {"G", "C"}), ("A2143", 2144, {"G"})):
            call = pileup_one(args.samtools, rrna_reference, pooled_bam, "U27270.1_23S_rRNA", position)
            call.update({"sample_id": sample_id, "assay": "23S_POOLED", "copy": "pooled", "label": label, "position": position, "resistant_bases": "".join(sorted(resistant))})
            pooled_calls.append(call)
        gyr_calls = []
        for residue, positions in gyr_map.items():
            for coding_index, (position, reference_coding_base) in enumerate(positions, start=1):
                call = pileup_one(args.samtools, reference, whole_bam, gyr_sequence, position)
                genomic_major = str(call["major_base"])
                coding_major = str(Seq(genomic_major).complement()) if genomic_major and gyr_strand == -1 else genomic_major
                call.update({
                    "sample_id": sample_id, "assay": "GYRA_CODON", "residue": residue,
                    "coding_base_index": coding_index, "position": position,
                    "reference_coding_base": reference_coding_base, "coding_major_base": coding_major,
                })
                gyr_calls.append(call)

        for name, rows in (("23s_copy_calls.csv", copy_calls), ("23s_pooled_calls.csv", pooled_calls), ("gyrA_calls.csv", gyr_calls)):
            with (per_sample / name).open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
                writer.writeheader()
                writer.writerows(rows)
        predictions = summarize_calls(sample_id, copy_calls, pooled_calls, gyr_calls)
        with (per_sample / "predictions.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(predictions[0].keys()))
            writer.writeheader()
            writer.writerows(predictions)

        assembly = root / "data/assemblies/ZENODO_10369064" / f"{sample_id}.fna"
        assembly.parent.mkdir(parents=True, exist_ok=True)
        # High-depth H. pylori read sets can each use most of a 16-GB worker
        # during SKESA graph construction. Keep assembly phenotype-blind and
        # parameter-identical, but serialize this memory-bound support step.
        with ASSEMBLY_LOCK:
            run([
                args.skesa, "--fastq", f"{clean1},{clean2}", "--cores", str(args.threads_per_sample),
                "--contigs_out", str(assembly),
            ], capture_output=True)
        if not assembly.exists() or assembly.stat().st_size == 0:
            raise RuntimeError("SKESA produced no assembly")

    payload = {"dataset_id": "ZENODO_10369064", "isolate_id": sample_id, "processing_status": "PASS", "assembly_path": assembly.relative_to(root).as_posix()}
    done.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def merge_outputs(root: Path, statuses: list[dict[str, object]]) -> None:
    predictions = []
    for status in statuses:
        path = root / "data/variants/zenodo_10369064" / str(status["isolate_id"]) / "predictions.csv"
        if path.exists():
            predictions.extend(read_csv(path))
    output = root / "results/panels/zenodo_read_marker_predictions.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    if predictions:
        predictions.sort(key=lambda row: (row["isolate_id"], row["antibiotic"]))
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(predictions[0].keys()))
            writer.writeheader()
            writer.writerows(predictions)
    status_output = root / "results/qc/zenodo_processing_status.csv"
    with status_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["dataset_id", "isolate_id", "processing_status", "assembly_path", "error"])
        writer.writeheader()
        for status in sorted(statuses, key=lambda row: str(row["isolate_id"])):
            writer.writerow({key: status.get(key, "") for key in writer.fieldnames})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--acquisition-manifest", type=Path, default=Path("metadata/phase2/acquisition_reads.csv"))
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--threads-per-sample", type=int, default=3)
    parser.add_argument("--fastp", default="fastp")
    parser.add_argument("--bwa", default="bwa")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--skesa", default="skesa")
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = args.acquisition_manifest if args.acquisition_manifest.is_absolute() else root / args.acquisition_manifest
    rows = [row for row in read_csv(manifest) if row["validation_status"] == "PASS"]
    reference = root / "data/interim/reference/ncbi_dataset/data/GCF_000008525.1/GCF_000008525.1_ASM852v1_genomic.fna"
    rrna_reference = root / "data/interim/reference/markers/U27270.1_23S_rRNA.fasta"
    for reference_path in (reference, rrna_reference):
        if not Path(str(reference_path) + ".bwt").exists():
            run([args.bwa, "index", str(reference_path)], capture_output=True)
    gyr_sequence, gyr_strand, gyr_map = gyr_positions(root / "data/interim/reference/ncbi_dataset/data/GCF_000008525.1/genomic.gbff")
    statuses: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_sample, row, root, args, reference, rrna_reference, gyr_sequence, gyr_strand, gyr_map): row for row in rows}
        for index, future in enumerate(as_completed(futures), start=1):
            row = futures[future]
            try:
                status = future.result()
            except Exception as exc:
                status = {"dataset_id": "ZENODO_10369064", "isolate_id": row["isolate_id"], "processing_status": "FAIL", "assembly_path": "", "error": str(exc)}
            statuses.append(status)
            print(f"[{index}/{len(rows)}] {status['processing_status']} {row['isolate_id']}", flush=True)
    merge_outputs(root, statuses)
    failures = [status for status in statuses if status["processing_status"] == "FAIL"]
    if failures:
        raise SystemExit(f"{len(failures)} Zenodo samples failed; see results/qc/zenodo_processing_status.csv")


if __name__ == "__main__":
    main()
