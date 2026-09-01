#!/usr/bin/env python3
"""Phase 0 smoke caller for canonical 23S rRNA and gyrA AMR markers.

This is deliberately a targeted extraction audit, not the final AMR model.
Common 23S labels A2142/A2143 map to U27270.1 feature positions 2143/2144
(the one-base convention difference used explicitly by the HpGP code). gyrA
numbering follows the translated 26695 CDS from GCF_000008525.1.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


DEFAULT_MIN_IDENTITY = 0.90


def extract_reference_features(u27270_gb: Path, hp26695_gb: Path, output_dir: Path) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    u_record = SeqIO.read(u27270_gb, "genbank")
    rrna = next(
        feature for feature in u_record.features
        if feature.type == "rRNA" and "23S" in " ".join(feature.qualifiers.get("product", []))
    )
    rrna_record = rrna.extract(u_record)
    rrna_record.id = "U27270.1_23S_rRNA"
    rrna_record.description = "23S rRNA feature; common A2142/A2143 labels map to feature positions 2143/2144"
    rrna_path = output_dir / "U27270.1_23S_rRNA.fasta"
    SeqIO.write(rrna_record, rrna_path, "fasta")

    hp_record = SeqIO.read(hp26695_gb, "genbank")
    gyra = next(
        feature for feature in hp_record.features
        if feature.type == "CDS" and "gyrA" in feature.qualifiers.get("gene", [])
    )
    gyra_record = gyra.extract(hp_record)
    gyra_record.id = "GCF_000008525.1_26695_gyrA"
    gyra_record.description = "gyrA CDS; translated residue numbering"
    gyra_path = output_dir / "GCF_000008525.1_26695_gyrA.fasta"
    SeqIO.write(gyra_record, gyra_path, "fasta")
    return rrna_path, gyra_path


def blast(query: Path, subject: Path, blastn: str, task: str = "megablast") -> list[dict[str, str]]:
    columns = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "qseq", "sseq"]
    command = [
        blastn, "-task", task, "-query", str(query), "-subject", str(subject),
        "-dust", "no", "-max_target_seqs", "20", "-outfmt", "6 " + " ".join(columns),
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    return [dict(zip(columns, line.split("\t"))) for line in completed.stdout.splitlines() if line]


def query_to_subject_map(hit: dict[str, str]) -> dict[int, str]:
    query_position = int(hit["qstart"]) - 1
    mapping: dict[int, str] = {}
    for query_base, subject_base in zip(hit["qseq"], hit["sseq"]):
        if query_base != "-":
            query_position += 1
            if subject_base != "-":
                mapping[query_position] = subject_base.upper()
    return mapping


def query_to_subject_position_map(hit: dict[str, str]) -> dict[int, int]:
    query_position = int(hit["qstart"]) - 1
    subject_position = int(hit["sstart"])
    subject_step = 1 if int(hit["send"]) >= int(hit["sstart"]) else -1
    mapping: dict[int, int] = {}
    for query_base, subject_base in zip(hit["qseq"], hit["sseq"]):
        if query_base != "-":
            query_position += 1
            if subject_base != "-":
                mapping[query_position] = subject_position
        if subject_base != "-":
            subject_position += subject_step
    return mapping


def best_nonredundant_hits(
    hits: list[dict[str, str]],
    minimum_coverage: float,
    query_length: int,
    minimum_identity: float = DEFAULT_MIN_IDENTITY,
) -> list[dict[str, str]]:
    eligible = [
        hit for hit in hits
        if int(hit["length"]) / query_length >= minimum_coverage
        and float(hit["pident"]) / 100 >= minimum_identity
    ]
    eligible.sort(key=lambda hit: (float(hit["pident"]), int(hit["length"])), reverse=True)
    retained: list[dict[str, str]] = []
    for hit in eligible:
        locus = (hit["sseqid"], min(int(hit["sstart"]), int(hit["send"])), max(int(hit["sstart"]), int(hit["send"])))
        if any(
            locus[0] == kept[0] and min(locus[2], kept[2]) - max(locus[1], kept[1]) > 0.8 * min(locus[2] - locus[1], kept[2] - kept[1])
            for kept in [(item["sseqid"], min(int(item["sstart"]), int(item["send"])), max(int(item["sstart"]), int(item["send"]))) for item in retained]
        ):
            continue
        retained.append(hit)
    return retained


def call_one(
    assembly: Path,
    rrna_path: Path,
    gyra_path: Path,
    blastn: str,
    minimum_identity: float = DEFAULT_MIN_IDENTITY,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    rrna_reference = str(SeqIO.read(rrna_path, "fasta").seq).upper()
    gyra_reference = str(SeqIO.read(gyra_path, "fasta").seq).upper()

    rrna_hits = best_nonredundant_hits(
        blast(rrna_path, assembly, blastn), 0.05, len(rrna_reference), minimum_identity
    )
    rrna_hits = [hit for hit in rrna_hits if all(position in query_to_subject_map(hit) for position in (2143, 2144))]
    for copy_index, hit in enumerate(rrna_hits, start=1):
        mapping = query_to_subject_map(hit)
        position_mapping = query_to_subject_position_map(hit)
        for position, query_position in ((2142, 2143), (2143, 2144)):
            observed = mapping.get(query_position, "")
            reference = rrna_reference[query_position - 1]
            known = (position == 2142 and observed in {"G", "C"}) or (position == 2143 and observed == "G")
            rows.append({
                "assembly": assembly.name, "gene": "23S_rRNA", "copy": copy_index,
                "coordinate_system": "common_H_pylori_label; U27270_feature_position=label+1", "position": position,
                "subject_sequence": hit["sseqid"], "subject_position": position_mapping.get(query_position, ""),
                "reference": reference, "observed": observed,
                "change": f"{reference}{position}{observed}" if observed else "uncallable",
                "known_resistance_marker": "yes" if known else "no",
                "hit_identity_percent": hit["pident"], "hit_coverage": round(int(hit["length"]) / len(rrna_reference), 5),
                "status": "PASS" if observed else "UNCALLABLE",
            })
    if not rrna_hits:
        for position in (2142, 2143):
            rows.append({
                "assembly": assembly.name, "gene": "23S_rRNA", "copy": "", "coordinate_system": "common_H_pylori_label; U27270_feature_position=label+1",
                "position": position, "subject_sequence": "", "subject_position": "", "reference": "A", "observed": "", "change": "uncallable", "known_resistance_marker": "",
                "hit_identity_percent": "", "hit_coverage": "", "status": "UNCALLABLE_NO_TARGET_SPANNING_HIT",
            })

    gyra_hits = best_nonredundant_hits(
        blast(gyra_path, assembly, blastn), 0.90, len(gyra_reference), minimum_identity
    )
    if gyra_hits:
        hit = gyra_hits[0]
        mapping = query_to_subject_map(hit)
        position_mapping = query_to_subject_position_map(hit)
        for amino_acid_position in (87, 88, 91):
            nucleotide_positions = range((amino_acid_position - 1) * 3 + 1, amino_acid_position * 3 + 1)
            observed_codon = "".join(mapping.get(position, "") for position in nucleotide_positions)
            reference_codon = gyra_reference[(amino_acid_position - 1) * 3:amino_acid_position * 3]
            reference_aa = str(Seq(reference_codon).translate())
            observed_aa = str(Seq(observed_codon).translate()) if len(observed_codon) == 3 else ""
            known_changes = {87: {"K", "I"}, 88: {"V", "P"}, 91: {"G", "N", "Y"}}
            known = observed_aa in known_changes[amino_acid_position]
            rows.append({
                "assembly": assembly.name, "gene": "gyrA", "copy": 1,
                "coordinate_system": "26695_gyrA_protein", "position": amino_acid_position,
                "subject_sequence": hit["sseqid"],
                "subject_position": f"{position_mapping.get(min(nucleotide_positions), '')}-{position_mapping.get(max(nucleotide_positions), '')}",
                "reference": reference_aa, "observed": observed_aa,
                "change": f"{reference_aa}{amino_acid_position}{observed_aa}" if observed_aa else "uncallable",
                "known_resistance_marker": "yes" if known else "no",
                "hit_identity_percent": hit["pident"], "hit_coverage": round(int(hit["length"]) / len(gyra_reference), 5),
                "status": "PASS" if observed_aa else "UNCALLABLE",
            })
    else:
        for amino_acid_position, reference_aa in ((87, "N"), (88, "A"), (91, "D")):
            rows.append({
                "assembly": assembly.name, "gene": "gyrA", "copy": "", "coordinate_system": "26695_gyrA_protein",
                "position": amino_acid_position, "subject_sequence": "", "subject_position": "", "reference": reference_aa, "observed": "", "change": "uncallable", "known_resistance_marker": "",
                "hit_identity_percent": "", "hit_coverage": "", "status": "UNCALLABLE_NO_FULL_LENGTH_HIT",
            })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--u27270-genbank", type=Path, required=True)
    parser.add_argument("--hp26695-genbank", type=Path, required=True)
    parser.add_argument("--reference-output-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--minimum-identity", type=float, default=DEFAULT_MIN_IDENTITY)
    parser.add_argument("assemblies", nargs="+", type=Path)
    args = parser.parse_args()

    rrna_path, gyra_path = extract_reference_features(args.u27270_genbank, args.hp26695_genbank, args.reference_output_dir)
    rows = [
        row for assembly in args.assemblies
        for row in call_one(assembly, rrna_path, gyra_path, args.blastn, args.minimum_identity)
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["assembly", "gene", "copy", "coordinate_system", "position", "subject_sequence", "subject_position", "reference", "observed", "change", "known_resistance_marker", "hit_identity_percent", "hit_coverage", "status"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
