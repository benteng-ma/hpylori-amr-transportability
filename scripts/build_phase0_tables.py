#!/usr/bin/env python3
"""Build Phase 0 dataset, isolate, phenotype, and mapping audit tables."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

import pdfplumber


ASSEMBLY_HEADER = [
    "dataset_id", "strain", "biosample", "assembly_accession", "genbank_accession",
    "refseq_accession", "wgs_accession", "assembly_level", "coverage", "contig_n50",
    "total_length", "ftp_genbank", "ftp_refseq",
]


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def ncbi_rows(path: Path, dataset_id: str) -> list[dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    result = payload["summary"]["result"]
    rows: list[dict[str, object]] = []
    for uid in result.get("uids", []):
        record = result[uid]
        infraspecies = record.get("biosource", {}).get("infraspecieslist", [])
        strain = next((x.get("sub_value", "") for x in infraspecies if x.get("sub_type") == "strain"), "")
        synonym = record.get("synonym", {})
        meta = record.get("meta", "")
        length_match = re.search(r'category="total_length"[^>]*>(\d+)<', meta)
        rows.append(
            {
                "dataset_id": dataset_id,
                "strain": strain,
                "biosample": record.get("biosampleaccn", ""),
                "assembly_accession": record.get("assemblyaccession", ""),
                "genbank_accession": synonym.get("genbank", ""),
                "refseq_accession": synonym.get("refseq", ""),
                "wgs_accession": record.get("wgs", ""),
                "assembly_level": record.get("assemblystatus", ""),
                "coverage": record.get("coverage", ""),
                "contig_n50": record.get("contign50", ""),
                "total_length": length_match.group(1) if length_match else "",
                "ftp_genbank": record.get("ftppath_genbank", ""),
                "ftp_refseq": record.get("ftppath_refseq", ""),
            }
        )
    return rows


def extract_china60_mic(path: Path) -> list[dict[str, str]]:
    with pdfplumber.open(path) as document:
        lines = [line for index in range(3) for line in (document.pages[index].extract_text() or "").splitlines()]
    header = ["isolate_id", "gender", "age_years", "eradication_history", "clarithromycin", "levofloxacin", "amoxicillin", "furazolidone", "tetracycline", "metronidazole"]
    rows = []
    for line in lines:
        if not line.startswith("SHZY"):
            continue
        clean = [value.replace("＞", ">").replace("＜", "<").strip() for value in line.split()]
        if len(clean) != len(header):
            raise ValueError(f"unexpected China MIC row with {len(clean)} fields: {line}")
        rows.append(dict(zip(header, clean)))
    return rows


def extract_china60_published_calls(path: Path) -> list[dict[str, str]]:
    """Extract the publication's phenotype-resistance (P-R) calls from Table S9."""
    drug_columns = {
        "clarithromycin": 2,
        "levofloxacin": 4,
        "amoxicillin": 6,
        "tetracycline": 8,
        "metronidazole": 10,
        "furazolidone": 12,
        "multidrug_resistance": 14,
    }
    rows: list[dict[str, str]] = []
    with pdfplumber.open(path) as document:
        for page_index in (13, 14, 15):
            for table in document.pages[page_index].extract_tables():
                for values in table:
                    if not values or not re.fullmatch(r"SHZY\d{2}", values[0] or ""):
                        continue
                    if len(values) != 15:
                        raise ValueError(f"unexpected China Table S9 row: {values}")
                    row = {"isolate_id": values[0]}
                    row.update({drug: values[index] for drug, index in drug_columns.items()})
                    rows.append(row)
    if len(rows) != 60 or len({row["isolate_id"] for row in rows}) != 60:
        raise ValueError(f"expected 60 unique China Table S9 phenotype rows, found {len(rows)}")
    return rows


def numeric_mic(value: str) -> tuple[str, str]:
    value = value.strip()
    match = re.fullmatch(r"([<>]=?)?\s*(\d+(?:\.\d+)?)", value)
    return (match.group(2), match.group(1) or "=") if match else ("", "")


def breakpoint_call(drug: str, value: str) -> tuple[str, str]:
    numeric, operator = numeric_mic(value)
    if not numeric:
        return "", ""
    mic = float(numeric)
    resistant_threshold = {
        "clarithromycin": 0.5,
        "levofloxacin": 1.0,
        "amoxicillin": 0.125,
        "metronidazole": 8.0,
        "tetracycline": 1.0,
        "furazolidone": 4.0,
    }[drug]
    if drug == "clarithromycin" and mic == 0.5 and operator == "=":
        call = "I"
    elif drug == "furazolidone" and mic == resistant_threshold:
        call = "R"
    else:
        call = "R" if operator.startswith(">") or mic > resistant_threshold else "S"
    return call, "yes" if mic == resistant_threshold else "no"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.root.resolve()
    phase0 = root / "metadata" / "phase0"

    assembly_specs = [
        ("HPGP_GLOBAL", root / "metadata/raw/PRJNA529500_assembly_2026-08-30.json"),
        ("US_MAYO_2025", root / "metadata/raw/PRJNA1242368_assembly_2026-08-30.json"),
        ("CHINA_NINGXIA_2022", root / "metadata/raw/PRJNA745492_assembly_2026-08-30.json"),
    ]
    assemblies = [row for dataset_id, path in assembly_specs for row in ncbi_rows(path, dataset_id)]
    write_csv(phase0 / "ncbi_assembly_index.csv", assemblies, ASSEMBLY_HEADER)

    china_supplement = root / "literature/supplements/PMC9211431.1/aac.02188-21-s0001.pdf"
    china_mic = extract_china60_mic(china_supplement)
    write_csv(phase0 / "china_ningxia_60_mic.csv", china_mic, list(china_mic[0]))
    china_published_calls = extract_china60_published_calls(china_supplement)
    write_csv(phase0 / "china_ningxia_60_published_calls.csv", china_published_calls, list(china_published_calls[0]))
    china_published_by_id = {row["isolate_id"]: row for row in china_published_calls}

    hpgp_source = next((root / "literature/code_snapshots/hpgp_antimicrobialresistance_main_2026-08-30").rglob("Phenotypes_binary.tsv"))
    hpgp = read_csv(hpgp_source, delimiter="\t")
    hpgp_rows = [{"isolate_id": r["ID"], "clarithromycin": r["Clari"], "levofloxacin": r["Levo"]} for r in hpgp]
    write_csv(phase0 / "hpgp_binary_phenotype.csv", hpgp_rows, ["isolate_id", "clarithromycin", "levofloxacin"])

    zenodo_source = root / "metadata/extracted_supplements/PMC10878306.1/Table_1.xlsx/Table_1.csv"
    zenodo = read_csv(zenodo_source)
    zenodo_rows = []
    for r in zenodo:
        zenodo_rows.append(
            {
                "isolate_id": r["PID"],
                "age_years": r["Age"],
                "gender": r["Gender"],
                "clarithromycin": r["Clarithromycin"],
                "levofloxacin": r["Levofloxacin"],
                "amoxicillin": r["Amoxicillin"],
                "furazolidone": r["Furazolidone"],
                "tetracycline": r["Tetracycline"],
                "metronidazole": r["Metronidazole"],
            }
        )
    write_csv(phase0 / "zenodo_10369064_phenotype.csv", zenodo_rows, list(zenodo_rows[0]))

    assembly_by_dataset_strain = {(str(r["dataset_id"]), str(r["strain"]).upper()): r for r in assemblies if r["strain"]}
    for r in assemblies:
        strain = str(r["strain"]).upper()
        if r["dataset_id"] == "HPGP_GLOBAL" and strain.startswith("HPGP-"):
            assembly_by_dataset_strain[("HPGP_GLOBAL", strain.removeprefix("HPGP-"))] = r
    hpgp_mapped = sum(("HPGP_GLOBAL", r["isolate_id"].upper()) in assembly_by_dataset_strain for r in hpgp_rows)
    china_mapped = sum(("CHINA_NINGXIA_2022", r["isolate_id"].upper()) in assembly_by_dataset_strain for r in china_mic)
    zenodo_meta = json.loads((root / "metadata/raw/zenodo_10369064_2026-08-30.json").read_text(encoding="utf-8"))
    zenodo_files = {Path(item["key"]).stem.upper() for item in zenodo_meta["files"] if item["key"].lower().endswith(".zip")}
    zenodo_mapped = sum(r["isolate_id"].upper() in zenodo_files for r in zenodo_rows)
    mapping_rows = [
        {"dataset_id": "HPGP_GLOBAL", "phenotype_rows": len(hpgp_rows), "sequence_rows": sum(r["dataset_id"] == "HPGP_GLOBAL" for r in assemblies), "exact_id_links": hpgp_mapped, "mapping_status": "PASS" if hpgp_mapped == len(hpgp_rows) else "PARTIAL"},
        {"dataset_id": "CHINA_NINGXIA_2022", "phenotype_rows": len(china_mic), "sequence_rows": sum(r["dataset_id"] == "CHINA_NINGXIA_2022" for r in assemblies), "exact_id_links": china_mapped, "mapping_status": "PASS" if china_mapped == len(china_mic) else "PARTIAL"},
        {"dataset_id": "ZENODO_10369064", "phenotype_rows": len(zenodo_rows), "sequence_rows": len(zenodo_files), "exact_id_links": zenodo_mapped, "mapping_status": "PASS" if zenodo_mapped == len(zenodo_rows) else "PARTIAL"},
    ]
    write_csv(phase0 / "mapping_summary.csv", mapping_rows, ["dataset_id", "phenotype_rows", "sequence_rows", "exact_id_links", "mapping_status"])

    isolate_rows: list[dict[str, object]] = []
    for r in hpgp_rows:
        assembly = assembly_by_dataset_strain.get(("HPGP_GLOBAL", r["isolate_id"].upper()), {})
        isolate_rows.append({"dataset_id": "HPGP_GLOBAL", "isolate_id": r["isolate_id"], "patient_id": "", "sample_id": assembly.get("biosample", ""), "run_id": "", "assembly_id": assembly.get("genbank_accession", assembly.get("assembly_accession", "")), "country": r["isolate_id"].split("-")[0], "site": "", "collection_year": "", "clinical_diagnosis": "", "primary_or_post_treatment": "unknown", "lineage_reported": "available in HpGP metadata", "lineage_recomputed": "pending Phase 1", "near_clone_group": "pending Phase 1", "included": "yes" if assembly else "no", "exclusion_reason": "" if assembly else "no exact assembly-name link", "notes": "binary phenotype from frozen GitLab commit 4b2e5d00"})
    for r in china_mic:
        assembly = assembly_by_dataset_strain.get(("CHINA_NINGXIA_2022", r["isolate_id"].upper()), {})
        isolate_rows.append({"dataset_id": "CHINA_NINGXIA_2022", "isolate_id": r["isolate_id"], "patient_id": "", "sample_id": assembly.get("biosample", ""), "run_id": "", "assembly_id": assembly.get("genbank_accession", assembly.get("assembly_accession", "")), "country": "China", "site": "Ningxia", "collection_year": "", "clinical_diagnosis": "", "primary_or_post_treatment": "mixed; eradication history field retained", "lineage_reported": "reported in publication", "lineage_recomputed": "pending Phase 1", "near_clone_group": "pending smoke test", "included": "yes" if assembly else "no", "exclusion_reason": "" if assembly else "no exact assembly-name link", "notes": f"eradication_history={r['eradication_history']}"})
    for r in zenodo_rows:
        linked = r["isolate_id"].upper() in zenodo_files
        isolate_rows.append({"dataset_id": "ZENODO_10369064", "isolate_id": r["isolate_id"], "patient_id": r["isolate_id"], "sample_id": r["isolate_id"], "run_id": "", "assembly_id": "", "country": "China", "site": "multicentre", "collection_year": "", "clinical_diagnosis": "", "primary_or_post_treatment": "unknown", "lineage_reported": "study phylogeny only", "lineage_recomputed": "pending Phase 1", "near_clone_group": "pending smoke test", "included": "yes" if linked else "no", "exclusion_reason": "" if linked else "no exact Zenodo file-name link", "notes": "raw paired-end reads in per-isolate ZIP"})
    isolate_header = ["dataset_id", "isolate_id", "patient_id", "sample_id", "run_id", "assembly_id", "country", "site", "collection_year", "clinical_diagnosis", "primary_or_post_treatment", "lineage_reported", "lineage_recomputed", "near_clone_group", "included", "exclusion_reason", "notes"]
    write_csv(root / "metadata/isolate_manifest.csv", isolate_rows, isolate_header)

    phenotype_rows: list[dict[str, object]] = []
    for r in hpgp_rows:
        for drug in ("clarithromycin", "levofloxacin"):
            phenotype_rows.append({"dataset_id": "HPGP_GLOBAL", "isolate_id": r["isolate_id"], "antibiotic": drug, "mic_raw": "", "mic_numeric": "", "mic_operator": "", "unit": "mg/L", "ast_method": "centralized phenotyping; exact per-isolate MIC not in code snapshot", "medium": "see publication", "incubation": "see publication", "breakpoint_standard": "publication-defined", "breakpoint_version": "under audit", "susceptibility_original": "R" if r[drug] == "1" else "S", "susceptibility_recomputed": "", "borderline_mic": "unknown", "phenotype_quality": "A_binary_public_exact_link"})
    for r in china_mic:
        for drug in ("clarithromycin", "levofloxacin", "amoxicillin", "furazolidone", "tetracycline", "metronidazole"):
            numeric, operator = numeric_mic(r[drug])
            call, borderline = breakpoint_call(drug, r[drug])
            original = china_published_by_id[r["isolate_id"]][drug]
            agreement = original == call or (call == "I" and original == "S")
            phenotype_rows.append({"dataset_id": "CHINA_NINGXIA_2022", "isolate_id": r["isolate_id"], "antibiotic": drug, "mic_raw": r[drug], "mic_numeric": numeric, "mic_operator": operator, "unit": "mg/L", "ast_method": "Etest (BIO-KONT)", "medium": "Mueller-Hinton blood agar", "incubation": "72 h; microaerophilic; 37 C", "breakpoint_standard": "EUCAST" if drug != "furazolidone" else "study literature rule", "breakpoint_version": "publication says EUCAST but does not state version" if drug != "furazolidone" else "4 mg/L resistance breakpoint", "susceptibility_original": original, "susceptibility_recomputed": call, "borderline_mic": borderline, "phenotype_quality": "A_mic_public_exact_link" if agreement else "A_mic_public_exact_link_call_discordance"})
    for r in zenodo_rows:
        for drug in ("clarithromycin", "levofloxacin", "amoxicillin", "furazolidone", "tetracycline", "metronidazole"):
            raw = r[drug]
            normalized = {"Resistant": "R", "Sensitive": "S", "Intermediate": "I", "Lowly-sensitive": "I"}.get(raw, raw)
            phenotype_rows.append({"dataset_id": "ZENODO_10369064", "isolate_id": r["isolate_id"], "antibiotic": drug, "mic_raw": "", "mic_numeric": "", "mic_operator": "", "unit": "disk diameter mm", "ast_method": "Kirby-Bauer disk diffusion", "medium": "Columbia blood agar", "incubation": "3-4 days; microaerophilic; 37 C", "breakpoint_standard": "study-defined", "breakpoint_version": "Supplementary Figure 1 zone-diameter thresholds", "susceptibility_original": normalized, "susceptibility_recomputed": "", "borderline_mic": "unknown", "phenotype_quality": "B_categorical_public_exact_link_custom_threshold"})
    phenotype_header = ["dataset_id", "isolate_id", "antibiotic", "mic_raw", "mic_numeric", "mic_operator", "unit", "ast_method", "medium", "incubation", "breakpoint_standard", "breakpoint_version", "susceptibility_original", "susceptibility_recomputed", "borderline_mic", "phenotype_quality"]
    write_csv(root / "metadata/phenotype_manifest.csv", phenotype_rows, phenotype_header)

    balance_rows = []
    for dataset_id in sorted({str(r["dataset_id"]) for r in phenotype_rows}):
        for drug in ("clarithromycin", "levofloxacin"):
            values = [str(r["susceptibility_original"]) for r in phenotype_rows if r["dataset_id"] == dataset_id and r["antibiotic"] == drug]
            counts = Counter(values)
            balance_rows.append({"dataset_id": dataset_id, "antibiotic": drug, "n": len(values), "susceptible": counts["S"], "intermediate": counts["I"], "resistant": counts["R"], "other": len(values) - counts["S"] - counts["I"] - counts["R"]})
    write_csv(phase0 / "primary_drug_class_balance.csv", balance_rows, ["dataset_id", "antibiotic", "n", "susceptible", "intermediate", "resistant", "other"])

    china_comparison = []
    for drug in ("clarithromycin", "levofloxacin", "amoxicillin", "furazolidone", "tetracycline", "metronidazole"):
        subset = [row for row in phenotype_rows if row["dataset_id"] == "CHINA_NINGXIA_2022" and row["antibiotic"] == drug]
        exact = sum(row["susceptibility_original"] == row["susceptibility_recomputed"] for row in subset)
        intermediate_as_s = sum(row["susceptibility_original"] == "S" and row["susceptibility_recomputed"] == "I" for row in subset)
        china_comparison.append({"antibiotic": drug, "n": len(subset), "exact_agreement": exact, "intermediate_as_s_in_publication": intermediate_as_s, "other_discordance": len(subset) - exact - intermediate_as_s})
    write_csv(phase0 / "china_ningxia_60_call_comparison.csv", china_comparison, ["antibiotic", "n", "exact_agreement", "intermediate_as_s_in_publication", "other_discordance"])


if __name__ == "__main__":
    main()
