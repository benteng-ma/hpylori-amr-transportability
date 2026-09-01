#!/usr/bin/env python3
"""Audit Ningxia source-reported genotype calls against deposited assemblies.

The publication supplement reports per-isolate genotypic resistance (G-R) and
phenotypic resistance (P-R) in Table S9. These calls are retained as a source
audit only and never substituted for independently recalled frozen-panel calls.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from benchmark_frozen_panels import exact_interval, point_metrics, read_csv, write_csv


DRUG_COLUMNS = {
    "clarithromycin": (0, 1),
    "levofloxacin": (2, 3),
}
SOURCE_RULES = {
    "clarithromycin": "source Table S2/full text: 23S rRNA A2143G",
    "levofloxacin": "source full text: gyrA N87K/I/Y or D91N/G",
}
GYRA_START_1_BASED = 752512
GYRA_REFERENCE_AA = {87: "N", 88: "A", 91: "D"}
GYRA_FROZEN_MARKERS = {"N87K", "N87I", "A88V", "A88P", "D91G", "D91N", "D91Y"}


def parse_table_s9(text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for raw_line in text.splitlines():
        tokens = raw_line.strip().split()
        if len(tokens) != 15 or not tokens[0].startswith("SHZY"):
            continue
        if not all(value in {"S", "R"} for value in tokens[1:]):
            continue
        row = {"isolate_id": tokens[0]}
        values = tokens[1:]
        for drug, (genotype_index, phenotype_index) in DRUG_COLUMNS.items():
            row[f"{drug}_published_genotype"] = values[genotype_index]
            row[f"{drug}_published_phenotype"] = values[phenotype_index]
        rows.append(row)
    rows.sort(key=lambda row: int(row["isolate_id"].removeprefix("SHZY")))
    if len(rows) != 60 or len({row["isolate_id"] for row in rows}) != 60:
        raise ValueError(f"expected 60 unique Table S9 isolates, found {len(rows)}")
    return rows


def recalled_substitutions(calls: list[dict[str, str]]) -> dict[tuple[str, str], str]:
    grouped: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in calls:
        if row["dataset_id"] != "CHINA_NINGXIA_2022" or row["status"] != "PASS":
            continue
        if row["reference"] == row["observed"]:
            continue
        if row["gene"] == "23S_rRNA" and row["position"] in {"2142", "2143"}:
            grouped[(row["isolate_id"], "clarithromycin")].add(row["change"])
        if row["gene"] == "gyrA" and row["position"] in {"87", "88", "91"}:
            grouped[(row["isolate_id"], "levofloxacin")].add(row["change"])
    return {key: ";".join(sorted(values)) for key, values in grouped.items()}


def qrdr_substitutions(sequence: str, gene_start_1_based: int = GYRA_START_1_BASED) -> tuple[str, str]:
    changes: list[str] = []
    for position, reference_aa in GYRA_REFERENCE_AA.items():
        codon_start = gene_start_1_based - 1 + (position - 1) * 3
        codon = sequence[codon_start:codon_start + 3].upper()
        if len(codon) != 3 or any(base not in {"A", "C", "G", "T"} for base in codon):
            return "", "UNCALLABLE"
        observed_aa = str(Seq(codon).translate())
        if observed_aa != reference_aa:
            changes.append(f"{reference_aa}{position}{observed_aa}")
    change_string = ";".join(changes)
    prediction = "R" if set(changes) & GYRA_FROZEN_MARKERS else "S"
    return change_string, prediction


def snippy_qrdr_calls(root: Path) -> dict[str, tuple[str, str]]:
    output: dict[str, tuple[str, str]] = {}
    for directory in sorted((root / "data/phylogeny/snippy").glob("CHINA_NINGXIA_2022__SHZY*")):
        alignment = directory / "snps.aligned.fa"
        if not alignment.exists():
            continue
        records = list(SeqIO.parse(alignment, "fasta"))
        if len(records) != 1:
            raise ValueError(f"expected one Snippy aligned record in {alignment}, found {len(records)}")
        output[directory.name.split("__", 1)[1]] = qrdr_substitutions(str(records[0].seq))
    return output


def performance_row(
    drug: str,
    estimate_source: str,
    scope: str,
    rows: list[dict[str, str]],
    prediction_field: str,
    interpretation: str,
) -> dict[str, object]:
    truth = [1 if row["phenotype"] == "R" else 0 for row in rows]
    prediction = [1 if row[prediction_field] == "R" else 0 for row in rows]
    metrics = point_metrics(truth, prediction)
    sensitivity_low, sensitivity_high = exact_interval(int(metrics["tp"]), int(metrics["n_resistant"]))
    specificity_low, specificity_high = exact_interval(int(metrics["tn"]), int(metrics["n_susceptible"]))
    return {
        "row_type": "PHENOTYPE_PERFORMANCE",
        "antibiotic": drug,
        "estimate_source": estimate_source,
        "scope": scope,
        **metrics,
        "sensitivity_ci_low": sensitivity_low,
        "sensitivity_ci_high": sensitivity_high,
        "specificity_ci_low": specificity_low,
        "specificity_ci_high": specificity_high,
        "n_agree": "",
        "call_concordance": "",
        "source_rule": SOURCE_RULES[drug] if estimate_source == "SOURCE_REPORTED_G_R" else "frozen HpGP panel",
        "interpretation": interpretation,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    root = args.root.resolve()

    source_text = (root / "metadata/extracted_supplements/PMC9211431.1/aac.02188-21-s0001.pdf.txt").read_text(encoding="utf-8")
    published = parse_table_s9(source_text)
    phenotype = {
        (row["isolate_id"], row["antibiotic"]): row["susceptibility_original"]
        for row in read_csv(root / "metadata/phenotype_manifest.csv")
        if row["dataset_id"] == "CHINA_NINGXIA_2022" and row["antibiotic"] in DRUG_COLUMNS
    }
    qc = {
        row["isolate_id"]: row["final_qc_status"]
        for row in read_csv(root / "results/qc/assembly_qc_with_checkm2.csv")
        if row["dataset_id"] == "CHINA_NINGXIA_2022"
    }
    recalled = {
        (row["isolate_id"], row["antibiotic"]): row
        for row in read_csv(root / "results/external_validation/sample_level_predictions.csv")
        if row["dataset_id"] == "CHINA_NINGXIA_2022" and row["antibiotic"] in DRUG_COLUMNS
    }
    substitutions = recalled_substitutions(read_csv(root / "results/panels/assembly_marker_calls.csv"))
    snippy_calls = snippy_qrdr_calls(root)

    samples: list[dict[str, str]] = []
    for source in published:
        isolate = source["isolate_id"]
        for drug in DRUG_COLUMNS:
            source_phenotype = source[f"{drug}_published_phenotype"]
            manifest_phenotype = phenotype[(isolate, drug)]
            if source_phenotype != manifest_phenotype:
                raise ValueError(f"phenotype mismatch for {isolate} {drug}: {source_phenotype} != {manifest_phenotype}")
            current = recalled[(isolate, drug)]
            source_call = source[f"{drug}_published_genotype"]
            current_call = current["prediction"]
            both_callable = source_call in {"S", "R"} and current_call in {"S", "R"}
            snippy_substitutions, snippy_call = snippy_calls.get(isolate, ("", "")) if drug == "levofloxacin" else ("", "")
            samples.append({
                "dataset_id": "CHINA_NINGXIA_2022",
                "isolate_id": isolate,
                "antibiotic": drug,
                "phenotype": manifest_phenotype,
                "published_genotype_call": source_call,
                "recalled_frozen_panel_call": current_call,
                "recalled_target_substitutions": substitutions.get((isolate, drug), ""),
                "snippy_qrdr_substitutions": snippy_substitutions,
                "snippy_frozen_panel_call": snippy_call,
                "independent_recallers_agree": (
                    "yes" if current_call in {"S", "R"} and snippy_call in {"S", "R"} and current_call == snippy_call
                    else "no" if current_call in {"S", "R"} and snippy_call in {"S", "R"}
                    else ""
                ),
                "final_qc_status": qc[isolate],
                "both_calls_available": "yes" if both_callable else "no",
                "calls_concordant": "yes" if both_callable and source_call == current_call else "no" if both_callable else "",
                "source_rule": SOURCE_RULES[drug],
                "interpretation": "source-reported G-R is an audit comparator, not an independently recalled prediction",
            })

    summary: list[dict[str, object]] = []
    for drug in DRUG_COLUMNS:
        drug_rows = [row for row in samples if row["antibiotic"] == drug]
        for scope, eligible in (
            ("ALL_60_SOURCE_ROWS", drug_rows),
            ("FINAL_GENOME_QC_PASS", [row for row in drug_rows if row["final_qc_status"] == "PASS"]),
        ):
            summary.append(performance_row(
                drug,
                "SOURCE_REPORTED_G_R",
                scope,
                eligible,
                "published_genotype_call",
                "published per-isolate genotype audit against the independently transcribed Table S9 phenotype",
            ))
        recalled_eligible = [
            row for row in drug_rows
            if row["final_qc_status"] == "PASS" and row["recalled_frozen_panel_call"] in {"S", "R"}
        ]
        summary.append(performance_row(
            drug,
            "INDEPENDENT_RECALL",
            "FINAL_GENOME_QC_PASS_AND_TARGET_CALLABLE",
            recalled_eligible,
            "recalled_frozen_panel_call",
            "primary independently recalled frozen-panel estimate",
        ))
        comparable = [row for row in drug_rows if row["final_qc_status"] == "PASS" and row["both_calls_available"] == "yes"]
        agreements = sum(row["calls_concordant"] == "yes" for row in comparable)
        summary.append({
            "row_type": "CALL_REPRODUCIBILITY",
            "antibiotic": drug,
            "estimate_source": "SOURCE_REPORTED_G_R_VS_INDEPENDENT_RECALL",
            "scope": "FINAL_GENOME_QC_PASS_AND_BOTH_CALLABLE",
            "n": len(comparable),
            "n_resistant": "", "n_susceptible": "", "tp": "", "tn": "", "fp": "", "fn": "",
            "prevalence": "", "sensitivity": "", "specificity": "", "false_susceptible_rate": "",
            "false_resistant_rate": "", "balanced_accuracy": "", "mcc": "", "ppv": "", "npv": "",
            "auroc_binary_score": "", "auprc_binary_score": "", "brier_binary_score": "",
            "sensitivity_ci_low": "", "sensitivity_ci_high": "", "specificity_ci_low": "", "specificity_ci_high": "",
            "n_agree": agreements,
            "call_concordance": agreements / len(comparable) if comparable else "",
            "source_rule": SOURCE_RULES[drug],
            "interpretation": "descriptive reproducibility audit; discordance is not resolved by substituting source calls",
        })

    snippy_comparable = [
        row for row in samples
        if row["antibiotic"] == "levofloxacin"
        and row["final_qc_status"] == "PASS"
        and row["recalled_frozen_panel_call"] in {"S", "R"}
        and row["snippy_frozen_panel_call"] in {"S", "R"}
    ]
    snippy_agreements = sum(row["independent_recallers_agree"] == "yes" for row in snippy_comparable)
    summary.append({
        "row_type": "INDEPENDENT_CALLER_REPRODUCIBILITY",
        "antibiotic": "levofloxacin",
        "estimate_source": "TARGET_BLAST_VS_SNIPPY_WHOLE_GENOME_ALIGNMENT",
        "scope": "FINAL_GENOME_QC_PASS_AND_BOTH_CALLABLE",
        "n": len(snippy_comparable),
        "n_resistant": "", "n_susceptible": "", "tp": "", "tn": "", "fp": "", "fn": "",
        "prevalence": "", "sensitivity": "", "specificity": "", "false_susceptible_rate": "",
        "false_resistant_rate": "", "balanced_accuracy": "", "mcc": "", "ppv": "", "npv": "",
        "auroc_binary_score": "", "auprc_binary_score": "", "brier_binary_score": "",
        "sensitivity_ci_low": "", "sensitivity_ci_high": "", "specificity_ci_low": "", "specificity_ci_high": "",
        "n_agree": snippy_agreements,
        "call_concordance": snippy_agreements / len(snippy_comparable) if snippy_comparable else "",
        "source_rule": "frozen HpGP gyrA panel",
        "interpretation": "independent targeted-BLAST and whole-genome-alignment callers",
    })

    write_csv(root / "results/external_validation/ningxia_published_genotype_audit_samples.csv", samples)
    write_csv(root / "results/external_validation/ningxia_published_genotype_audit_summary.csv", summary)


if __name__ == "__main__":
    main()
