from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parents[1]


EXPECTED = {
    "metadata/dataset_manifest.csv": [
        "dataset_id", "publication", "pmid", "doi", "repository", "accession",
        "country", "study_period", "n_reported", "n_genomes_available",
        "n_isolates_with_phenotype", "antibiotics", "mic_available", "binary_only",
        "phenotype_method", "breakpoint_standard", "breakpoint_version",
        "raw_reads_available", "assembly_available", "lineage_metadata_available",
        "patient_id_available", "published_model_source", "primary_external_candidate",
        "verified", "limitations",
    ],
    "metadata/isolate_manifest.csv": [
        "dataset_id", "isolate_id", "patient_id", "sample_id", "run_id", "assembly_id",
        "country", "site", "collection_year", "clinical_diagnosis",
        "primary_or_post_treatment", "lineage_reported", "lineage_recomputed",
        "near_clone_group", "included", "exclusion_reason", "notes",
    ],
    "metadata/phenotype_manifest.csv": [
        "dataset_id", "isolate_id", "antibiotic", "mic_raw", "mic_numeric",
        "mic_operator", "unit", "ast_method", "medium", "incubation",
        "breakpoint_standard", "breakpoint_version", "susceptibility_original",
        "susceptibility_recomputed", "borderline_mic", "phenotype_quality",
    ],
}


def test_required_manifest_headers_are_frozen():
    for relative, expected in EXPECTED.items():
        with (ROOT / relative).open(newline="", encoding="utf-8") as handle:
            observed = next(csv.reader(handle))
        assert observed == expected


def test_canonical_root_is_e_drive():
    text = (ROOT / "config/project.yaml").read_text(encoding="utf-8")
    assert 'canonical_root: "E:/1 Codex project/hpylori-amr-transportability"' in text
