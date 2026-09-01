import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_rows(relative):
    with (ROOT / relative).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_phase0_decision_and_mapping_counts_are_frozen():
    decision = json.loads((ROOT / "reports/phase0_decision.json").read_text(encoding="utf-8"))
    assert decision["decision"] == "CONDITIONAL_GO_CLR_LVX_BINARY"
    mapping = {row["dataset_id"]: row for row in read_rows("metadata/phase0/mapping_summary.csv")}
    assert int(mapping["HPGP_GLOBAL"]["exact_id_links"]) == 414
    assert int(mapping["CHINA_NINGXIA_2022"]["exact_id_links"]) == 60
    assert int(mapping["ZENODO_10369064"]["exact_id_links"]) == 52


def test_ningxia_public_calls_are_reconciled():
    comparison = {row["antibiotic"]: row for row in read_rows("metadata/phase0/china_ningxia_60_call_comparison.csv")}
    assert int(comparison["clarithromycin"]["intermediate_as_s_in_publication"]) == 6
    for drug, row in comparison.items():
        assert int(row["other_discordance"]) == 0, drug


def test_isolate_and_phenotype_keys_are_unique():
    isolates = read_rows("metadata/isolate_manifest.csv")
    isolate_keys = [(row["dataset_id"], row["isolate_id"]) for row in isolates]
    assert len(isolate_keys) == len(set(isolate_keys))
    phenotypes = read_rows("metadata/phenotype_manifest.csv")
    phenotype_keys = [(row["dataset_id"], row["isolate_id"], row["antibiotic"]) for row in phenotypes]
    assert len(phenotype_keys) == len(set(phenotype_keys))


def test_reference_numbering_and_23s_copies_are_explicit():
    variants = read_rows("results/smoke/known_variant_calls.csv")
    reference_name = "GCF_000008525.1_ASM852v1_genomic.fna"
    reference = [row for row in variants if row["assembly"] == reference_name]
    rrna_copies = {row["copy"] for row in reference if row["gene"] == "23S_rRNA"}
    assert rrna_copies == {"1", "2"}
    gyr_a = {int(row["position"]): row["reference"] for row in reference if row["gene"] == "gyrA"}
    assert gyr_a == {87: "N", 88: "A", 91: "D"}


def test_frozen_panels_are_not_refitted():
    panels = read_rows("metadata/published_panel_manifest.csv")
    eligible = [row for row in panels if row["phase0_eligibility"] == "eligible_frozen_baseline"]
    assert {row["panel_id"] for row in eligible} == {"HPGP_CLR_PANEL", "HPGP_LEV_PANEL"}
    assert all(row["reproducible"] == "yes" for row in eligible)


def test_smoke_outputs_are_complete_but_not_overinterpreted():
    relatedness = read_rows("results/smoke/pairwise_relatedness.csv")
    assert len(relatedness) == 36
    assert not any(row["near_clone_phase0"] == "yes" for row in relatedness)
    variants = read_rows("results/smoke/known_variant_calls.csv")
    assert any(row["status"].startswith("UNCALLABLE") for row in variants)
    assert any(row["known_resistance_marker"] == "yes" for row in variants)
    pileup = read_rows("results/smoke/S27_marker_pileup.csv")
    assert len(pileup) == 7
    assert all(row["status"] in {"PASS", "LOW_COVERAGE"} for row in pileup)
    assert any(row["status"] == "LOW_COVERAGE" for row in pileup)
