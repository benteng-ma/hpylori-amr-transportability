from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_ningxia_published_genotypes import parse_table_s9, qrdr_substitutions


def test_table_s9_parser_requires_complete_unique_table():
    text = "\n".join(
        f"SHZY{i:02d} " + " ".join(["R", "S"] * 7)
        for i in range(1, 61)
    )
    rows = parse_table_s9(text)
    assert len(rows) == 60
    assert rows[0]["clarithromycin_published_genotype"] == "R"
    assert rows[0]["clarithromycin_published_phenotype"] == "S"
    assert rows[-1]["isolate_id"] == "SHZY60"


def test_table_s9_parser_fails_closed_on_partial_table():
    with pytest.raises(ValueError, match="expected 60"):
        parse_table_s9("SHZY01 " + " ".join(["R", "S"] * 7))


def test_qrdr_translation_uses_frozen_amino_acid_positions():
    codons = ["AAA"] * 100
    codons[86] = "TAT"  # N87Y
    codons[87] = "GCA"  # A88A
    codons[90] = "GAT"  # D91D
    changes, prediction = qrdr_substitutions("".join(codons), gene_start_1_based=1)
    assert changes == "N87Y"
    assert prediction == "S"
