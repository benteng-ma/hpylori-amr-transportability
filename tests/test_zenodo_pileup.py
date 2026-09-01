from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from process_zenodo_reads import parse_pileup_bases


def test_parse_pileup_bases_counts_reference_and_variants():
    counts = parse_pileup_bases(".,..Gg^F.$+2AA,-1t", "A")
    assert counts == {"A": 6, "C": 0, "G": 2, "T": 0}
