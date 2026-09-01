from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from process_phase2_assemblies import audit_existing, fasta_metrics


def test_fasta_metrics_on_reference():
    path = ROOT / "data/interim/reference/ncbi_dataset/data/GCF_000008525.1/GCF_000008525.1_ASM852v1_genomic.fna"
    metrics = fasta_metrics(path)
    assert 1_600_000 < metrics["assembly_size_bp"] < 1_700_000
    assert metrics["contigs"] >= 1
    assert 35 < metrics["gc_percent"] < 45


def test_existing_assembly_keeps_dataset_and_isolate_provenance(tmp_path):
    dataset = tmp_path / "ZENODO_10369064"
    dataset.mkdir()
    assembly = dataset / "S01.fna"
    assembly.write_text(">contig\n" + "ACGT" * 375000 + "\n", encoding="utf-8")
    result = audit_existing(assembly, tmp_path)
    assert result["dataset_id"] == "ZENODO_10369064"
    assert result["isolate_id"] == "S01"
    assert result["size_gate"] == "PASS"
