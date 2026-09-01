from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from cluster_lineages import variant_matrix


def test_variant_matrix_uses_reference_and_excludes_gaps(tmp_path):
    alignment = tmp_path / "toy.aln"
    alignment.write_text(">Reference\nAACCGG\n>A__1\nAATCGG\n>B__2\nAAC-GA\n", encoding="utf-8")
    names, matrix, variants = variant_matrix(alignment)
    assert names == ["A__1", "B__2"]
    assert variants == 2
    assert matrix.nnz == 2


def test_variant_matrix_feature_space_is_training_defined(tmp_path):
    alignment = tmp_path / "toy.aln"
    alignment.write_text(">Reference\nAACCGG\n>DEV__1\nAATCGG\n>DEV__2\nAACCGG\n>EXT__2\nAAC-GA\n", encoding="utf-8")
    names, matrix, variants = variant_matrix(alignment, "DEV")
    assert names == ["DEV__1", "DEV__2", "EXT__2"]
    assert variants == 1
    assert matrix.nnz == 1


def test_variant_matrix_keeps_distinct_nonreference_alleles(tmp_path):
    alignment = tmp_path / "toy.aln"
    alignment.write_text(">Reference\nAAA\n>DEV__1\nCAA\n>DEV__2\nGAA\n>DEV__3\nAAA\n", encoding="utf-8")
    _names, matrix, variants = variant_matrix(alignment, "DEV")
    assert variants == 1
    assert matrix.shape == (3, 2)
    assert matrix.nnz == 2


def test_variant_matrix_locates_reference_independent_of_record_order(tmp_path):
    alignment = tmp_path / "toy.aln"
    alignment.write_text(">DEV__1\nAATCGG\n>Reference\nAACCGG\n>EXT__2\nAAC-GA\n", encoding="utf-8")
    names, matrix, variants = variant_matrix(alignment)
    assert names == ["DEV__1", "EXT__2"]
    assert variants == 2
    assert matrix.shape == (2, 2)
    assert matrix.nnz == 2
