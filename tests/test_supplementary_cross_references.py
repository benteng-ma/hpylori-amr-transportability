import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAIN_TEXT = (ROOT / "manuscript/manuscript.qmd").read_text(encoding="utf-8")


def cited_after_label(label: str, number: int) -> bool:
    pattern = rf"{label}(?:s)? [^.;\n]*S{number}(?!\d)"
    return re.search(pattern, MAIN_TEXT) is not None


def test_every_supplementary_figure_is_cited_in_the_main_text() -> None:
    assert [number for number in range(1, 9) if not cited_after_label("Supplementary Figure", number)] == []


def test_every_reader_facing_supplementary_table_group_is_cited_in_the_main_text() -> None:
    reader_facing = [*range(1, 12), *range(13, 21)]
    assert [number for number in reader_facing if not cited_after_label("Table", number)] == []


def test_worksheet_s12_is_not_presented_as_a_reader_facing_result_table() -> None:
    generated = (ROOT / "manuscript/supplementary_tables_generated.md").read_text(encoding="utf-8")
    supplement = (ROOT / "manuscript/supplementary_methods.qmd").read_text(encoding="utf-8")
    assert "Table S12. Main-figure source-data manifest" not in generated
    assert "Worksheet S12" in supplement
