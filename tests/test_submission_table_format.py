from pathlib import Path
import zipfile
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "submission/Audit_first_Hpylori_AMR_manuscript.docx"
SUPPLEMENT = ROOT / "submission/Audit_first_Hpylori_AMR_supplement.docx"
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def read_xml(archive: zipfile.ZipFile, member: str) -> ET.Element:
    return ET.fromstring(archive.read(member))


def test_article_files_have_no_running_headers_or_running_title() -> None:
    for path in (MAIN, SUPPLEMENT):
        with zipfile.ZipFile(path) as archive:
            document = read_xml(archive, "word/document.xml")
            assert "Running title:" not in " ".join(document.itertext())
            section_properties = document.findall(f".//{W}sectPr")
            assert section_properties
            assert all(section.find(f"{W}headerReference") is None for section in section_properties)
            assert not any(name.startswith("word/header") for name in archive.namelist())


def test_main_table_uses_strict_three_line_rules() -> None:
    with zipfile.ZipFile(MAIN) as archive:
        document = read_xml(archive, "word/document.xml")
    tables = document.findall(f".//{W}tbl")
    assert len(tables) == 1
    table = tables[0]
    table_borders = table.find(f"{W}tblPr/{W}tblBorders")
    assert table_borders is not None
    assert table_borders.find(f"{W}top").get(f"{W}val") == "single"
    assert table_borders.find(f"{W}bottom").get(f"{W}val") == "single"
    for edge in ("left", "right", "insideH", "insideV"):
        assert table_borders.find(f"{W}{edge}").get(f"{W}val") == "nil"

    rows = table.findall(f"{W}tr")
    assert rows
    for cell in rows[0].findall(f"{W}tc"):
        borders = cell.find(f"{W}tcPr/{W}tcBorders")
        assert borders is not None
        assert borders.find(f"{W}top").get(f"{W}val") == "single"
        assert borders.find(f"{W}bottom").get(f"{W}val") == "single"
        assert cell.find(f"{W}tcPr/{W}shd") is None
    for row_index, row in enumerate(rows[1:], start=1):
        for cell in row.findall(f"{W}tc"):
            borders = cell.find(f"{W}tcPr/{W}tcBorders")
            assert borders is not None
            assert all(
                borders.find(f"{W}{edge}").get(f"{W}val") == "nil"
                for edge in ("top", "left", "right", "insideH", "insideV")
            )
            expected_bottom = "single" if row_index == len(rows) - 1 else "nil"
            assert borders.find(f"{W}bottom").get(f"{W}val") == expected_bottom


def test_supplement_prints_reader_facing_groups_as_strict_three_line_tables() -> None:
    with zipfile.ZipFile(SUPPLEMENT) as archive:
        document = read_xml(archive, "word/document.xml")
    text = " ".join(document.itertext())
    for number in [*range(1, 12), *range(13, 21)]:
        assert f"Table S{number}" in text
    assert "Table S12. Main-figure source-data manifest" not in text
    assert "Worksheet S12" in text

    tables = document.findall(f".//{W}tbl")
    assert len(tables) >= 35
    for table in tables:
        table_borders = table.find(f"{W}tblPr/{W}tblBorders")
        assert table_borders is not None
        assert table_borders.find(f"{W}top").get(f"{W}val") == "single"
        assert table_borders.find(f"{W}bottom").get(f"{W}val") == "single"
        for edge in ("left", "right", "insideH", "insideV"):
            assert table_borders.find(f"{W}{edge}").get(f"{W}val") == "nil"

        rows = table.findall(f"{W}tr")
        assert rows
        for cell in rows[0].findall(f"{W}tc"):
            borders = cell.find(f"{W}tcPr/{W}tcBorders")
            assert borders is not None
            assert borders.find(f"{W}top").get(f"{W}val") == "single"
            assert borders.find(f"{W}bottom").get(f"{W}val") == "single"
        for row_index, row in enumerate(rows[1:], start=1):
            for cell in row.findall(f"{W}tc"):
                borders = cell.find(f"{W}tcPr/{W}tcBorders")
                assert borders is not None
                expected_bottom = "single" if row_index == len(rows) - 1 else "nil"
                assert borders.find(f"{W}bottom").get(f"{W}val") == expected_bottom
