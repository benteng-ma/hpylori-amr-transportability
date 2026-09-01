import json
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_COFIRST = ["Benteng Ma", "Bing Chen"]
EXPECTED_ORDER = ["Benteng Ma", "Bing Chen", "Ting Cai", "Xiao-ming Liu", "Fen Wang"]
EXPECTED_CORRESPONDING = ["Xiao-ming Liu", "Fen Wang"]
TING_ONLY_AFFILIATIONS = [5, 6, 7]
TING_AFFILIATION_TEXT = [
    "Hunan Provincial University Key Laboratory of the Fundamental and Clinical Research on Neurodegenerative Diseases",
    "Hunan Provincial University Key Laboratory of the Fundamental and Clinical Research on Functional Nucleic Acid",
    "Hunan Provincial Key Laboratory of the Research and Development of Novel Pharmaceutical Preparations",
]
FUNDING_TEXT = [
    '"The 14th Five-Year Plan" Application Characteristic Discipline of Hunan Province (Clinical Medicine)',
    "Aid Program for Science and Technology Innovative Research Team in Higher Educational Institutions of Hunan Province, China",
]
REMOVED_AUTHORS = [
    "Yin Zhao",
    "Shuo-Yi Yao",
    "Xin-meng Li",
    "Minglin Zhang",
    "Aojian Deng",
    "Qian Yin",
    "Chunbing Zheng",
]


def docx_part_text(path: Path, part: str) -> str:
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read(part))
    return " ".join(text.strip() for text in root.itertext() if text.strip())


def xlsx_text(path: Path) -> str:
    with zipfile.ZipFile(path) as archive:
        members = [name for name in archive.namelist() if name.startswith("xl/") and name.endswith(".xml")]
        return " ".join(archive.read(name).decode("utf-8", errors="ignore") for name in members)


def test_authorship_registry_locks_exact_co_first_order() -> None:
    registry = json.loads((ROOT / "metadata/authorship.json").read_text(encoding="utf-8"))
    authors = sorted(registry["authors"], key=lambda item: item["order"])
    assert [item["name"] for item in authors] == EXPECTED_ORDER
    assert [item["name"] for item in authors if item["co_first"]] == EXPECTED_COFIRST
    assert registry["co_first_author_names"] == EXPECTED_COFIRST
    assert registry["other_first_authors"] == []
    assert [item["name"] for item in authors if item["corresponding"]] == EXPECTED_CORRESPONDING
    assert registry["corresponding_author_names"] == EXPECTED_CORRESPONDING
    ting = next(item for item in authors if item["name"] == "Ting Cai")
    assert ting["affiliation_ids"] == TING_ONLY_AFFILIATIONS
    assert all(
        not set(item["affiliation_ids"]) & set(TING_ONLY_AFFILIATIONS)
        for item in authors
        if item["name"] != "Ting Cai"
    )
    assert [registry["affiliations"][str(key)].split(", Changsha Medical University", 1)[0] for key in TING_ONLY_AFFILIATIONS] == TING_AFFILIATION_TEXT


def test_versioned_sources_preserve_locked_order_and_equal_contribution() -> None:
    main = (ROOT / "manuscript/manuscript.qmd").read_text(encoding="utf-8")
    supplement = (ROOT / "manuscript/supplementary_methods.qmd").read_text(encoding="utf-8")
    cover = (ROOT / "submission/cover_letter_Microbial_Genomics.md").read_text(encoding="utf-8")
    title_page = (ROOT / "submission/title_page_template.md").read_text(encoding="utf-8")
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    for text in (main, supplement, cover, title_page):
        positions = [text.index(name) for name in EXPECTED_ORDER]
        assert positions == sorted(positions)
        assert all(name not in text for name in REMOVED_AUTHORS)
    assert all(name not in citation for name in REMOVED_AUTHORS)
    author_string = "; ".join(EXPECTED_ORDER)
    assert f'author: "{author_string}"' in main
    assert f'author: "{author_string}"' in supplement
    assert "contributed equally" in main
    assert all(value in main for value in FUNDING_TEXT)
    assert all(value in title_page for value in TING_AFFILIATION_TEXT)
    assert "No other author is designated as a first author." in cover
    assert "Xiao-ming Liu" in cover and "Fen Wang" in cover


def test_submission_docx_files_preserve_locked_authorship() -> None:
    files = [
        ROOT / "submission/Audit_first_Hpylori_AMR_manuscript.docx",
        ROOT / "submission/Audit_first_Hpylori_AMR_supplement.docx",
        ROOT / "submission/Audit_first_Hpylori_AMR_cover_letter.docx",
    ]
    for path in files:
        body = docx_part_text(path, "word/document.xml")
        core = docx_part_text(path, "docProps/core.xml")
        positions = [body.index(name) for name in EXPECTED_ORDER]
        assert positions == sorted(positions)
        assert "; ".join(EXPECTED_ORDER) in core
        assert all(name not in body and name not in core for name in REMOVED_AUTHORS)
        assert "Author names and affiliations required" not in body
    assert "contributed equally" in docx_part_text(files[0], "word/document.xml")
    main_body = docx_part_text(files[0], "word/document.xml")
    supplement_body = docx_part_text(files[1], "word/document.xml")
    assert all(value in main_body for value in FUNDING_TEXT)
    assert all(value in main_body and value in supplement_body for value in TING_AFFILIATION_TEXT)
    assert "No other author is designated as a first author." in docx_part_text(
        files[2], "word/document.xml"
    )


def test_supplementary_workbook_preserves_five_author_lock() -> None:
    text = xlsx_text(ROOT / "submission/Supplementary_Data_S1-S20.xlsx")
    positions = [text.index(name) for name in EXPECTED_ORDER]
    assert positions == sorted(positions)
    for token in ("#REF!", "#DIV/0!", "#VALUE!", "#NAME?", "#N/A"):
        assert token not in text
    assert all(name not in text for name in REMOVED_AUTHORS)
