import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = (ROOT / "manuscript" / "manuscript.qmd").read_text(encoding="utf-8")
COVER = (ROOT / "submission" / "cover_letter_Microbial_Genomics.md").read_text(encoding="utf-8")
CHECKLIST = (ROOT / "submission" / "MGEN_submission_checklist.md").read_text(encoding="utf-8")
PORTAL = (ROOT / "submission" / "MGEN_portal_metadata.md").read_text(encoding="utf-8")
TITLE_PAGE = (ROOT / "submission" / "title_page_template.md").read_text(encoding="utf-8")
CONFIRMATION = (ROOT / "submission" / "MGEN_author_confirmation_form.md").read_text(encoding="utf-8")
CITATION = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
EXPECTED_TITLE = "External validation of genomic resistance prediction in Helicobacter pylori across independent cohorts"


def section(text: str, heading: str, next_heading: str | None = None) -> str:
    start = text.index(f"## {heading}") + len(f"## {heading}")
    if next_heading:
        end = text.index(f"## {next_heading}", start)
    else:
        match = re.search(r"\n## ", text[start:])
        end = start + match.start() if match else len(text)
    return text[start:end].strip()


def word_count(text: str) -> int:
    clean = re.sub(r"\[[^\]]+\]\([^\)]+\)", " ", text)
    clean = re.sub(r"[`*_#{}]", "", clean)
    return len(re.findall(r"\b[\w'-]+\b", clean))


def test_article_title_is_plain_and_synchronized():
    title = re.search(r'^title: "(.+)"$', MANUSCRIPT, flags=re.MULTILINE).group(1)
    assert title == EXPECTED_TITLE
    assert all(mark not in title for mark in ("-", "–", "—", ":"))
    assert "Audit-first" not in title
    assert EXPECTED_TITLE in COVER.replace("*", "")
    assert EXPECTED_TITLE in PORTAL
    assert EXPECTED_TITLE in TITLE_PAGE.replace("*", "")
    assert EXPECTED_TITLE in CONFIRMATION.replace("*", "")
    assert f'title: "{EXPECTED_TITLE}"' in CITATION


def test_required_front_matter_order_and_lengths():
    headings = ["## Abstract", "## Impact Statement", "## Data Summary", "## Introduction"]
    positions = [MANUSCRIPT.index(item) for item in headings]
    assert positions == sorted(positions)
    abstract = section(MANUSCRIPT, "Abstract", "Impact Statement")
    impact = section(MANUSCRIPT, "Impact Statement", "Data Summary")
    assert 150 <= word_count(abstract) <= 350
    assert word_count(impact) <= 200
    keyword_line = re.search(r"\*\*Keywords:\*\* (.+)", abstract).group(1)
    assert 3 <= len([item for item in keyword_line.split(";") if item.strip()]) <= 6


def test_main_text_is_within_research_article_indicative_range():
    body = MANUSCRIPT[MANUSCRIPT.index("## Introduction") : MANUSCRIPT.index("## Author statements")]
    assert 3000 <= word_count(body) <= 7000


def test_open_data_and_author_statements_are_present():
    for identifier in ["PRJNA529500", "PRJNA745492", "10.5281/zenodo.10369064", "NC_000915.1", "U27270.1"]:
        assert identifier in MANUSCRIPT
    for heading in [
        "## Author statements",
        "### Author contributions",
        "### Conflicts of interest",
        "### Funding information",
        "### Ethical approval",
        "### Consent for publication",
        "### Acknowledgements",
    ]:
        assert heading in MANUSCRIPT
    assert "permanent repository URL and version DOI required before submission" in MANUSCRIPT


def test_ai_disclosures_cover_methods_and_acknowledgements():
    assert "### Generative AI assistance" in MANUSCRIPT
    acknowledgements = MANUSCRIPT[MANUSCRIPT.index("### Acknowledgements") :]
    assert "OpenAI Codex" in acknowledgements
    assert "30-31 August 2026" in MANUSCRIPT


def test_authorship_lock_is_unchanged():
    author_line = 'author: "Benteng Ma; Bing Chen; Ting Cai; Xiao-ming Liu; Fen Wang"'
    assert author_line in MANUSCRIPT
    assert "Benteng Ma and Bing Chen contributed equally" in MANUSCRIPT
    assert "No other author is designated as a first author" in MANUSCRIPT
    assert "Xiao-ming Liu" in MANUSCRIPT and "Fen Wang" in MANUSCRIPT


def test_credit_funder_and_all_author_statements_are_finalized():
    contribution = section(MANUSCRIPT, "Author statements")
    for author in ["Benteng Ma", "Bing Chen", "Ting Cai", "Xiao-ming Liu", "Fen Wang"]:
        assert f"{author}:" in contribution
    for role in ["Conceptualization", "Data curation", "Formal analysis", "Funding acquisition", "Software", "Supervision", "Writing - original draft", "Writing - review and editing"]:
        assert role in contribution
    assert "All authors meet the authorship criteria, approved the final manuscript" in contribution
    assert "The funders had no role in study design" in contribution
    assert "roles have not been inferred" not in MANUSCRIPT
    assert "roles have not been inferred" not in PORTAL


def test_figure_panel_style_is_lowercase_parenthesized():
    legends = "\n".join(line for line in MANUSCRIPT.splitlines() if line.startswith("![Figure"))
    assert not re.search(r"\([A-Z](?:[,–-][A-Z])?\)", legends)
    for number in range(1, 10):
        assert f"![Figure {number}." in legends
    figure_code = (ROOT / "scripts" / "make_manuscript_figures.py").read_text(encoding="utf-8")
    assert 'journal_label = f"({label.lower()})"' in figure_code


def test_no_running_title_or_page_header_and_three_line_tables():
    assert "Running title" not in MANUSCRIPT
    assert "Transportability of H. pylori resistance genotyping" not in MANUSCRIPT
    builder = (ROOT / "scripts" / "build_submission_docx.py").read_text(encoding="utf-8")
    assert "configure_three_line_table" in builder
    main_call = builder[builder.index('root / "manuscript/manuscript.qmd"') : builder.index('root / "manuscript/supplementary_methods.qmd"')]
    assert "None," in main_call


def test_journal_specific_submission_materials_are_complete():
    assert "Genomics, Epidemiology and Evolution of Campylobacter, Helicobacter and Related Organisms" in COVER
    assert "Open Data policy" in COVER
    assert "Microbial Genomics portal metadata worksheet" in PORTAL
    assert "Hard stops before submission" in PORTAL
    for requirement in ["Data Summary", "CRediT", "three-line", "Open Data", "£2,203"]:
        assert requirement in CHECKLIST
