import ast
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def title_calls(path: Path, function_names: set[str]) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    violations: list[str] = []
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef) or node.name not in function_names:
            continue
        for call in ast.walk(node):
            if not isinstance(call, ast.Call):
                continue
            if isinstance(call.func, ast.Attribute) and call.func.attr in {"set_title", "suptitle"}:
                violations.append(f"{node.name}:{call.lineno}:{call.func.attr}")
    return violations


def test_main_figure_panels_do_not_use_embedded_subtitles() -> None:
    expanded = ROOT / "scripts/make_expanded_manuscript_figures.py"
    transport = ROOT / "scripts/make_transport_shift_figures.py"
    expanded_functions = {
        "figure1",
        "figure2",
        "forest",
        "figure3",
        "figure4",
        "confusion_matrix",
        "figure5",
        "gate_grid",
        "figure6",
        "figure7",
        "figure8",
    }
    assert title_calls(expanded, expanded_functions) == []
    assert title_calls(transport, {"figure9"}) == []


def test_supplementary_figure_panels_do_not_use_embedded_subtitles() -> None:
    expanded = ROOT / "scripts/make_expanded_manuscript_figures.py"
    transport = ROOT / "scripts/make_transport_shift_figures.py"
    callability = ROOT / "scripts/make_23s_callability_sensitivity_figure.py"
    coverage = ROOT / "scripts/make_coverage_selective_stress_figure.py"
    assert title_calls(expanded, {"supplementary_figures"}) == []
    assert title_calls(transport, {"supplementary_s4", "supplementary_s5", "supplementary_s6"}) == []
    assert title_calls(callability, {"heatmap", "main"}) == []
    assert title_calls(coverage, {"main"}) == []


def test_figure_manifest_contains_only_submission_figure_formats() -> None:
    manifest = ROOT / "results/source_data/main_figure_manifest.csv"
    with manifest.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    assert rows
    assert all(Path(row["file"]).suffix.lower() in {".pdf", ".png", ".tiff"} for row in rows)
