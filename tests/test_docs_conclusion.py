from pathlib import Path

import pytest

DOCS = Path(__file__).resolve().parents[1] / "docs"
REPORTS = ["static-focusing", "virtual-microphone", "field-map", "parametric-study", "broadband-study",
           "three-dimensional-field", "array-optimization"]


@pytest.mark.parametrize("name", REPORTS)
def test_every_study_report_ends_with_a_conclusion(name):
    text = (DOCS / f"{name}.md").read_text(encoding="utf-8")
    headings = [line for line in text.split("\n") if line.startswith("## ")]
    assert headings[-1] == "## Conclusion", "the Conclusion must be the last section"
    assert len(text.split("## Conclusion")[1].strip()) > 300
