from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CREDIT = "*[maujabur/bean-forming-lab](https://github.com/maujabur/bean-forming-lab) · made with Kiro*"
DOCS = [ROOT / "README.md"] + sorted((ROOT / "docs").glob("*.md"))


@pytest.mark.parametrize("path", DOCS, ids=lambda p: p.name)
def test_credit_line_follows_the_title(path):
    lines = path.read_text(encoding="utf-8").split("\n")
    assert lines[0].startswith("# "), "first line must be the level-1 title"
    assert lines[1] == "" and lines[2] == CREDIT and lines[3] == "", "credit must sit right below the title"
