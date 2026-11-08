---
inclusion: always
---

# Documentation credit convention

Every Markdown document of the project (`README.md` and every `docs/*.md`) carries a credit line **immediately below the level-1 title**, separated from the title and from the following text by one blank line:

```markdown
# Document title

*[maujabur/bean-forming-lab](https://github.com/maujabur/bean-forming-lab) · made with Kiro*

First paragraph...
```

Rules:

- Use exactly the line above (italic, repository link, middle dot, "made with Kiro").
- Put it right after the `# ` title, not at the end of the file and not after the first paragraph.
- Apply it to every new `.md` file in `docs/` and to the README when they are created, and keep it when a document is regenerated or rewritten (including documents produced from templates or scripts).
- Do not add it to generated data files, JSON, CSV or images.
- `tests/test_docs_credit.py` checks the convention; run the full suite before finishing.
