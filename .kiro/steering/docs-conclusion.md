---
inclusion: always
---

# Conclusion section in study reports

Every study report in `docs/` (static focusing, virtual microphone, field maps, parametric, broadband, three-dimensional and array-optimization studies, and any new one) ends with a `## Conclusion` section, as the **last** section of the file.

- One short paragraph (about 5 to 8 lines): the answer to the study's main questions with the few numbers that matter, what remains uncertain, and the model's assumptions (free field, ideal point sources, predictions and not measurements).
- Do not introduce new results in the conclusion; every number must already appear in the report.
- Reference and convention documents (`array-configuration.md`, `grid-geometry.md`, `hardware/*.md`, `index.md`) do not need one.
- When adding a new report, add its name to the list in `tests/test_docs_conclusion.py`.
