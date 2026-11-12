# Ticket 04: docs-packaging

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] `pyproject.toml` gets an optional `examples = ["matplotlib"]` extra (plus Jupyter only if needed).
- [ ] The README points to the paper reproduction notebook (examples section) and to the golden re-record command.
- [ ] A changelog fragment for this session lists the notebook, `paper_runs.py`, the golden data and test.
- [ ] `ruff check radiopyo`, `mypy radiopyo`, `python -m pytest` pass.

**Files to Touch:**
- `pyproject.toml`
- `README.md`
- `.work/changelog/2026-11-12__paper-reproduction-notebook.md`

**Verification Step:**

Run:
```bash
ruff check radiopyo && mypy radiopyo && python -m pytest -q
```

Expected:
Zero ruff/mypy errors; all tests pass.

**Notes:**

README.md already has uncommitted user edits: touch only the examples section and do not commit the user's other pending changes (stage hunks or coordinate). Check how the previous session wrote its changelog fragment (`.work/changelog/2026-11-12__chunked-run-scaling.md`).
