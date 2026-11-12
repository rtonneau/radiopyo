# Ticket 04: docs-packaging

**Status:** ✅ Done

## Local Test Result

- Verification step: `ruff check radiopyo && mypy radiopyo && python -m pytest -q` (radiobio_env, `PYTHONPATH=.`) exits 0: ruff all checks passed, mypy no issues in 42 source files, 802 passed in 37.9 s (log: `.scratch/tests/2026-11-12__paper-reproduction-notebook/04-verify.log`).

## Review Notes

Checked against the acceptance criteria:
- `pyproject.toml`: new extra `examples = ["matplotlib"]`. Jupyter is not added: the notebook runner (Jupyter, VS Code) brings its own kernel, and `paper_runs.py` / the golden test need nothing beyond the core dependencies.
- README: the Documentation list links the paper reproduction notebook (what it draws, `RECOMPUTE`, the `examples` extra); the Development section explains `tests/golden/` and gives the re-record command of the paper grid (`PYTHONHASHSEED=0 python tests/golden/record_golden.py --force --only paper --workers 8`).
- Changelog fragment: not written by hand in this ticket. `/gps finish` drafts `.work/changelog/2026-11-12__paper-reproduction-notebook.md` through `changelog-prepare.js` / `changelog-apply.js`, so writing it here would duplicate or conflict with that step; the finish payload lists the notebook, `paper_runs.py`, the golden data and test, and the `detail_consumption_of` fix.
- The README edits of the user noted in the plan were committed by the user in b396bd3 before this ticket, so only this ticket's hunks are in the working tree.

## Blockers / Challenges

None.

## Commits

- 3cfa505 docs: link the paper reproduction notebook, examples extra (ticket 04)

## Time Spent

2m (ticket-start.js to ticket-complete.js)
