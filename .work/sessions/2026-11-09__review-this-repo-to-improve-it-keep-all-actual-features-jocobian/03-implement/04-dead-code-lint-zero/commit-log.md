# Ticket 04: dead-code-lint-zero

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`: 262 passed (incl. tests/golden); `All checks passed!`; `Success: no issues found in 36 source files`. Also `python -W error -c "import radiopyo"` is clean (fixed invalid `\ ` escapes in docstrings).

## Review Notes

- Parent: base d0b0f53bf3b612921d61984ae379bdfa4e90147d; pre-existing changes: gps session state only.

All listed removals done (grep showed no remaining references in package, tests, README, notebook). `_set_mass` kept as a no-op with a comment (dead body removed). CHANGELOG `[Unreleased] / Removed (breaking)` added; CLAUDE.md RON/run_multi_doses mentions removed. Lint/mypy fixed by deletion, `ruff --fix` (import/`__all__` sorting), narrow `tp.cast` for numpy Polynomial returns and two targeted `type: ignore` (pandas/TypedDict) with reasons. `jacobian()` kept (type-checks). No known behaviour bugs touched; `initial_cc` now returns `float(...)` (same value).

- Parent acceptance: independent reviewer (sonnet): 1 P3 (blank lines before CstReaction) fixed inline by parent; reviewer gaps closed: parent re-ran pytest (262 passed), ruff clean, mypy clean (evidence/04-parent-verify.txt); HDF5 read path covered by tests/golden SimCollection round-trip. Model sonnet (requested medium).

## Blockers / Challenges

`physics/exceptions.py` still imports `RadioPyoError` from `simulation.exceptions` (layering inversion): not trivial to move, left as is.

## Commits

- ad88152 refactor: remove dead code and RON/run_multi_doses/e_LET, ruff and mypy at zero (ticket 04)

## Time Spent

8m (ticket-start.js to ticket-complete.js)
