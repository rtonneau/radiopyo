# Ticket 05: docs

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q -p no:cacheprovider && ruff check radiopyo && mypy radiopyo`: 585 passed; ruff and mypy report no error.

## Review Notes

README (use_jac default/methods, enzymatic depletion block), CHANGELOG (Added/Changed/Fixed, stale NotImplementedError lines removed from Known issues and Removed), CLAUDE.md (solver and driver sections) updated; no code touched.

Parent: base 4fbae6cd49aced0ee58871e9e09ac9eb6e95527d, no pre-existing changes. Independent review skipped: docs-only diff read in full by the parent; no stale NotImplementedError/'no jacobian' statements left. Accepted (implementer: haiku, low).

## Blockers / Challenges

None

## Commits

- 8750a82 docs: analytic Jacobian, use_jac default and enzymatic depletion block (ticket 05)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
