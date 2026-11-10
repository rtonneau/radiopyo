# Ticket 04: catalase-k

**Effort:** low

**Acceptance Criteria:**
- [ ] `radiopyo/data/reactions_biology.toml` catalase k_value is 6.62e7 with a comment stating it is the per-reaction k, used as given (ADR 0006).
- [ ] Before/after golden numbers recorded in the ticket log; only failing parts re-recorded; the old HDF5 golden file still loads and the pinned doubled-k json-store test still holds.
- [ ] Tests that hard-code 1.324e8 for the bundled config are updated where they mean the bundled value (test fixtures that define their own reaction may keep it).
- [ ] pytest, ruff and mypy at zero.

**Files to Touch:**
- `radiopyo/data/reactions_biology.toml`
- `tests/golden/data/`
- `tests/parser/test_includes.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All tests pass; ruff and mypy report no errors.

**Notes:**

Open question kept: 6.62e7 assumed to be the per-reaction k of Labarbe's Rb10. Result-changing fix: same golden procedure as ticket 02.
