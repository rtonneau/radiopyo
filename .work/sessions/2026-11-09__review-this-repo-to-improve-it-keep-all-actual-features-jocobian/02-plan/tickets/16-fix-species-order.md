# Ticket 16: fix-species-order

**Model:** fast
**Effort:** medium

**Acceptance Criteria:**
- [ ] `ListSpecies.sort` and any other `set` iteration that decides indices use a deterministic order (insertion or sorted); species indices identical across `PYTHONHASHSEED` values (test via subprocess).
- [ ] Before/after report (expected: round-off-level differences only) shown to the user.

**Files to Touch:**
- `radiopyo/chemistry/species_collection.py`
- `tests/chemistry/test_species_collection.py`
- `CHANGELOG.md`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass.

**Notes:**

Goldens are label-keyed, so they should pass without regeneration; if run results move beyond tolerance, report it.
