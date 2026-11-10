# Ticket 02: port-tests

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `tests/chemistry` collects and passes against the current API (`BasicSpecies`/`SimSpecies`, `stoi=`, `RadiolyticReaction.from_product`, species objects in `ListSpecies`), keeping each test's original intent.
- [ ] New tests pin current include-merge behaviour (override of a k-reaction rate from an included file) and current timer behaviour (`Timer`/`MacroTimer` state and transition times on a few hand-computed cases).
- [ ] Tests for behaviour that ticket 12-17 will change are marked with a comment naming that ticket, not encoded as "correct".
- [ ] Full `python -m pytest` passes.

**Files to Touch:**
- `tests/chemistry/conftest.py`
- `tests/chemistry/test_*.py`
- `tests/parser/test_includes.py`
- `tests/physics/test_timers.py`

**Verification Step:**

Run:
```bash
python -m pytest -q
```

Expected:
All tests pass, no collection errors.

**Notes:**

The old tests used `add_as_*(..., multiple=2)` and adding a species twice incremented; now it raises. Adapt the assertion to the current contract rather than changing code.
