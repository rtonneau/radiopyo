# Ticket 01: remove-force-positive

**Effort:** medium

**Acceptance Criteria:**
- [ ] `force_positive` no longer exists in `RunOptions`, `UnitCell.run`, `UnitCell.prepare_chunked_run`, `UnitCell.prepare_single_pulse` (and any other run helper), `derive` or `jacobian`; passing it raises TypeError.
- [ ] `zero_depleted_consumption` is deleted and removed from `__all__`.
- [ ] Tests referring to force_positive are removed or rewritten (test_run_options.py, test_jacobian.py, tests/golden/test_golden.py); a test asserts `RunOptions(force_positive=True)` raises TypeError.
- [ ] Default results unchanged: golden tests pass without re-recording.
- [ ] pytest, ruff check radiopyo and mypy radiopyo report zero errors.

**Files to Touch:**
- `radiopyo/simulation/ode_functions_solver.py`
- `radiopyo/simulation/unit_cell.py`
- `radiopyo/utils/sim_types.py`
- `tests/simulation/test_run_options.py`
- `tests/simulation/test_jacobian.py`
- `tests/golden/test_golden.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All tests pass (no golden re-record), "All checks passed!", "Success: no issues found".

**Notes:**

ADR 0004. Docs (README, CLAUDE.md) are updated in ticket 06; grep for remaining `force_positive` in radiopyo/ and tests/ must be empty. Do not touch .scratch/ or .work/ history files.
