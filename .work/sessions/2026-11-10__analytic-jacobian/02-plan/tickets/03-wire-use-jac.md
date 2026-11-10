# Ticket 03: wire-use-jac

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `_raw_run` no longer raises NotImplementedError; `jac=jacobian` with the same args as `derive` is passed only for LSODA, BDF and Radau.
- [ ] `use_jac=True` with RK45 runs with no warning about jac.
- [ ] Tests: use_jac True vs False agree within tolerance on LSODA, BDF, Radau, and on a chunked run. Default stays False; docstrings updated.

**Files to Touch:**
- `radiopyo/simulation/unit_cell.py`
- `tests/simulation/test_use_jac.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass, zero lint/type errors.

**Notes:**

Remove or update existing tests expecting NotImplementedError.
