# Ticket 04: benchmark-and-default

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] Benchmark script and report (wall time, nfev, njev, success, analytic vs FD) on the bundled configs saved under `.scratch/tests/2026-11-10__analytic-jacobian/`.
- [ ] Gate: if the analytic Jacobian is not at least as fast and robust, stop BLOCKED and report to the user without changing the default.
- [ ] Otherwise `use_jac` defaults to True in `RunOptions`, `run`, `prepare_chunked_run`, `prepare_single_pulse`, `run_till_max_dose` paths; goldens pass against unchanged data within a tolerance written in the test.

**Files to Touch:**
- `radiopyo/simulation/unit_cell.py`
- `tests/golden/`
- `.scratch/tests/2026-11-10__analytic-jacobian/benchmark_jacobian.py`

**Verification Step:**

Run:
```bash
python .scratch/tests/2026-11-10__analytic-jacobian/benchmark_jacobian.py && python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
Benchmark report written; all tests pass with the new default.

**Notes:**

Do not modify tests/golden/data. Tolerance chosen from observed differences, stated in the test.
