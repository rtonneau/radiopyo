# Ticket 01: solver-failure-error

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `SolverFailedError(RadioPyoError)` exists in `radiopyo/simulation/exceptions.py`, exported.
- [ ] `UnitCell._raw_run` raises it when `sol.success` is False, message holding t_span, the last reached t and scipy's message; `run`, chunked and single-pulse runs therefore raise it.
- [ ] Test: a run forced to fail (e.g. BDF on the reproduction chunk that fails, or a tiny-step monkeypatched solver) raises; successful runs unchanged.
- [ ] CHANGELOG Unreleased bullet.

**Files to Touch:**
- `radiopyo/simulation/exceptions.py`
- `radiopyo/simulation/unit_cell.py`
- `tests/simulation/test_solver_failure.py`
- `CHANGELOG.md`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
all tests pass, ruff and mypy report zero errors.

**Notes:**

Reproduction of a real failure: `.scratch/tests/2026-11-10__lsoda-stalls-fenton-high-concentration/fail.py` (BDF, fe.toml, pulse_dose 1, 1e5 Gy/s). Prefer a fast deterministic test (e.g. monkeypatch solve_ivp to return success False).
