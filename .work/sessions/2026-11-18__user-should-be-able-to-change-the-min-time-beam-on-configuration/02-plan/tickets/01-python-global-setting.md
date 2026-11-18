# Ticket 01: python-global-setting

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `set_min_beam_on_time(value)` accepts any positive finite float, raises `ValueError` otherwise; `None` restores 0.4e-6 s; `min_beam_on_time()` returns the current value; both exported from `radiopyo`.
- [ ] `ParticleBeam._check_on_times` uses the current value (beam built after a change follows it); the `MIN_BEAM_ON_TIME` constant stays as the default.
- [ ] `UnitCell.MIN_BEAM_ON_TIME` is a read-only view of the global value (instance and class access).
- [ ] Run-time short ON span test uses max(setting, T_ZERO); setting a value < T_ZERO logs one warning; a run with such a value and a span <= T_ZERO completes with the analytic path.
- [ ] `backend="rust"` with a non-default value raises `RustBackendError` (interim guard); default results unchanged (golden tests).

**Files to Touch:**
- `radiopyo/physics/beam.py`
- `radiopyo/simulation/unit_cell.py`
- `radiopyo/rust.py`
- `radiopyo/__init__.py`
- `tests/physics/test_min_beam_on_time.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All tests pass (new ones included), ruff "All checks passed!", mypy "Success: no issues found".

**Notes:**

Run from the session checkout with PYTHONPATH set to it and the radiobio_env python. `is_short_on_span(t_span, min_on_time)` already takes the threshold as an argument.
