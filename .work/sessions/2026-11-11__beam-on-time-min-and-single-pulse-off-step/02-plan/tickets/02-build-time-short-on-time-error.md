# Ticket 02: build-time-short-on-time-error

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `ShortBeamOnError(BeamError, ValueError)` in `radiopyo/physics/exceptions.py`, exported.
- [ ] A module constant `MIN_BEAM_ON_TIME = 1e-6` in `physics` and `UnitCell.MIN_BEAM_ON_TIME` (same value, documented next to `T_ZERO`).
- [ ] `PulsedBeam` (constructor, classmethods, `from_dict`, `update_beam_param`) and `MacroPulsedBeam` (`on_time`, `micro_on_time`) raise `ShortBeamOnError` for an ON time < 1e-6 s; exactly 1e-6 s is accepted.
- [ ] `ConstantBeam` / `SinglePulseBeam` with a finite `max_dose` raise it when `time_at_max_dose - start_time` < 1e-6 s.
- [ ] The message gives the ON time and the threshold.
- [ ] Tests for each beam type (raise below, accept at the threshold), and `pytest.raises(ValueError)` still catches it.

**Files to Touch:**
- `radiopyo/physics/exceptions.py`
- `radiopyo/physics/beam.py` (and `timer.py` if the check fits better there)
- `radiopyo/simulation/unit_cell.py`
- `tests/physics/`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
Only the two known golden time-grid failures; ruff and mypy report zero errors.

**Notes:**

Grep tests/ and examples/ for on_time values < 1e-6 first (e.g. `examples/bench_fenton_stall.py` skips 1e7 Gy/s points with on_time <= T_ZERO; make it skip on_time < MIN_BEAM_ON_TIME). The partial last pulse cut by max_dose is NOT an error (ticket 3).
