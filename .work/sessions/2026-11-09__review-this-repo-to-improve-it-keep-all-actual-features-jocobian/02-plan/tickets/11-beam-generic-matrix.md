# Ticket 11: beam-generic-matrix

**Model:** heavy
**Effort:** xhigh

**Acceptance Criteria:**
- [ ] One generic `make_matrix` (and time-interval builder) walks timer transitions up to `min(stop, time_at_max_dose)`; the four hand-written versions are removed.
- [ ] Beam-matrix goldens match exactly for Constant, SinglePulse, Pulsed and MacroPulsed beams (including the current quirks to be fixed in ticket 15, e.g. duplicate row at train start).
- [ ] Chunked-run goldens pass.

**Files to Touch:**
- `radiopyo/physics/beam.py`
- `radiopyo/physics/timer.py`
- `tests/physics/test_beams.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass, zero lint/type errors.

**Notes:**

If reproducing a quirk exactly makes the generic version contorted, stop and report rather than bending the design; it may be better to merge this ticket into 15.
