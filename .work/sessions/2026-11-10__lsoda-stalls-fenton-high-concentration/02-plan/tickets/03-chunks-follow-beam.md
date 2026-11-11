# Ticket 03: chunks-follow-beam

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] `make_time_intervals` no longer shifts pulse-end rows by GROWTH_TIME: ON intervals end exactly at the pulse end (GROWTH_TIME removed or unused, documented).
- [ ] In chunked runs (`prepare_chunked_run`, `prepare_single_pulse`) the solver sees a constant dose rate per chunk (the beam dose rate inside that interval); `run()` keeps `beam.at(t)`.
- [ ] The initial burst still uses the chunk's dose rate; `ODEResult` reconstruction of rates (which reads the beam) stays consistent with what was integrated.
- [ ] Tests for every beam type (constant with max_dose, pulsed, single pulse, macro-pulsed): chunk boundaries equal the beam transitions and the per-chunk dose rate equals the beam dose rate.

**Files to Touch:**
- `radiopyo/physics/beam.py`
- `radiopyo/simulation/unit_cell.py`
- `tests/physics/test_beam_intervals.py`
- `tests/simulation/test_chunked_dose_rate.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
all tests pass (golden drift reported, not re-recorded), zero ruff/mypy errors.

**Notes:**

Prototype: `patch_const.py` (constant dose rate read at chunk midpoint) + `GROWTH_TIME = 0`. Without the exact cut the constant dose rate overdoses the 1 µs tail (ROOH +9 %): both changes go together. Check how `ODEResult` reads the beam for reconstruct_reactions at chunk boundaries.
