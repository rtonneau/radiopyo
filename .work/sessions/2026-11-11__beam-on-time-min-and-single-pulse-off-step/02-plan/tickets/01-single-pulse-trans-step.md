# Ticket 01: single-pulse-trans-step

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] In `UnitCell.prepare_single_pulse`, the beam OFF chunks before tm (leading OFF chunk and gaps between pulses) carry `max_step_size_trans`; the last OFF chunk keeps `max_step_size_final`; the ON -> OFF transition keeps `max_step_size_trans`.
- [ ] Docstring of `max_step_size_trans` / `max_step_size_final` updated accordingly.
- [ ] A test checks the max step of each chunk for a pulsed beam starting after t_span[0].

**Files to Touch:**
- `radiopyo/simulation/unit_cell.py`
- `tests/simulation/` (existing prepare_single_pulse tests)

**Verification Step:**

Run:
```bash
python -m pytest tests/simulation -q
```

Expected:
All tests pass (new chunk max-step test included).

**Notes:**

Line ~920: `max_step_size_on if on else max_step_size_final` -> `max_step_size_trans`.
