# Ticket 15: fix-macro-beam

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] Walrus precedence fixed (`(r := x % p) > 0`): train on_time 0.205, micro period 0.1, micro on_time 0.01 gives 2.5 pulses per train and duty 0.025; the no-op `tot_on_time + res` fixed.
- [ ] `MacroPulsedBeam` without `max_dose` terminates and respects `stop`; micro start time passed to the macro timer; `update_beam_param` keeps both timers in sync; no duplicate row at train start.
- [ ] `tot_pulses` and average dose rate account for `start_time`; `max_dose=inf` gives no `nan`; `copy()` keeps `start_time`.
- [ ] Float comparisons in `n_pulse` vs `max_dose` checks use a tolerance.
- [ ] Before/after report for the macro-pulsed golden case shown to the user before golden regeneration.

**Files to Touch:**
- `radiopyo/physics/beam.py`
- `radiopyo/physics/timer.py`
- `tests/physics/test_beams.py`
- `tests/physics/test_timers.py`
- `tests/golden/data/*`
- `CHANGELOG.md`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass after approved golden regeneration.

**Notes:**

Stop for user approval before regenerating goldens.
