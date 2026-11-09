# Ticket 10: beam-split

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] `physics/timer.py` (`TimerState`, `BaseTimer`, `Timer`, `MacroTimer`), `physics/beam.py` (`ParticleBeam` + 4 beams), `physics/beam_collection.py` (`BeamCollection`); `physics/beam.py` re-exports timers and `BeamCollection` so old imports work.
- [ ] `BEAM_TYPES` registry replaces the type if-chain and the `split("beam")` hack.
- [ ] `from_dict`/`from_config_dict` never mutate the caller's dict.
- [ ] Duplicated timer methods removed (`time_till_next_period`, `period_index`, `time_left_in_train`/`time_till_next_train`, repeated phase expression); unused methods removed (`delivery`, `_compute_dose_rate_peak`, `within_train`, `micro_pulses_left_in_train`, `micro_pulse_index`, `GROWTH_FACTOR`).
- [ ] Crash fixes: `PulsedBeam.update_beam_param(period=...)` KeyError; unlabelled bad config raises `ParameterError`.
- [ ] `state_at` docstring states that a time exactly at pulse end is ON.
- [ ] Positional argument order of every public constructor unchanged (test `from_dose_rate_peak(1000, 1e-3, 5e-6, 10)` gives `n_pulse=10`).
- [ ] Golden tests pass.

**Files to Touch:**
- `radiopyo/physics/beam.py`
- `radiopyo/physics/timer.py`
- `radiopyo/physics/beam_collection.py`
- `radiopyo/physics/exceptions.py`
- `radiopyo/simulation/unit_cell.py` (imports)
- `tests/physics/test_beams.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass, zero lint/type errors.

**Notes:**

Do not fix the walrus precedence, `tot_pulses`, start-time or Macro hang bugs here (ticket 15).
