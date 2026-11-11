# Ticket 03: chunks-follow-beam

**Status:** ✅ Done

## Local Test Result

Command: `python -m pytest -q && ruff check radiopyo && mypy radiopyo` (run as three commands, output in `.scratch/tests/2026-11-10__lsoda-stalls-fenton-high-concentration/evidence/ticket03-verification.txt`).

- pytest: 692 passed, 2 failed. The 2 failures are the expected golden drift, not re-recorded:
  - `tests/golden/test_golden.py::test_beam_matrices_match_golden`: `constant__1__intervals` has `0.205` instead of `0.205001` (no GROWTH_TIME shift of the pulse end).
  - `tests/golden/test_golden.py::test_chunked_pulsed_run_matches_golden`: the chunk boundary times moved (`2e-6` instead of `3e-6`, `1.0002e-2` instead of `1.0003e-2`, ...). Checked separately (`evidence`-less script `cmp_golden.py` in the scratch dir): species at t = 1e-2 s and t = 5e-2 s agree with the golden values to <= 8e-9 relative, so only the time grid drifts.
- ruff check radiopyo: All checks passed.
- mypy radiopyo: Success, no issues in 42 source files.

## Review Notes

Review rounds: 0 (implementer self-check only).

- AC1: `ParticleBeam.GROWTH_TIME` removed; `make_time_intervals` is `make_matrix` with column 1 turned into 0/1, so ON intervals end exactly at the pulse end or at `time_at_max_dose` (docstring cites ADR 0007). Tests: `tests/physics/test_beam_intervals.py` (end equality, no `GROWTH_TIME` attribute); `tests/physics/test_beams.py` updated (`intervals[1, 0] == 2e-6`).
- AC2: new `ConstantDoseRate` (beam interface `at(t).dose_rate()`) and `DoseRateSource` protocol in `radiopyo/physics/beam.py`; solver annotations in `ode_functions_solver.py` widened to `DoseRateSource`. `SimChunk` gets `dose_rate` (default None), filled by `_make_chunks` from the beam read at the chunk midpoint (`chunk_dose_rate`); `_raw_run(..., dose_rate=)` passes `ConstantDoseRate` to the solver. `run()` is unchanged (still reads `beam.at(t)`; `step_dose_rates == [None]`).
- AC3: `_apply_initial_burst` takes the chunk dose rate (falls back to the beam at t0 + T_ZERO when None, as before). `ODEResult` reads the beam at the stored time points, which breaks at a chunk boundary: the last point of an ON chunk lies on the pulse end where `beam.at` is OFF, so half a step of production would be lost. Added `ODEResult.step_dose_rates` (parallel to `n_points`, filled by `append_sol(dose_rate=)`, concatenated in `__add__`) and `_stacked_rates`, which reconstructs a constant-rate step with its own dose rate. Results read from files have an empty list and keep the old behaviour; the HDF5 layer is untouched.
- AC4: `tests/simulation/test_chunked_dose_rate.py` covers constant with max_dose, single pulse, pulsed and macro-pulsed beams (chunk boundaries equal the `make_matrix` transitions, chunk dose rate equals the beam dose rate at several points inside it), `prepare_single_pulse`, the beam object seen by `derive` (spy), the burst dose rate, and the reconstructed rates/integral per ON chunk.
- Changed after self-review: existing spies of `_raw_run` in `tests/simulation/test_run_options.py` now accept the new `dose_rate` keyword; `_timed_solve` annotation widened (mypy).
- Correction pass 1 (parent review round 1):
  - P2: `prepare_single_pulse` now inserts a leading beam OFF chunk `[t_span[0], start_time]` when `t_span[0] < start_time < time_at_max_dose`, so the ON chunk starts at the beam start time and every constant-dose-rate chunk lies inside one beam interval. Test `test_single_pulse_chunks_with_beam_starting_after_t_span` (SinglePulseBeam and ConstantBeam, start_time 1e-3, t_span (0, 1e-2)): boundaries, chunk dose rates `[0, peak, 0, 0]`, beam dose rate inside every chunk, burst at the peak dose rate, `step_dose_rates`.
  - P3: no behaviour change; CHANGELOG Unreleased now has bullets for the pulse-shorter-than-T_ZERO ValueError (under Fixed), GROWTH_TIME removal / exact chunk cuts and constant per-chunk dose rate (under Changed), and `ODEResult.step_dose_rates` (under Added).
  - Re-run: pytest 694 passed, 2 failed (the same 2 expected golden time-grid failures); ruff and mypy clean. Output: `evidence/ticket03-verification-round1.txt`.
- Correction pass 2 (parent review round 2):
  - `prepare_single_pulse` now builds its pre-tm chunks from `beam.make_time_intervals(start=t_span[0], stop=min(tm, t_span[1]))` (as `prepare_chunked_run`: leading OFF chunk, one chunk per pulse, OFF gaps; it no longer reads the timer start time, so macro-pulsed beams work), then the transition and final chunks. Every boundary is clipped to `t_span[1]` and empty chunks are dropped. ON chunks use `max_step_size_on`; the OFF chunks before tm use `max_step_size_final` (the method has no `max_step_size_off`; docstring updated).
  - Tests added to `tests/simulation/test_chunked_dose_rate.py`: multi-pulse `PulsedBeam` (with and without start_time 1e-4) and `MacroPulsedBeam` (start_time 1e-4, max_dose 10) with t_span (0, 3e-3): boundaries equal the beam intervals, chunk dose rate equals `beam.at(t)` sampled inside every chunk, final state of `prepare_single_pulse` equals `prepare_chunked_run` (rtol 1e-5); t_span[1] clip (inside the pulse, inside the transition, exactly at tm).
  - CHANGELOG bullet reworded for the single-pulse chunks.
  - Re-run: pytest 701 passed, 2 failed (same expected golden time-grid failures); ruff and mypy clean. Output: `evidence/ticket03-verification-round2.txt`.
- Gap: the golden data (`beams.npz`, chunked run in `runs.npz`) still hold the GROWTH_TIME grid; re-recording is left to the parent. The README/CHANGELOG are not updated (outside the ticket's file list).

## Blockers / Challenges

None blocking. Notes: the trapezoid integral across two chunks still counts half of the skipped T_ZERO gap between the end of an OFF chunk and the first point of the next ON chunk (pre-existing, unchanged); the per-chunk integral test avoids it. Some files in the working tree are CRLF while the index is LF; the edited files were kept CRLF.

**Blocked (2026-11-11T08:23:24.294Z):** re-review after the one allowed correction pass: prepare_single_pulse gives a wrong constant dose rate for multi-pulse beams (3.5x dose regression) and for MacroPulsedBeam starting after t_span[0]

## Commits

- 93dfc06 fix: chunks follow the beam exactly with a constant dose rate per chunk (ticket 03)

## Time Spent

49m (ticket-start.js to ticket-complete.js)
