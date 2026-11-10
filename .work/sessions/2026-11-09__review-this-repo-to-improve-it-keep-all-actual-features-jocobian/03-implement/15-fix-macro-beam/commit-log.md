# Ticket 15: fix-macro-beam

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo` (output in `.scratch/tests/<session>/evidence/15-verify.txt`):

- pytest: 527 passed, 1 failed. The only failure is the expected golden one, `tests/golden/test_golden.py::test_beam_matrices_match_golden` (first mismatch `constant__2__matrix`: ConstantBeam ON at `stop` now closes with `[stop, 0]`). Goldens in `tests/golden/data/*` are untouched and NOT regenerated.
- `test_chunked_pulsed_run_matches_golden` still passes (PulsedBeam matrices and intervals are bit-identical to base for every golden window).
- `ruff check radiopyo`: All checks passed. `mypy radiopyo`: no issues in 42 source files.
- Baseline was 502 passed; 25 new tests (tests/physics/test_beams.py, test_timers.py), 4 TODO(ticket 15) tests in test_timers.py updated to the fixed behaviour.
- Before/after report: `.scratch/tests/<session>/evidence/15-before-after.md` (row counts, first/last rows and one-line causes for every golden beam case). Re-record later with `PYTHONHASHSEED=1 python tests/golden/record_golden.py --force --only beams` (runs not needed; never hdf).

## Review Notes

- Parent: base fdb3b6e4d683f7c1ad862537a64dc547cdb45843; pre-existing changes: gps session state only.

- Ticket 11 is merged here: its AC1 (one generic `make_matrix` walking timer transitions up to `min(stop, time_at_max_dose)`, the four hand-written versions removed; `make_time_intervals` shared) is DONE in this ticket. Its AC2 (bit-identical quirks) was dropped (the quirks are the bugs fixed here); its AC3 (chunked golden) holds.
- AC1 (walrus, duty): `MacroTimer.micro_pulses_per_train` / `duty_cycle` rewritten on a tolerant `split_periods` helper (timer.py); train 0.205 / micro period 0.1 / on_time 0.01 gives 2.5 pulses and duty 0.025 (macro period 1.0), tested.
- AC2: macro without `max_dose` terminates and respects `stop` (generic walk); micro start time is passed to the macro timer (and `timer` is the macro subtimer, one object); `MacroPulsedBeam.update_beam_param` rebuilds both (macro names period/frequency/on_time/start_time and micro_* names); no duplicate row at train start (test checks strictly increasing times). `MacroTimer.time_till_next_ON` no longer shifts `start_time` twice (also fixed `time_till_next_train` before start, `time_till_next_OFF` cut by train end).
- AC3: `tot_pulses` and average dose rate use `time_at_max_dose - start_time`; `max_dose=inf` gives `peak * duty_cycle` / inf (no nan); `Timer.copy` and `MacroTimer.copy` keep `start_time`. Macro `at()` used the micro start as beam start; it now relies on the macro timer.
- AC4: tolerance (`COUNT_RTOL` = 1e-9 relative) in pulse/train counting, in the macro whole-train max dose test and `math.isclose` in `PulsedBeam.from_dose_rate_peak` (n_pulse vs max_dose).
- AC5: before/after report written, goldens not regenerated; parent to show it to the user.
- Seed: `BeamCollection.merge` raised TypeError (`s[name][0] |= beam` on a tuple item); now `.update`. Tests added (same type override, other type replaced, macro, operands untouched).
- ConstantBeam / SinglePulseBeam changes (reported per the brief): ON at `stop` < time_at_max_dose now closes with `[stop, 0]` (was `[stop, dose_rate]`): golden constant windows 2 and 5, singlepulse window 4. Also, zero-length ON rows `[T, dr]`/`[T, 0]` duplicates at stop == start_time or stop == T are no longer emitted (constant, e.g. window (0, 0.2) gives 2 rows instead of 3).
- `make_time_intervals` rewritten generically (OFF rows that are neither first nor last are shifted by GROWTH_TIME when the next row is further than GROWTH_TIME). Second column is now strictly 0/1 (old: raw dose rate for 2-row matrices; first row flagged ON in 3-row matrices even when OFF). PulsedBeam output on all 156 compared windows + all golden windows: matrices bit-identical; intervals identical for golden windows, and different only for short windows in the old buggy 2-/3-row cases (extra comparison in `.scratch/tests/<session>/diff15.py`).
- Macro `time_at_max_dose` for whole trains is now the end of the last micro pulse (was the end of the train, 0.1 ms+ later). Not hit by the golden macro case (2.5 trains), whose scalars change only through the average dose rate (A1).
- Public API kept: constructors, classmethods, `at(t).dose_rate()`, `make_matrix`, `make_time_intervals`, `copy`, `update_beam_param`, `BEAM_TYPES`. Added private `_transitions`, `_first_on_timer`, `_update_timers`; `timer.split_periods` and `COUNT_RTOL` (not in `__all__`).
- CHANGELOG: Unreleased / Fixed entry added.

- Parent acceptance: implementer opus/high (also did ticket 11 AC1: one generic make_matrix); reviewer opus/high found one P3 (micro_start_time >= train on_time looped forever), fixed by parent with a ValueError in MacroTimer.new_macropulse plus test. USER APPROVED before/after (evidence/15-before-after.md); beams.npz re-recorded with PYTHONHASHSEED=1 --only beams (runs/hdf untouched). Full suite 529 passed, ruff 0, mypy 0. Gap: no chunked MacroPulsedBeam run golden.

## Blockers / Challenges

- The Bash tool rejected heredocs holding triple-quoted Python; edits were done through files written with the Write tool and spliced in. No impact on the result.
- Remaining gap: golden beam data must be re-recorded after user approval (one expected failing test). The chunked-run golden for a MacroPulsedBeam does not exist, so macro chunked runs are not covered by goldens.

## Commits

- 10978f3 fix: generic beam make_matrix, macro beam and timer fixes, average dose rate from start_time (tickets 15, 11)

## Time Spent

51m (ticket-start.js to ticket-complete.js)
