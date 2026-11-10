# Ticket 10: beam-split

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo` -> 450 passed (baseline 410: -1 test_within_train removed with the deleted method, +41 new in tests/physics/test_beams.py), `All checks passed!`, `Success: no issues found in 42 source files`. Output kept in `.scratch/tests/2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian/evidence/ticket10-verify.txt`. Golden tests (incl. beams.npz) and the 54 timer tests pass without any change to golden data/tolerances.

## Review Notes

- Parent: base 9227d70c8c3f1c77cf32f09c89652f97bfa1965c; pre-existing changes: gps session state only.

Implementer self-check against the Acceptance Criteria (no review rounds yet):

- Split: `physics/timer.py` (TimerState, BaseTimer, Timer, MacroTimer), `physics/beam.py` (TimeMessage, ParticleBeam + 4 beams, `BEAM_TYPES`), `physics/beam_collection.py` (BeamCollection). `beam.py` re-exports timers and BeamCollection (import at the bottom, after BEAM_TYPES; `beam_collection` imports `BEAM_TYPES` lazily inside `build_beam_from_dict` to avoid the cycle). `unit_cell.py` imports from the new modules; other importers (`parser/utils.py`, `ode_result.py`, `radiopyo/__init__.py`) still use `physics.beam` and are covered by the re-export.
- Registry: `BEAM_TYPES` lookup on `beam_type.lower()`; the `split("beam")` hack removed in `as_dict` and `merge` (`beam_type().lower()` gives the same keys). Unknown type still raises NotImplementedError. Added `ParticleBeam.from_dict` (raises NotImplementedError) so the registry is typed.
- No mutation: `PulsedBeam.from_dict`, `MacroPulsedBeam.from_dict`, `build_beam_from_dict` and `from_config_dict` work on copies (tests deep-compare the input).
- Timer cleanup: removed `Timer.time_till_next_period` (identical to base), `MacroTimer.period_index` (identical to base), `MacroTimer.time_left_in_train` (identical to `time_till_next_train`, kept as it is used), phase expression factored into `BaseTimer._phase` (used by `time_till_next_period`, `Timer.time_till_next_OFF`, `_relative_time_pulse`). Removed `ParticleBeam.delivery`, both `_compute_dose_rate_peak`, `MacroTimer.within_train`, `micro_pulses_left_in_train`, `micro_pulse_index`, `GROWTH_FACTOR` (grep over radiopyo/, tests/, examples/, README.md, notebooks: only test_timers.py used within_train / micro_pulse_index). `GROWTH_TIME` is now a `ClassVar` instead of being set in `__post_init__`.
- Existing test edit (tests/physics/test_timers.py): deleted `test_within_train` and the three `micro_pulse_index` assertions in `test_period_indexes`, because those methods were removed as unused; no surviving equivalent exists. Nothing else changed there (imports from `radiopyo.physics.beam` kept, which also exercises the re-export).
- Crash fixes: `PulsedBeam.update_beam_param(period=...)` (`kwargs.pop("frequency", None)`); unlabelled bad config in `PulsedBeam.from_dict` now raises ParameterError (`param.get('label', 'default')`).
- `state_at` docstrings (BaseTimer, Timer, MacroTimer) state that a time exactly at pulse end is ON; also fixed the Timer class docstring example that said `state_at(0.50)` is isOff.
- Constructor positional order unchanged (tests: `from_dose_rate_peak(1000, 1e-3, 5e-6, 10)` -> n_pulse=10, plus PulsedBeam/MacroPulsedBeam constructors).
- Layering inversion: added `radiopyo/exceptions.py` holding `RadioPyoError`; `simulation/exceptions.py` re-imports it (still importable from there and in its `__all__`); `physics/exceptions.py` and `chemistry/exceptions.py` (one-line change, outside the ticket's file list) now import from the root module.
- Left as is on purpose (ticket 15): walrus precedence, no-op `tot_on_time + res`, `tot_pulses`, `copy()` dropping start_time, Macro hang.

- Parent acceptance: independent reviewer (sonnet/high): no findings, old-vs-new behaviour script on 10 beam configs identical (3017 lines); parent re-ran pytest (450 passed), ruff 0, mypy 0. RadioPyoError moved to radiopyo/exceptions.py (layering fixed). Seeds for ticket 15: BeamCollection.merge TypeError on same label+type; MacroTimer.time_till_next_ON double start_time shift. test_timers: removed within_train/micro_pulse_index asserts (methods deleted).

## Blockers / Challenges

Circular import between `beam.py` (re-export) and `beam_collection.py` (needs `BEAM_TYPES`): solved with a bottom-of-module import in beam.py and a lazy import in `build_beam_from_dict`.

Pre-existing bugs found and NOT fixed (out of scope, suggest ticket 15): (1) `BeamCollection.merge` raises `TypeError: 'tuple' object does not support item assignment` when both collections define the same label with the same beam type (`s[name][0] |= beam` on a tuple); (2) `MacroTimer.time_till_next_ON` passes an already relative time to `time_till_next_train`, which shifts it by start_time again (wrong when start_time != 0). `chemistry` and `physics` exceptions now depend only on `radiopyo.exceptions`.

## Commits

- cd77b2b refactor: split beams into timer, beam, beam_collection with BEAM_TYPES registry (ticket 10)

## Time Spent

9m (ticket-start.js to ticket-complete.js)
