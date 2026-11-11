# Ticket 02: build-time-short-on-time-error

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`: 732 passed in 103.98s (the two golden time-grid tests listed in CLAUDE.md pass too); ruff: one RUF022 (`__all__` order) fixed, then all checks passed; mypy: no issues in 42 files. New file `tests/physics/test_min_beam_on_time.py`: 27 passed. Output in .scratch/tests/<session>/02-verify.txt.

## Review Notes

Base 7e1eaaac1bd9e8b9215b80b893d08246a21444e2, clean tree apart from the gps session record, inline. AC1: `ShortBeamOnError(BeamError, ValueError)` in physics/exceptions.py, exported. AC2: `MIN_BEAM_ON_TIME = 1e-6` in physics/beam.py (in `__all__`), `UnitCell.MIN_BEAM_ON_TIME` documented next to T_ZERO. AC3/AC4: check in `ParticleBeam.__post_init__` (`_configured_on_times` / `_check_on_times`), so every constructor, classmethod, from_dict and update_beam_param path goes through it: pulsed timer on_time; Macro overrides to check micro_on_time and the train on_time; constant / single pulse (period inf) with finite max_dose check time_at_max_dose - start_time. Relative tolerance 1e-9 so a computed on_time of exactly 1e-6 (1 Gy / 1e6 Gy/s) is accepted. A pulse cut by max_dose (partial last or first pulse) is not checked (ticket 03). AC5: message gives name, value and threshold. AC6: raise/accept tests per beam type, update_beam_param, ValueError compatibility. Also: examples/bench_fenton_stall.py grid now skips points with on_time < MIN_BEAM_ON_TIME (every 1e7 Gy/s point; previously only those <= T_ZERO). Diff inspected against the ACs; inline ticket, no independent review (localized to beam construction, covered by tests). Note: update_beam_param that raises leaves the beam with the new timer (same as the existing errors raised from __post_init__).

## Blockers / Challenges

None

## Commits

- 0ce58c8 feat: raise ShortBeamOnError for beam ON times shorter than 1e-6 s (ticket 02)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
