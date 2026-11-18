# Ticket 01: python-global-setting

**Status:** ✅ Done

## Local Test Result

`pytest -q && ruff check radiopyo && mypy radiopyo`: 866 passed; ruff all checks passed; mypy no issues in 43 files.

## Review Notes

- `rp.set_min_beam_on_time(value|None)` / `rp.min_beam_on_time()` in `physics/beam.py`: ValueError for non-positive/non-finite, warning below T_ZERO.
- `_check_on_times` and `is_short_on_span` (default None) read the current value.
- `UnitCell.MIN_BEAM_ON_TIME` is a read-only descriptor view (setting raises AttributeError); the Rust guard in `rust.py` therefore refuses a non-default value until ticket 03.
- `_raw_run` takes the analytic path for spans <= T_ZERO (with rounding slack, shared `_span_slack`).
- New tests: `tests/physics/test_min_beam_on_time_setting.py` (12).

## Blockers / Challenges

A 0.1 us span computed from chunk boundaries rounds slightly above T_ZERO and hit `_apply_initial_burst`'s ValueError: the T_ZERO floor now uses the same rounding slack as `is_short_on_span`.

## Commits

- 267f36e feat: configurable minimum beam ON time, rp.set_min_beam_on_time (ticket 01)

## Time Spent

5m (ticket-start.js to ticket-complete.js)
