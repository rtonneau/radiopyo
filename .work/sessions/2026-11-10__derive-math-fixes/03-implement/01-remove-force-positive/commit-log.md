# Ticket 01: remove-force-positive

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`: 580 passed (golden tests included, no re-record), "All checks passed!", "Success: no issues found in 42 source files".

## Review Notes

Parent: base bfba771fb90adc1086b6d7ce2d740845a10961f4; pre-existing changes: only gps session state (.session-config.json, 03-implement/), owned by the parent.

- AC1: `force_positive` removed from RunOptions, run, prepare_chunked_run, prepare_single_pulse, derive, jacobian (also the unused `SimParameters.force_positive` key and its default in sim_types.py). run_till_max_dose takes **kwargs, so it now raises TypeError too (tested). `grep force_positive radiopyo/` is empty; in tests/ only the new TypeError tests mention it.
- AC2: `zero_depleted_consumption` deleted, removed from `__all__`, its test removed.
- AC3: test_run_options.py, test_jacobian.py (parametrization dropped, depleted-row test rewritten as non-zero-row + FD match), test_golden.py rewritten; new tests assert RunOptions(force_positive=True) and the run methods raise TypeError.
- AC4: golden derive.json not re-recorded; the `beam_on_force_positive` case is dropped from DERIVE_CASES and test_derive_matches_golden now checks actual cases are a subset of the recorded ones (the stale recorded entries are ignored). Parent may want to decide whether to prune them at a later re-record.
- AC5: see Local Test Result.
- Parent review: fresh reviewer (base bfba771) REVIEWED, no findings; confirmed that no load path (TOML parser, SimEnv.from_dict, ODEResult.from_json_config) validates parameter keys, so old configs/stores with force_positive still load. P3 notes accepted (test_reaction_rates.py touched legitimately; stale derive.json force_positive entries pruned at the next approved derive re-record). Parent re-ran verification: 580 passed, ruff clean, mypy clean. Accepted.

## Blockers / Challenges

None

## Commits

- 3375410 refactor!: remove force_positive (ticket 01)

## Time Spent

7m (ticket-start.js to ticket-complete.js)
