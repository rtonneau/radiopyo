# Ticket 03: runtime-short-on-spans-docs

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`: 746 passed in 110.92s (the two golden time-grid tests expected to fail pass, the golden data was already re-recorded); ruff: All checks passed; mypy: no issues in 42 source files.

## Review Notes

Base 91daae24e6765b901130b17827f48c091bb833e9, inline. AC1: `_raw_run` sends a beam ON span shorter than `MIN_BEAM_ON_TIME` (`is_short_on_span`, slack for chunk-boundary rounding) to `_short_on_span_step`: yield over the span added analytically, one point at the span end; covered for chunked (last/first pulse cut by max_dose, micro pulse cut by its train end, t_span clipping), single pulse and plain run. AC2: warning logged; `ODEResult.short_on_spans` concatenated by `__add__`, empty after an HDF5 round trip. AC3: e_aq radiolytic production equals the expected dose (rel 1e-9) for max_dose 10.1, 10.5 and 0.1; the step dose rate of a short span is 0 so it is counted once. AC4: `reconstruct_reactions` / `to_pandas` tested. AC5: golden tests pass unchanged. AC6: glossary *Short ON span*, ADR 0008, changelog fragment (the Known issue line removed from the lsoda fragment), README and CLAUDE.md updated; the outdated CLAUDE.md note on two failing golden tests dropped. Deviation: `ShortOnSpan` also carries the dose rate (start, stop, dose_rate), needed to recompute the production integral. `test_run_options.py`: spans of 1e-6 s moved to 2e-6 s (a 1e-6 - 1e-8 span would now be short). New tests: `tests/simulation/test_short_on_span.py` (14).

## Blockers / Challenges

`test_raw_run_uses_options` needed a beam-ON empty span longer than the threshold after widening (2e-6). `_apply_initial_burst` keeps its ValueError guard, now unreachable for beam ON spans.

## Commits

- 60c4aac feat: short beam ON spans added analytically at run time, docs (ticket 03)

## Time Spent

13m (ticket-start.js to ticket-complete.js)
