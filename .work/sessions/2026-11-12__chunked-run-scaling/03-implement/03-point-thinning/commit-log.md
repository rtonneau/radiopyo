# Ticket 03: point-thinning

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`:
786 passed in 39.68s (golden included; 18 new in tests/simulation/test_point_thinning.py); ruff "All checks passed!"; mypy "Success: no issues found in 42 source files".

## Review Notes

Mode inline, base a08e16b.
- `prepare_chunked_run` / `prepare_single_pulse` take `points_per_chunk_on` / `points_per_chunk_off` (default None: unchanged, every solver step stored; checked by test_default_unchanged, bit-identical). Values validated by `check_points_per_chunk` (None or int >= 2, bools rejected -> ValueError).
- `SimChunk.n_points` carries the setting; `_raw_run(n_points=)` replaces `t_eval` by `thinned_times(start, stop, n)` after the initial burst: start, then n-1 log-spaced offsets from span*1e-6 to span (stop exact), np.unique against float collisions; an empty span returns [start] (geomspace cannot take 0).
- Solver steps unchanged, so the chained y0 matches the full run at every chunk end (dense output, tested at atol 1e-6 relative to each species max).
- `ODEResult.thinned` (ORed by `__add__`), set by `SimChunksIter` unless the chunk was a short ON span; production/consumption/species integrals log one warning per result (seed `exact-integrals-thinned-runs` for exact integrals later).
- tests/simulation/test_run_options.py: the signature guard lists the two new keyword parameters and the `_raw_run` spy accepts `n_points` (asserts None by default).
- Self-review (opt-in feature, default path covered by golden + unchanged-output test).

## Blockers / Challenges

A first patch through a shell heredoc failed on quoting (nothing applied); applied through a script file instead.

## Commits

- 7de28dc feat: opt-in point thinning of chunked runs (ticket 03)

## Time Spent

8m (ticket-start.js to ticket-complete.js)
