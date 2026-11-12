# Ticket 01: map-helper

**Status:** ✅ Done

## Local Test Result

Ran `cd docs/issues/fenton && python -c "import fenton_map as m; r=m.compute_map(n=2, period=60, cache=None); ..."`.
Output: shapes {'dose_per_pulse': (2,), 'avg_dose_rate': (2,), 'rooh_max_dose': (2, 2), 'rooh_end': (2, 2)}; rooh_end = [[8.38e-4, 8.0e-4], [8.31e-4, 8.0e-4]] (finite), rooh_max_dose finite (~0 for the single-pulse column, physical). About 1 min for 4 points (period 60 s, dominated by the 668-chunk point). Also ran with cache='.map_cache': .npz written (fenton_map_5f727a4b2b25.npz); the second-call load path was not timed separately. Cache dir removed afterwards. ruff not installed in this env (not run).

## Review Notes

compute_map(n, dpp_range, dr_range, period, fe, cache, max_workers) returns the axes plus rooh_max_dose and rooh_end (n_dr x n_dpp, umol/L), NaN when period<=on_time, max dose not reached within period, or solver failure. ProcessPoolExecutor and sha1-keyed .npz cache implemented. compute_point(...) exposed with fe override (Fe_plus_plus via config dict). Banner layout and logging used. Changed after first run: ROOH-at-max-dose index now the last sample at t_max (duplicate times at chunk boundaries). Added repo-root sys.path fallback because radiopyo is not pip-installed in the env.

## Blockers / Challenges

radiopyo not importable from docs/issues/fenton (not installed): fixed by sys.path fallback. Per-point runtime is long for high frequency (hundreds of chunks); a full 12x12 map at 1200 s was not run.

## Commits

- 098fd9b docs: add cached parallel ROOH map helper for the Fenton case (ticket 01)

## Time Spent

7m (ticket-start.js to ticket-complete.js)
