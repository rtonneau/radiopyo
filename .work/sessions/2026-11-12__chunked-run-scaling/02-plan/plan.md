# Implementation Plan

**Session:** chunked-run-scaling
**Date:** 2026-11-12T10:58:33.327Z
**Estimated effort:** 1 day

## Strategy

First remove the quadratic term (lazy, list-backed storage in ODEResult, bit-identical), then trim the per-call overhead of derive/jacobian (identical to rounding), then add the opt-in point thinning on top of the new storage, and finish with the 10k-pulse benchmark and the docs. Each ticket keeps pytest, ruff and mypy green.

## Tickets Overview

- **Ticket 1:** ODEResult stores per-step t/y arrays in lists, concatenated lazily and cached; SimChunksIter takes y0 from the last chunk; O(chunk) append test.
- **Ticket 2:** Trim per-call allocations in derive / jacobian / to_full_molar / mass-action helpers, results identical to rounding.
- **Ticket 3:** Opt-in point thinning (`points_per_chunk_on` / `points_per_chunk_off`) on prepare_chunked_run and prepare_single_pulse, thinned flag and integral warning.
- **Ticket 4:** examples/bench_chunked_scaling.py (up to 10k pulses), README, CHANGELOG, CLAUDE.md.

## Sequencing Rationale

Ticket 1 changes how ODEResult stores points and is the precondition for flat per-chunk cost; ticket 3 builds its per-chunk t_eval on that storage. Ticket 2 is independent but placed before 3 so the benchmark of ticket 4 measures both. Ticket 4 documents the final API and numbers.

## Risks & Mitigation

- **Risk:** code (internal or user notebooks) assigns or mutates `res.t` / `res.y` in place and breaks with properties → **Mitigation:** keep setters (assignment replaces the lists by one block); the cached array is the one returned, so in-place edits of a read array persist until the next append; tests cover `__add__`, `from_dataframe`, HDF5 round trip.
- **Risk:** the derive trim changes results beyond rounding → **Mitigation:** golden tests and test_overshoot_equivalence (rtol 1e-12) must pass unchanged; no golden re-recording allowed.
- **Risk:** t_eval with dense output differs at chunk ends from the solver's last step → **Mitigation:** the last t_eval time equals the chunk stop exactly (solve_ivp then returns the exact end state); the test checks the next chunk's y0 equals the full-resolution run's chunk end to solver tolerance.
- **Risk:** the 10k-pulse benchmark takes too long to run routinely → **Mitigation:** benchmark has a `--pulses` argument and is outside pytest; the pytest scaling check is structural (no copy of earlier arrays), not timing-based.

## Assumptions

- SciPy solve_ivp with `t_eval` keeps the solver steps unchanged and interpolates with dense output (LSODA, BDF, Radau).
- Log-spaced times: chunk start (after the initial-burst shift), then N-1 times geomspace(span·1e-6, span) offset from the start, last time = chunk stop exactly; N >= 2 (validated, ValueError otherwise).
- The thinned flag is in memory only, like short_on_spans.
