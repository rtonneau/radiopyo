# Session: chunked-run-scaling

**Date:** 2026-11-12T10:04:38.416Z
**Status:** Grill phase complete

## Problem Statement

Chunked runs (`prepare_chunked_run`, `prepare_single_pulse`) do not scale to conventional pulsed linacs (~1e4 pulses, ~2e4 chunks). Every `ODEResult.append_sol` copies the whole stored history (`np.append` / `np.hstack`), so the storage cost is quadratic in the number of chunks; on top of that every solver step is stored (~450 points per chunk on config_Labarbe), so memory grows to GB, and `derive` has a sizeable per-call NumPy overhead.

## Context & Constraints

- **Current behavior:** Measured on config_Labarbe (23 species), PulsedBeam 1e4 Gy/s peak, 10 µs ON, 10 ms period (scratch script `.scratch/tests/2026-11-12__chunked-run-scaling/scale.py`): 100 pulses / 200 chunks / 96 k points -> 23 s (117 ms per chunk); 400 / 800 / 365 k -> 95 s (118 ms); 1600 / 3200 / 1.43 M -> 455 s (142 ms), of which `hstack` + `append_sol` = 114 s (25 %). In derive, `to_full_molar`, `np.ones` and `broadcast_to` take about 25 % of the solve. `SimChunksIter.__next__` reads `self.sol.y[:, -1]` for the next y0 and yields the growing `ODEResult` on every chunk.
- **Pain point:** extrapolated to 10k pulses: ~70 min of quadratic copying on top of ~35 min of linear solve, and ~9 M stored points (~1.6 GB for `y` alone, several times that for `reconstruct_reactions`).
- **Dependencies:** results must stay bit-identical for the storage change and identical to rounding for the derive change (golden tests `tests/golden/`, `tests/simulation/test_overshoot_equivalence.py` rtol 1e-12); `ODEResult.t` / `.y` stay readable/assignable as ndarrays (persistence/HDF5 loading, `from_dataframe`, `__add__`, `to_pandas`, `reconstruct_reactions`); files written by older versions keep loading; ruff and mypy at zero.
- **Tech stack:** Python 3.11+, NumPy, SciPy `solve_ivp` (LSODA, `t_eval` uses dense-output interpolation), pandas, h5py.

## Success Metrics

- Appending a chunk no longer copies the earlier history: a fast pytest check shows the append cost / copies do not depend on the stored history length; results of existing tests and golden outputs unchanged.
- `examples/bench_chunked_scaling.py` (out of pytest) runs up to 10k pulses on config_Labarbe and shows a flat per-chunk cost from a few hundred to thousands of chunks.
- `points_per_chunk_on` / `points_per_chunk_off` on `prepare_chunked_run` and `prepare_single_pulse` store exactly the requested number of log-spaced points per chunk (first and last kept), the next chunk's y0 is the exact chunk end, and integral methods on a thinned result log one warning; default None leaves results unchanged.

## Architecture & Approach

1. **Linear storage (bit-identical):** `ODEResult` keeps per-step lists of `t` / `y` arrays; `t` and `y` become properties that concatenate lazily and cache the result, invalidated on append; setters keep plain assignment working (persistence, `from_dataframe`, `__add__`). `SimChunksIter` takes the next y0 from the last appended chunk's own array, never from the concatenated history.
2. **derive overhead (identical to rounding):** remove per-call allocations (`np.ones`, `broadcast_to`, `to_full_molar` work) in `derive`, `jacobian` and `_mass_action`, guarded by the golden and overshoot-equivalence tests.
3. **Opt-in point thinning:** new `points_per_chunk_on` / `points_per_chunk_off` on `prepare_chunked_run` and `prepare_single_pulse` (there the OFF value covers the transition chunks; the final chunk follows its beam state). Each chunk gets a `t_eval` of N times log-spaced from the chunk start (after the initial-burst shift), first and last points always kept; solver steps unchanged. Default None = every solver step stored (current behavior). `ODEResult` records in memory that it is thinned; the integral methods (trapezoid over stored points) log one warning on a thinned result; a short ON span does not warn about the generated t_eval.
4. **Acceptance:** `examples/bench_chunked_scaling.py` benchmark up to 10k pulses; fast pytest O(chunk) append check; README, CHANGELOG fragment and CLAUDE.md updated.

## Assumptions & Trade-offs

- Solver cost per chunk (LSODA steps) is not reduced beyond the derive overhead trim; it stays the linear floor.
- Thinning trades integral accuracy for memory; it is opt-in, documented and warned about. Exact integrals at full solver resolution are deferred to seed `exact-integrals-thinned-runs`.
- The thinned flag is in memory only (results read from a file carry no flag, no warning).
- Not doing: the deepcopy of env/beam/matrices in `ODEResult.__init__` (once per run), periodic steady-state skipping (result-changing), thinning for plain `run()` (it already has `t_eval`), a `points_per_chunk_final` override.

## Open Questions

- Lower bound of the log spacing (first offset after the chunk start): implementation choice, set in the plan (e.g. tied to `FIRST_STEP` / a fraction of the span).

## Notes

Scout seed: architecture-review-2026-11-11 ("chunked-run-scaling", Worth exploring). The per-chunk cost measured here (~117 ms on config_Labarbe, profiled) is higher than the 37 ms in the seed (other config/settings). Glossary terms added: Stored point, Point thinning. New seed recorded: `exact-integrals-thinned-runs`.
