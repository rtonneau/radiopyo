# Session: overshoot-path-precompute

**Date:** 2026-11-11T23:31:37.726Z
**Status:** Grill phase complete

## Problem Statement

In a pulsed chunked run (1e4 Gy/s peak, 10 ms period, 10 µs ON, 10 Gy, 200 chunks: 7.5 s), the overshoot relaxation path (ADR 0007) takes 38% of the time (`_k_relaxation_terms` 2.85 s). With a species below 0, every `derive` / `jacobian` call rebuilds the same per-species sub-matrices, and the derivative helpers are Python loops.

## Context & Constraints

- **Current behavior:** `_overshoot` loops over the species below 0; for each, `_k_relaxation_terms` / `_ez_relaxation_terms` rebuild sub-`KMatrices` / `EzMatrices` with `dataclasses.replace` (re-running `_lock_arrays`: 23,640 calls) and `_relaxation_sum` builds one more to drop the species' own negative-product row. `depletion_switch`, `depletion_switch_derivatives` (`np.delete`, 22,872 calls) and `mass_action_derivatives` (`np.delete` per slot) loop in Python; they also run on every `jacobian` call.
- **Pain point:** the tables depend only on the matrices and the species index, not on the state, yet they are rebuilt at every solver call.
- **Dependencies:** `SimMatrices` stays frozen and its constructor unchanged (built only in `environment_build.py`; not persisted to HDF5). `derive`, `jacobian`, `overshoot_slopes`, `reaction_rates` keep their signatures and results. ruff and mypy stay at zero errors.
- **Tech stack:** numpy, scipy `solve_ivp` (LSODA), frozen dataclasses in `utils/sim_types.py`.

## Success Metrics

- `derive` and `jacobian` match the current implementation (a reference copy kept in a test) within rtol 1e-12 on random states with species below 0, with and without enzymatic reactions and negative products.
- The golden, FD-Jacobian, overshoot relaxation and Fenton stall tests pass unchanged; ruff and mypy report zero errors.
- `examples/bench_overshoot_pulsed.py` (the scout scenario) reports before/after wall time and the overshoot share of the profile; the numbers are written in the commit log (no pytest threshold).

## Architecture & Approach

- A frozen, read-only `OvershootTables` (in `utils/sim_types.py`) held by `SimMatrices` as an `init=False` field built in `SimMatrices.__post_init__`, for every dynamic species: the k and enzymatic reaction columns that change it, the sub-arrays of those columns, the reactant-slot replacement mask and the "held once" flags for the mass-action d law / d c_j, and the `negative_products` of the sub-reactions without the species' own row.
- `_k_relaxation_terms`, `_ez_relaxation_terms` and `_relaxation_sum` read the tables instead of calling `dataclasses.replace`.
- Vectorise `depletion_switch`, `depletion_switch_derivatives` and `mass_action_derivatives`: no `np.delete`, no per-column Python loop.

## Assumptions & Trade-offs

- Tables are built for the dynamic species only (only they can be below 0 in `derive`); the build cost is paid once per `SimMatrices` (every rebuild: temperature, G values) even for runs that never go negative, assumed negligible.
- Vectorised products may change the last bits, hence rtol 1e-12 rather than bit-identical.
- Not doing: changing the overshoot math (ADR 0007), the rate laws or the solver options.

## Open Questions

None.

## Notes

From scout report `architecture-review-2026-11-11T13-37-44-572Z.html` (seed Strong). Size: Standard. No glossary term or ADR (internal, reversible refactor).
