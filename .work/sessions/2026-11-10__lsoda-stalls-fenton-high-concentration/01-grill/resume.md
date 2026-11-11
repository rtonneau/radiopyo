# Session: lsoda-stalls-fenton-high-concentration

**Date:** 2026-11-10T21:30:13.840Z
**Status:** Grill phase complete

## Problem Statement

With the Fenton reaction (`Fe_plus_plus + H2O2 -> Fe_plus_plus_plus + OH_minus + OH_r`, Fe2+ fixed at 2.2e-8 M) and RH = 1 M, O2 = 50 µM, pulsed chunked runs "freeze" (minutes per chunk), restart, freeze again. The same happened earlier with RH = 1 M without Fenton. Dose-per-pulse × dose-rate maps of ROOH show isolated "islands" of suddenly lower ROOH.

Reproduction (`.scratch/tests/2026-11-10__lsoda-stalls-fenton-high-concentration/`, `fe.toml` = config_Labarbe + Fe2+ 2.2e-8 fixed + O2 50e-6): `PulsedBeam.from_dose_per_pulse(pulse_dose=1, period=0.1, on_time=1e-4, n_pulse=10)`, LET 0.201, `prepare_chunked_run(t_span=(0, 3), max_step_size_on=on_time/10, max_step_size_off=0.01).run()` with defaults: chunk 6 runs > 5 min; the full run was at chunk 9/20 after 36 min. On a 4 × 3 grid (pulse dose 0.1–3 Gy, 1e3–1e7 Gy/s) 10/12 points exceed 180 s.

## Context & Constraints

- **Current behavior:** LSODA, atol 1e-12 µmol/L (1e-18 M), rtol 1e-10. During beam OFF it steps at ~1e-8 s with a Jacobian on ~40 % of evaluations while OH_r, e_aq, H_r, O_minus hover at ~-1e-13 µmol/L. Rate laws read max(c, 0) (ADR 0004): below 0 the slope of a stiff radical (OH_r + RH, 1e9 s⁻¹) is 0, above 0 it is −1e9; Newton fails at that kink. Changing atol only moves the stall (1e-10 fixes Fe 2.2e-8 but breaks Fe 8.9e-7). A second stall sits at the pulse end: the ON chunk extends GROWTH_TIME (1 µs) past the pulse, so the dose-rate jump to 0 is inside the integrated interval, and `beam.at(t_end)` already reads 0. Failed solves (BDF/Radau: "Required step size is less than spacing between numbers") are stored in `ODEResult.success` but the next chunk silently starts from the truncated point: ROOH 15–40× too low, no message.
- **Pain point:** runs take hours, and results can be silently wrong (likely the map islands).
- **Dependencies:** keep clipping (ADR 0004: no reaction runs backwards), the analytic Jacobian exact for `derive` (FD tests), reaction index order, old HDF5 files loading. Golden run data is re-recorded only after the user approves before/after differences.
- **Tech stack:** scipy `solve_ivp` (LSODA/BDF/Radau), numpy; `radiopyo/simulation/ode_functions_solver.py`, `unit_cell.py`, `physics/beam.py`.

## Success Metrics

- The reproduction (fe.toml, pulse dose ∈ {0.1, 1, 2} Gy × dose rate ∈ {1e3, 1e5, 1e7} Gy/s where on_time ≥ T_ZERO, and the bundled Fe 8.9e-7 variant) finishes every point with default options, each in under ~60 s, all `success` True.
- ROOH, H2O2, O2 at the end agree with BDF / LSODA atol 1e-10 references within 1e-6 relative where those succeed (prototype: 1e-7).
- A failed solve raises `SolverFailedError` (chunked runs and `run`), with chunk, reached t and the scipy message.
- `python -m pytest`, `ruff check radiopyo`, `mypy radiopyo` stay at zero errors; golden run data re-recorded after approval of the before/after table.

## Architecture & Approach

1. **Overshoot relaxation** (ADR 0007, amends ADR 0004): in `derive`, for a dynamic species with c < 0 add J_jj(0⁺)·c to its own row (J_jj(0⁺) = diagonal of the Jacobian evaluated at the clipped state); other species keep reading it as 0. `jacobian` includes the term exactly (diagonal entry J_jj(0⁺) for negative columns plus c·∂J_jj/∂y_k). Computed efficiently (only the negative columns), no extra full Jacobian per call.
2. **Chunks cut exactly at the pulse end:** `make_time_intervals` no longer extends ON intervals by GROWTH_TIME; the OFF chunk starts at the pulse end.
3. **Constant dose rate per chunk:** in a chunked run each chunk gets the dose rate of its interval held constant (no `beam.at(t)` jump at the chunk boundary); `run()` over a span with several pulses keeps `beam.at(t)`.
4. **Solver failure:** `_raw_run` (hence `run`, chunks, single pulse) raises `SolverFailedError(RadioPyoError)` in `simulation/exceptions.py` when `sol.success` is False.

## Assumptions & Trade-offs

- dy/dt is no longer exactly S·r when a species is below 0: the imbalance is the refill of an atol-sized hole.
- Dropping the GROWTH_TIME extension removes the refined max step for the first µs after a pulse; prototype results are identical.
- Not doing: a looser default atol, a smooth (softplus) clip, unclipped first-order rates (see ADR 0007), changing the default method.
- Out of scope (separate ticket): pulses shorter than T_ZERO (100 ns) deliver no dose (burst read at t0 + T_ZERO when the beam is already off); a partial last pulse ≤ T_ZERO raises "Too Short beam pulse".

## Open Questions

- Whether the user's map islands come only from silent failures or also from LSODA accuracy loss while stalling: verified after the fix by re-running a map slice.

## Notes

Prototype evidence (monkeypatched, `proto3.py --g0 --const`): 14/14 grid points finish in 2–58 s at default tolerances vs 2/12 before; without relaxation 4/14; constant dose rate without the exact cut gives ROOH +9 % (wrong: the 1 µs tail is dosed). User sets up maps with defaults and Fe2+ fixed.
