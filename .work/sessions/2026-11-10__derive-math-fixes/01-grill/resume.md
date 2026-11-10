# Session: derive-math-fixes

**Date:** 2026-11-10T17:30:49.497Z
**Status:** Grill phase complete

## Problem Statement

A mathematical review of `derive` (`radiopyo/simulation/ode_functions_solver.py`) found that its core (mass action, Michaelis-Menten, radiolysis, units, n·k rule, Jacobian) is correct, but its edges are not: an overshoot below 0 can drive a species further negative or run a reaction backwards; `force_positive` patches the net derivative and breaks mass balance; the negative-product block is a hard step that makes the right-hand side discontinuous; the bundled catalase constant is likely 2x too large; a polynomial k on an `n A` reaction is silently not divided by n.

## Context & Constraints

- **Current behavior:** `reaction_rates` evaluates rate laws on the raw state (negative values included); a k or enzymatic reaction is forced to 0 while a negative product is at or below `DEPLETION_LIMIT` (1e-12 mol/L); `force_positive=True` zeroes the negative derivative of species below the limit after the stoichiometry is applied; the analytic Jacobian is piecewise (ADR 0003). The catalase reaction `2 H2O2 -- catalase >> O2 + 2H2O` has `k_value = 1.324e8` (changed from 6.62e7 in 2023-08, commit daef717) and Michaelis-Menten uses it as given with stoichiometry -2, i.e. d[H2O2]/dt = -4·6.62e7·E·S/(Km+S).
- **Pain point:** non-positivity-preserving rate laws (2A -> B with A < 0 diverges), mass balance broken by `force_positive`, a discontinuous RHS for the stiff solvers, a probable factor-2 error on catalase, and a silent convention trap for polynomial k.
- **Dependencies:** `derive` and `ODEResult.reconstruct_reactions` must keep sharing `reaction_rates`; the analytic Jacobian must keep matching `derive`; the old `tests/golden/data/golden_results.h5` must keep loading; ruff / mypy / pytest stay at zero.
- **Tech stack:** Python >= 3.11, numpy, scipy `solve_ivp` (LSODA default, analytic Jacobian for LSODA/BDF/Radau), pytest golden tests (`tests/golden/record_golden.py --only <part>`, parts derive / runs / beams / hdf, with `PYTHONHASHSEED` fixed).

## Success Metrics

- Every rate law (mass action and Michaelis-Menten) reads clipped concentrations max(c, 0); dy/dt = S·r holds exactly (mass-balance test), and an overshooting `2 A -> B` stays bounded (overshoot test).
- `force_positive` and `zero_depleted_consumption` no longer exist in the API (`RunOptions`, `run`, `prepare_*`, `derive`, `jacobian`); passing it raises TypeError; CHANGELOG has a "Removed" entry.
- The depletion block is replaced by the smooth depletion switch c⁺/(c⁺+ε), ε = DEPLETION_LIMIT = 1e-12 mol/L, for k and enzymatic reactions; the analytic Jacobian matches finite differences with the switch and with clipped species, including near ε.
- The bundled catalase k_value is 6.62e7.
- A non-scalar k on a single-reactant `n > 1` reaction logs one warning at build (caplog test); a scalar k does not.
- Each result-changing fix (clipping, switch, catalase) is its own ticket with before/after numbers; only the failing golden parts are re-recorded, after approval; pytest, ruff and mypy report zero errors.

## Architecture & Approach

1. Clipping (ADR 0004): `mass_action_rates`, `michaelis_menten_rates` and their derivatives read `np.maximum(y_full, 0)`; Jacobian: a species below 0 contributes 0, at exactly 0 the derivative of the c > 0 side is used (states starting at 0 keep an informative Jacobian).
2. Remove `force_positive` (ADR 0004): drop the parameter from `RunOptions`, `UnitCell.run` / `prepare_chunked_run` / `prepare_single_pulse`, `derive`, `jacobian`; delete `zero_depleted_consumption`; update tests and docs.
3. Smooth depletion switch (ADR 0005, supersedes the piecewise part of ADR 0003): replace `negative_product_block` by a factor = product over negative products of c⁺/(c⁺+ε), applied to k and enzymatic rates in `reaction_rates`; the Jacobian adds the switch derivative ε/(c+ε)² by the product rule. `DEPLETION_LIMIT` keeps its name and becomes ε. The `<=` vs `<` mismatch disappears with the block and `zero_depleted_consumption`.
4. Catalase (ADR 0006): Michaelis-Menten keeps its k as given (no n·k division); `reactions_biology.toml` k_value 1.324e8 -> 6.62e7 (README example already 6.62e7).
5. Polynomial-k warning: `logger.warning` once per reaction at build when a non-scalar k (polynomial, Arrhenius, equilibrium) sits on a single-reactant reaction with stoichiometry n > 1, saying it is used as given (ADR 0001).
6. Docs: README, CLAUDE.md solver section, CHANGELOG fragment; glossary already updated (Clipped concentration, Depletion switch).

## Assumptions & Trade-offs

- 6.62e7 is taken as the per-reaction k of Labarbe's Rb10 (rate halves); not verified against Labarbe's MatLab model.
- Clipping is always on, with no option to restore the raw behaviour; results only move where a step overshoots below 0.
- ε = 1e-12 mol/L is fixed (not configurable); the switch slows a reaction by under 1 % once its negative product is above 1e-10 mol/L.
- Hard removal of `force_positive` is a breaking API change, accepted; no deprecation period.
- Michaelis-Menten and k-reactions keep different k conventions (ADR 0006) rather than unifying.
- Not doing: changing the n·k rule for polynomial k (ADR 0001 stays), making ε configurable, changing the initial burst or G-values.

## Open Questions

- Labarbe's intended catalase H2O2 consumption rate (Rb10) is unchecked: if 6.62e7 turns out to be the H2O2 consumption rate itself, the value must be halved again (3.31e7).

## Notes

- Review findings: positivity invariance of exact mass action only holds when every consumed species is in its own rate law; negative products break it, hence the switch.
- The default `atol` (1e-12 µmol/L = 1e-18 mol/L) is six orders of magnitude below ε, so the switch is resolved by the solver.
- The HDF golden test pins rates rebuilt from the old json store (doubled catalase k 2.648e8 stored there); the store keeps loading and is not rewritten.
- New ADRs: 0004 (clipped concentrations, force_positive removed), 0005 (smooth depletion switch), 0006 (MM k as given).
