# Session: analytic-jacobian

**Date:** 2026-11-10T10:48:00.872Z
**Status:** Grill phase complete

## Problem Statement

radiopyo has no analytic Jacobian: `UnitCell.run(use_jac=True)` and the chunked runs raise `NotImplementedError` (the old 4-D `JRki` `jacobian()` was removed in the refactor), so LSODA/BDF/Radau always estimate it by finite differences (seed F1). Also, enzymatic (Michaelis-Menten) reactions ignore the depletion block of their negative products, unlike k-reactions (seed F9), so a Jacobian written now would be inconsistent with them.

## Context & Constraints

- **Current behavior:** `derive` -> `reaction_rates` (mass action of k-reactions, `negative_product_block` on k-reactions only, Michaelis-Menten, radiolytic) -> `species_rates`; optional `force_positive` zeroing; `RunOptions.use_jac` defaults to False and True raises `NotImplementedError` in `_raw_run`.
- **Pain point:** stiff solvers pay one finite-difference Jacobian (n_dyn extra `derive` calls) per Jacobian evaluation; `use_jac` is a dead option; enzymatic reactions with a depleted negative product keep running.
- **Dependencies:** keep `derive`, `reaction_rates`, `SimMatrices` and the public run API (signatures unchanged except the `use_jac` default); `tests/golden/data` is kept unchanged; global reaction index order unchanged; `ODEResult.reconstruct_reactions` must stay consistent with `derive` (shares `reaction_rates`).
- **Tech stack:** numpy, `scipy.integrate.solve_ivp` (LSODA default; LSODA needs a dense Jacobian), pytest, ruff, mypy.

## Success Metrics

- `jacobian(...)` matches a finite-difference Jacobian of `derive` (stated rtol/atol) on the bundled configs, with beam on and off, enzymatic reactions, fixed species, an active depletion block and `force_positive`.
- Enzymatic reactions are blocked while a negative product is depleted; no bundled or golden result changes from this.
- `use_jac=True` (new default) runs on LSODA, BDF and Radau; explicit methods get no `jac` and no warning; all golden tests pass with the new default against the unchanged recorded data within a tolerance written in the test.
- A benchmark (wall time, nfev, njev, success) of analytic vs finite-difference Jacobian on the bundled configs is saved under `.scratch/`; the default is True only if the analytic Jacobian is at least as fast and robust, otherwise stop and ask the user.
- `python -m pytest`, `ruff check radiopyo` and `mypy radiopyo` stay at zero failures/errors.

## Architecture & Approach

1. F9: `EzMatrices` gains `negative_products` (n_species, n_ez) bool built in `make_enzymatic_matrix` from the reaction's `ope`; `reaction_rates` applies the same depletion block to enzymatic rates as to k-reactions (shared helper).
2. New public pure function `jacobian(t, y, beam, mat, force_positive=False)` in `ode_functions_solver.py`, dense, built on `SimMatrices`: d(k-rate)/dy from `reactant_index` (product rule over reactant slots, repeated reactants handled), d(MM rate)/d[enzyme] and d/d[substrate], radiolytic and fixed species contribute nothing; J_full = S_k @ dr_k/dy + S_ez @ dr_ez/dy, blocked reactions contribute 0, rows zeroed by `force_positive` are 0 (piecewise, the current side of each switch, ADR 0003); returns the dynamic block (µmol/L in, µmol/L/s out: the 1e-6 / 1e6 factors cancel).
3. `UnitCell`: remove the `NotImplementedError`, pass `jac=jacobian` (same args as `derive`) to `solve_ivp` only for LSODA, BDF and Radau; `use_jac` defaults to True in `RunOptions`, `run`, `prepare_chunked_run`, `prepare_single_pulse` (after the benchmark gate).
4. Tests: FD comparison test, enzymatic block test, use_jac/method wiring test, goldens run with the default and compared within tolerance; benchmark script in `.scratch/`.
5. Docs: README, CHANGELOG (Unreleased), CLAUDE.md; ADR 0003; glossary terms negative product, depletion block, analytic Jacobian.

## Assumptions & Trade-offs

- Dense Jacobian: bundled systems are small (reactions_biology: 22 species, 15 dynamic, 12 k-reactions, 1 enzymatic) and LSODA only accepts dense; no sparse path.
- Piecewise Jacobian (ADR 0003): exact away from the switch thresholds, not at them; the smooth unblocked Jacobian was rejected because it disagrees with `derive` whenever a block is active.
- With the default True, default results move within solver tolerance; goldens are not re-recorded (only with approval if the tolerance check fails).
- Not doing: sparse Jacobian, d/dt terms (solve_ivp does not need them), changing the n·k rule, other seeds (H2O polynomial, G-value range, timers, persistence).

## Open Questions

- Whether the analytic Jacobian is actually faster/robust on the bundled configs: decided by the benchmark; the user is asked before keeping the default True if not.
- The golden tolerance value: chosen from the observed use_jac True vs recorded differences, written in the test.

## Notes

- Bundled configs: the only enzymatic reaction is catalase (`2 H2O2 -- catalase >> O2 + 2H2O`), no negative product, so the F9 fix changes no bundled result.
- Seed source: scout-reports/review-18-followups-2026-11-10T09-29-18-858Z.md (F1, F9).
- Decisions recorded: ADR 0003 (default True, piecewise Jacobian); glossary terms added.
