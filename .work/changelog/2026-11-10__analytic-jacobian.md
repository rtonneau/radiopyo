---
bump: minor
floor: minor
---
### Added

- `jacobian(t, y, beam, mat, force_positive=False)` (`radiopyo.simulation.ode_functions_solver`): the analytic Jacobian of the ODE system, used by the stiff solvers (blocked reactions and `force_positive` rows contribute 0).
- `EzMatrices.negative_products`: negative-product flags of the enzymatic reactions.

### Changed

- `use_jac` now defaults to `True` (`RunOptions`, `run`, `prepare_chunked_run`, `prepare_single_pulse`): LSODA, BDF and Radau get the analytic Jacobian, other methods ignore it without a warning. On the bundled configs LSODA runs 1.03-2.12x faster (BDF 0.97-1.49x); results agree with `use_jac=False` to about 1e-9 relative. Pass `use_jac=False` for the former finite-difference behaviour.
- Enzymatic reactions are now blocked while one of their negative products is below `DEPLETION_LIMIT`, like k-reactions.

### Fixed

- A negative product of an enzymatic reaction was added in the stoichiometry instead of subtracted. No bundled config is affected.
