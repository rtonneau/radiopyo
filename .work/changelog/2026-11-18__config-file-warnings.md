---
bump: minor
floor: minor
---
### Added

- Chemistry check (`radiopyo.simulation.chemistry_check`): `UnitCell.check_chemistry()` / `SimEnv.check_chemistry()` return a list of `ChemistryFinding` (`kind`, `species`, `reactions`, `message`) for species aliases (`OH` / `OH_r`, `O_minus` / `O_r_minus`, `e_minus` / `e_aq`, `H_plus` / `H3O_plus`), reactions whose charge is not balanced, dead reactants (needed but never present) and unused species. `SimEnv.from_dict` runs it and logs the findings (`CHEMISTRY::` warnings) once per distinct chemistry in a process; a manual call logs every time (`log=False` to only return them). The bundled `config_Labarbe.toml` gives two charge warnings (superoxide dismutation without its protons, `RH_r` as the radical anion of RH), kept on purpose.

### Changed

- A reaction whose scalar rate constant is 0 (`k_value`, `A`, or either constant of an acid/base or equilibrium pair, which removes the whole pair) is removed from the simulation with a warning (ADR 0011). Redefining an included reaction with `k_value = 0` is the way to switch it off. Results do not change, since such a reaction contributed nothing; species used only by removed reactions are no longer in the state vector.
- `H2O` in `[concentrations.fixed]` or `[concentrations.initial]` is ignored with a warning: it cannot replace the built-in H2O.
- A label of `[concentrations.initial]` that matches no species is reported by a warning instead of being ignored silently.
- An acid/base entry in the old syntax (`acid` / `base` keys) logs the converted entry to write before raising its error.
