# Changelog

All notable changes to this project. Versions follow [semantic versioning](https://semver.org/).

## [Unreleased]

### Fixed

- `UnitCell.run`: the arguments were forwarded to the solver by position, so `force_positive`
  never reached it and was used as the "beam ON" flag instead. The initial radiolytic burst
  (first `T_ZERO` seconds of a beam ON phase, added analytically to `y0`) is now applied by
  default (unless `use_beam=False`) and `force_positive` reaches `derive()`. Results of `run()`
  change: the trajectory starts at `t_span[0] + T_ZERO`, and `t_eval` times inside the skipped
  interval are dropped (with a warning). A beam ON span not longer than `T_ZERO` now raises
  `ValueError`.
- A `y0` given to `run` / `prepare_chunked_run` / `prepare_single_pulse` is never modified.
  In chunked runs the burst of a beam ON chunk no longer overwrites the last stored point of
  the previous chunk (it changed the reconstructed rates and production integrals).
- `run(use_beam=False)` / `prepare_*(use_beam=False)` add no initial burst. `ODEResult.integrate_production_of`
  accounts for the bursts actually applied (`ODEResult.burst_dose_rates`) instead of assuming one
  per beam pulse; results read from a file keep the former assumption
  (`burst_dose_rates` is not persisted).
- A missing beam now raises `NoBeamDefinedError` in `run` (the exception was created but not raised).
- Temperature now has a single source of truth, `env.parameters["temperature"]` (degC, default 25;
  `SimEnv.temperature` / `SimEnv.set_temperature`). `UnitCell.T`, `set_temperature` (and a new `T`
  setter) read and write it, and the `temperature` key of a config file is now honoured (it was
  ignored). The constant-species vector (e.g. H2O, a polynomial of the temperature in degC) was built
  without temperature, so it stayed at its 25 degC value whatever `set_temperature` said;
  `SimEnv.initial_values`, `ODEResult.to_pandas` and `ODEResult.__getitem__` also evaluated the
  species at the default temperature. All now use the simulation temperature. Results at 25 degC are
  unchanged; at another temperature the H2O concentration (hence the rates involving it) changes.
  Units: k polynomials take 1/K, pKa and initial concentrations take degC (conversion is internal).

### Removed (breaking)

- RON config format: `radiopyo.parser.ron`, its grammar and `UnitCell.from_ron`. Use TOML.
- `UnitCell.run_multi_doses` (was a non-functional stub raising `NotImplementedError`).
- `radiopyo.utils.e_LET` and the `data/energy_LET` data.
- `GValuesBoscolo._average_val_before_100ns`, `RawReader`, `umol_per_joule_to_ge`,
  `ODEResult.reconstruct_jacobian`, `H5Manager.push_result` / `read_result`,
  `KReaction.from_dict` / `Arrhenius.from_dict`, `Beam.make_beam_func`,
  `MichaelisMentenDict`, `parser.utils.dict_compare`, `radiopyo.test`, `radiopyo.physics.utils`.
