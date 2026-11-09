# Changelog

All notable changes to this project. Versions follow [semantic versioning](https://semver.org/).

## [Unreleased]

### Removed (breaking)

- RON config format: `radiopyo.parser.ron`, its grammar and `UnitCell.from_ron`. Use TOML.
- `UnitCell.run_multi_doses` (was a non-functional stub raising `NotImplementedError`).
- `radiopyo.utils.e_LET` and the `data/energy_LET` data.
- `GValuesBoscolo._average_val_before_100ns`, `RawReader`, `umol_per_joule_to_ge`,
  `ODEResult.reconstruct_jacobian`, `H5Manager.push_result` / `read_result`,
  `KReaction.from_dict` / `Arrhenius.from_dict`, `Beam.make_beam_func`,
  `MichaelisMentenDict`, `parser.utils.dict_compare`, `radiopyo.test`, `radiopyo.physics.utils`.
