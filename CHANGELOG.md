# Changelog

All notable changes to this project. Versions follow [semantic versioning](https://semver.org/).

## [Unreleased]

### Fixed

- Species indices (and the summation order in the solver) depended on `PYTHONHASHSEED`:
  `ListSpecies.sort` numbered the non-starter dynamic species by iterating a `set`. It now keeps the
  insertion order (order of appearance in the config file), so indices, matrices and the solver step
  sequence are reproducible from one process to the next. Results change only at round-off level
  (max relative difference 1.3e-7 on R_r, the Labarbe run now always gives 2432 solver points instead
  of 2294 to 2448 depending on the seed).
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
- Config `[includes]` merge (`radiopyo.parser.utils.combine_config`, `resolve_includes`) rewritten as
  one generic merge keyed on an order-insensitive reaction key (`radiopyo.parser.utils.reaction_key`).
  It no longer builds reaction objects, so `radiopyo.parser` does not import `radiopyo.simulation`
  anymore, and it does not modify its arguments. Fixed:
  - overriding an enzymatic reaction doubled its `k_value` (the stoichiometry was applied again);
  - overriding an acid/base (or equilibrium) reaction with a single constant failed or dropped the
    other one: only the fields written are replaced now (an `equilibrium` giving 2 or more of its 3
    constants still replaces them all);
  - a `constant` reaction override was added as a duplicate (it was compared with the k-reactions);
  - `B + A -> C` was a different reaction from `A + B -> C`, now it overrides it;
  - circular includes ended in a `RecursionError`, they now raise a `ValueError` naming the cycle.
  Behaviour changes to be aware of:
  - between two included files, the LAST one listed now wins (it was the first); the including file
    still overrides all its includes;
  - a relative include path is resolved against the folder of the including file, then the working
    directory (it was the working directory only); only an entry starting with exactly `radiopyo/`
    refers to the package data (anything starting with `radiopyo`, like `radiopyo_extra.toml`, used
    to);
  - the reactions of a merged config follow the order of the includes (the bundled
    `config_Labarbe.toml` now gives `advanced_water_radiolysis` then `reactions_biology`, it was the
    reverse), so global reaction indices change for configs with includes. Results do not change
    beyond floating-point summation order.
  `TOMLFileParser.parse` now starts from an empty-parameters config, so a default can never shadow an
  included value. The `ODEResult` json round trip of an enzymatic reaction still doubles its `k_value`
  (`MichaelisMenten.to_dict`), this is unrelated to includes and not changed here.
- Beams (`radiopyo.physics`): one generic `ParticleBeam.make_matrix` walks the timer transitions up to
  `min(stop, time_at_max_dose)`; the four hand-written versions (constant, single pulse, pulsed, macro
  pulsed) are removed. `make_time_intervals` is shared and its second column is now always 0 (OFF) or 1
  (ON) (it was the raw dose rate for 2-row matrices and the first row was flagged ON when the window
  started before the beam). `PulsedBeam` matrices are unchanged. Behaviour changes:
  - `ConstantBeam` / `SinglePulseBeam` ON at `stop` (before the max dose) now end with `[stop, 0]`
    like `PulsedBeam` (they ended with `[stop, dose_rate]`), and no zero-length ON row is added when
    `stop` is the start time or the max-dose time;
  - `MacroPulsedBeam` without `max_dose` hung forever in `make_matrix`; it now respects `stop`. Micro
    pulses are no longer duplicated at train start, the matrix honours `start` / `stop` windows and
    `micro_start_time` (which was not given to the macro timer), and the max dose of whole trains is
    reached at the end of the last micro pulse (it was the end of the train);
  - `MacroPulsedBeam.update_beam_param` accepts `period`/`frequency`/`on_time`/`start_time` (macro) and
    `micro_period`/`micro_frequency`/`micro_on_time`/`micro_start_time`, and keeps `timer` and
    `macro_timer.subtimer` as one object (before, the micro timer only was rebuilt, with the macro
    names);
  - average dose rate and `tot_pulses()` account for `start_time` (a pulsed beam delayed by 1 ms
    reported an average dose rate 144 times too large), and `max_dose=inf` gives `peak * duty cycle` / `inf`
    instead of `nan`;
  - `MacroTimer`: `micro_pulses_per_train` and `duty_cycle` had an operator precedence error in the
    remainder test (train of 0.205 s, micro period 0.1 s, micro on_time 0.01 s: 2.5 pulses, not 3)
    and `duty_cycle` dropped the partial pulse; `time_till_next_ON` shifted `start_time` twice;
  - `Timer.copy` / `MacroTimer.copy` keep the start time (they reset it to 0);
  - counting pulses and trains (`MacroTimer`, `MacroPulsedBeam` max dose) and the `n_pulse` /
    `max_dose` consistency check of `PulsedBeam.from_dose_rate_peak` use a relative tolerance instead
    of exact float comparisons.
  - a `micro_start_time` >= the train `on_time` (no micro pulse fits in a train) now raises
    `ValueError` (`make_matrix` looped forever).
- `BeamCollection.merge` always raised `TypeError` when both collections hold a beam of the same name
  and type.

### Removed (breaking)

- RON config format: `radiopyo.parser.ron`, its grammar and `UnitCell.from_ron`. Use TOML.
- `UnitCell.run_multi_doses` (was a non-functional stub raising `NotImplementedError`).
- `radiopyo.utils.e_LET` and the `data/energy_LET` data.
- `GValuesBoscolo._average_val_before_100ns`, `RawReader`, `umol_per_joule_to_ge`,
  `ODEResult.reconstruct_jacobian`, `H5Manager.push_result` / `read_result`,
  `KReaction.from_dict` / `Arrhenius.from_dict`, `Beam.make_beam_func`,
  `MichaelisMentenDict`, `parser.utils.dict_compare`, `radiopyo.test`, `radiopyo.physics.utils`.
