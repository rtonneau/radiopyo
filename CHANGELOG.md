# Changelog

All notable changes to this project. Versions follow [semantic versioning](https://semver.org/).

## [Unreleased]

### Added

- `radiopyo.exceptions.RadioPyoError`, the root of the exception hierarchy (the `chemistry`, `physics`
  and `simulation` exceptions derive from it; `radiopyo.simulation.exceptions.RadioPyoError` still works).
- `RunOptions` (`radiopyo.simulation.unit_cell`): the solver options (`atol`, `rtol`, `method`,
  `max_step`, `first_step`, `use_beam`, `force_positive`, `use_jac`, `t_eval`) behind the existing keyword
  arguments of `run`, `prepare_chunked_run` and `prepare_single_pulse` (their signatures are unchanged).
- `SimEnv.reactions_involving(species, role=, kinds=)`: one query behind the ten historical
  `*_involving*` methods, which stay as aliases. A `SimSpecies` is accepted as species.
- `UnitCell.T` setter, `SimEnv.temperature` / `SimEnv.set_temperature`.
- `MacroPulsedBeam` exported from `radiopyo`; `dev` (pytest, mypy, pandas-stubs) and `tqdm` extras in
  `pyproject.toml`; test suite (`tests/`, 551 tests) with golden outputs recorded from the code before the
  refactoring (`tests/golden/`).
- README: `[[beam.macro_pulsed]]`, `[[reactions.equilibrium]]`, the `[parameters]` keys, the n*k convention,
  `variable_GValues` and H3O_plus = e_aq (G(H3O+) is taken equal to G(e_aq), the `H3O_plus` column of the
  Boscolo data is not used).

### Changed

Internal refactoring: results are unchanged (bit-identical on the golden outputs) unless listed under
"Fixed"; the public names and call signatures are kept.

- Solver matrices: the positional tuples (`Ki_mat`, `Gi_mat`, `Ez_mat`, `cst_sp`) are replaced by the frozen
  dataclass `SimMatrices` (`radiopyo.utils.sim_types`) made of `KMatrices`, `GVector`, `EzMatrices` and
  `StateLayout` (read-only arrays). `derive(t, y, beam, mat, force_positive)` and
  `ODEResult.reconstruct_reactions` share one pure `reaction_rates()`; the depletion threshold is the single
  constant `DEPLETION_LIMIT`. `reconstruct_reactions` is vectorised over the time points. The conversion
  umol/L <-> mol/L and the fixed-species handling live in `StateLayout`.
- `UnitCell`: `_raw_run` is a sequence of small steps (`widen_empty_span`, `_solver_kwargs`,
  `_apply_initial_burst`, `_timed_solve`); `prepare_chunked_run` and `prepare_single_pulse` share their
  chunk-building code. The docstrings now give the real defaults (atol 1e-12, rtol 1e-10, first step 1e-15).
- `SimEnv`: table-driven readers for the reaction sections (`reactions_from_dict`, `_REACTION_BUILDERS`);
  `acid_base_from_dict` & co. are kept as thin wrappers. `pH = 0` is honoured. A config without the
  `concentrations` / `parameters` section raises a `ValueError` naming it.
- `radiopyo.simulation.utils` is split into `units.py` (unit converters), `result_collection.py`
  (`ResCollection`) and `persistence.py` (`H5Manager`, `SimIdentifier`, `SimCollection`); `utils.py` re-exports
  everything, so the old import paths keep working. The HDF5 format is unchanged (files written before
  still load). `ODEResult.to_hdf` delegates to `H5Manager.push_ODEResult`.
- Chemistry: the reactant/product population is shared (`BasicReaction._add_species_from_generic`);
  `BiReaction` owns its two k-reactions and their k values (`AcidBase.get_couple` deduplicated);
  `KReactionCollection` re-indexes after sorting and `update_label` renames in place;
  `SimSpecies.cc_init` is `None` (not a `Polynomial(-inf)`) without initial concentration;
  `BasicSpecies` has `__eq__` / `__hash__`; `RadiolyticReaction.label` is a property (the old
  `label()` call still works); `MichaelisMenten.__repr__` has its closing parenthesis (`make_hash` unchanged).
- Beams: `radiopyo.physics.beam` is split into `timer.py` (`Timer`, `MacroTimer`), `beam.py` (beams and the
  `BEAM_TYPES` registry) and `beam_collection.py` (`BeamCollection`); `beam.py` re-exports the others.
  `from_dict` of the beams does not modify the dictionary it is given. `PulsedBeam.update_beam_param(period=)`
  and a bad unlabelled beam config no longer raise `KeyError`.
- Packaging and logging: `mypy`, `pytest` and `pandas-stubs` moved to the `dev` extra, `tqdm` to the `tqdm`
  extra, `radiopyo/data/IBA` left out of the wheel, `dist/`, `.coverage` and `tests/output/` untracked,
  ruff settings under `[tool.ruff.lint]`, `__all__` fixed. Importing `radiopyo` configures its logger from
  `utils/default_log_config.toml` only, without disabling the loggers of the application
  (`disable_existing_loggers = false`); `[tool.logging]` is gone from `pyproject.toml`.
- `ruff check radiopyo` and `mypy radiopyo` report no error.
- README brought in line with the code (version badge, `from_toml`, `max_dose=` keyword in
  `PulsedBeam.from_dose_rate_peak`, bundled files, acid/base keys, dependencies, GitHub URL).

### Known issues (not fixed in this release)

- The analytic Jacobian is not available: `use_jac=True` raises `NotImplementedError` (rewrite planned on top
  of `reaction_rates` / `SimMatrices`). Enzymatic reactions ignore the negative-product flag.
- `GValuesBoscolo._scale_factor` (dose-rate factor of `variable_GValues`) has no validated range and is
  infinite at 0 Gy/s.
- The H2O concentration polynomial (`radiopyo/chemistry/species.py`) rises with the temperature (55.93 M at
  25 degC) whereas the density of water falls (about 55.35 M at 25 degC): to be checked against its source.
- `Timer.time_till_next_OFF` returns NaN for continuous timers (infinite period).
- `SimEnv.make_hash` always raises `TypeError` and hashes the `repr` of a `frozenset`
  (depends on `PYTHONHASHSEED`).
- `parse_dataset_name` strips characters (`.strip(".sim")`) instead of a suffix; `burst_dose_rates` is not
  stored in HDF5 files; `ODEResult.reconstruct_reactions` uses the beam of the result even for a run made with
  `use_beam=False`.
- No golden output for a chunked run with a `MacroPulsedBeam`.

### Fixed

- `[[reactions.arrhenius]]` could not be read from a config file (the reader asked for a `k_value`
  instead of `A` and `Ea`), and `Arrhenius.k_value` used `Ea` as a temperature
  (`A * exp(-Ea / T)`). `Ea` is now the activation energy in kJ/mol:
  `k = A * exp(-Ea * 1000 / (R * T))`, T in K from the simulation temperature (documented in the README).
  Any former `Ea` (in K) must be multiplied by 8.314e-3.
- Enzyme names with capitals (`-- SOD >>`, `-- GPx >>`) were parsed as chemical formulas (`BasicSpecies`
  holding a dict). In the reaction grammar an enzyme is now only a name.
- `ODEResult.__add__` modified and returned its left operand. It now returns a new `ODEResult` and leaves
  both operands unchanged.
- `MichaelisMenten.to_dict` wrote `k_value` times the substrate stoichiometry, while the reader takes the value
  as is: results saved to HDF5 and read back had their catalase `k_value` doubled (1.324e8 -> 2.648e8), hence
  2x enzymatic rates when rebuilt from the file. Files written before this fix keep their doubled value.
- `ODEResult.species_balance_evolution` only summed the k-reactions: the radiolytic and enzymatic terms were
  missing (e.g. H2O2 at 2.6 ms: +3.5e-9 mol/L/s reported, -1.15e-5 real) and `start`/`stop` were ignored. It
  now applies the solver stoichiometry matrices to all the reaction rates and honours `start`/`stop`.
- Enzymatic stoichiometry matrix: a species on both sides of the reaction (`Aa -- Ee >> Aa + Bb`) only got its
  product stoichiometry (`=`); it now gets the net one (`+=`). Catalase and the bundled configs are unchanged.
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
  included value.
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
- The stale `jacobian()` of the solver (its signature did not match `derive`, it was unused), `EPSILON`
  and the 4-D `JRki` tensor. `UnitCell.run(use_jac=True)` (and the chunked runs) now raises
  `NotImplementedError` until the Jacobian is rewritten.
- `SimEnv.pKw`, `pKH2O`, `update_temperature` and `update_flag_positive` (use `set_temperature` and the
  `force_positive` parameter); `ParticleBeam.delivery`,
  `MacroTimer.within_train`, `micro_pulse_index`, `micro_pulses_left_in_train`, `GROWTH_FACTOR`
  (unused).
- `UnitCell.run_multi_doses` (was a non-functional stub raising `NotImplementedError`).
- `radiopyo.utils.e_LET` and the `data/energy_LET` data.
- `GValuesBoscolo._average_val_before_100ns`, `RawReader`, `umol_per_joule_to_ge`,
  `ODEResult.reconstruct_jacobian`, `H5Manager.push_result` / `read_result`,
  `KReaction.from_dict` / `Arrhenius.from_dict`, `Beam.make_beam_func`,
  `MichaelisMentenDict`, `parser.utils.dict_compare`, `radiopyo.test`, `radiopyo.physics.utils`.
