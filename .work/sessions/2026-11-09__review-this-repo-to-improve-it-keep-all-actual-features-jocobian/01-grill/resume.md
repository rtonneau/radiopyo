# Session: review this repo to improve it. Keep all actual features. Jocobian will be one of the first next job. Just improve flow and readability. Don't break how code is used by user. Other things can be adapted if needed

**Date:** 2026-11-09T20:05:46.317Z
**Status:** Grill phase complete

## Problem Statement

Review radiopyo and improve its flow and readability while keeping every existing feature and not breaking how users call the code. The Jacobian rewrite is the next job, so the solver core must be restructured to make it easy. The review found a codebase with heavy duplication (rate maths written twice, 5 copies of run kwargs, 7 near-identical config builders, 6 copies of the include-merge block, 4 hand-written beam `make_matrix`), positional matrix tuples, a god-class `SimEnv`, a 2000-line `beam.py`, a test suite that fails at collection, and a set of real bugs (some crash/hang, some silently change results).

## Context & Constraints

- **Current behavior:** config file → `SimEnv` → matrices (`Ki_mat`/`Gi_mat`/`Ez_mat`/`cst_sp` positional tuples) → `derive()` → `ODEResult`. Tests fail at collection (`tests/chemistry/conftest.py` imports `RawSpecies`/`ChemSpecies`); ruff 28 errors, mypy 44 errors at 2261b92. Confirmed bugs: `UnitCell.run()` passes `force_positive` into `_raw_run`'s `beam_ON` slot (`unit_cell.py:433-444` vs `:294-306`), so a plain `run()` skips the initial burst; `PulsedBeam.from_dose_rate_peak`'s 4th positional argument is `n_pulse` while README calls it as `max_dose`.
- **Pain point:** hard to read and extend (especially for the Jacobian), no safety net, bugs hidden in duplicated code.
- **Dependencies:** must keep: the `rp.*` exports; `UnitCell.from_toml` / `from_parser` / `from_config_dict` / `run` / `set_beam(beam, LET=)` / `change_LET` / `set_temperature` / `update_env` / `prepare_chunked_run(...).run()` / `prepare_single_pulse` / `run_till_max_dose`; `uc.env.reactions.iter_all()`, `r.as_label()`, the `*_involving*` queries, `uc.env.species.get(...).set_initial_cc`; `ODEResult` (`[]`, `time`, `to_pandas`, `integrate_species`, `reconstruct_reactions`, integrate/detail/eq-yield methods); beam constructors with today's positional argument order; the TOML format (`[includes]`, reaction strings, all sections); deep import paths (kept through re-exports); HDF5 files written by the current version must still load (H5Manager / SimIdentifier / SimCollection format unchanged). Removed by decision (breaking, recorded in CHANGELOG): the RON parser and `UnitCell.from_ron`, `run_multi_doses`, `utils/e_LET.py` and `data/energy_LET/`.
- **Tech stack:** Python ≥3.11, numpy, scipy `solve_ivp` (LSODA), pandas, lark, h5py, hatchling; pytest, ruff, mypy (`disallow_untyped_defs`).

## Success Metrics

- Golden outputs recorded from the code at 2261b92 (derive at fixed states, `run()`, chunked pulsed run, `reconstruct_reactions`, beam matrices, an HDF5 file written by the old code) match after every Phase B ticket (bitwise where possible, otherwise rtol 1e-12).
- `python -m pytest` passes, including the ported `tests/chemistry` suite, include-merge tests, timer tests and the HDF5 old-file read test.
- `ruff check radiopyo` and `mypy radiopyo` report zero errors at the end.
- Every Phase C fix lands as its own ticket with before/after numbers shown to and approved by the user; golden outputs are regenerated only in those tickets.
- README examples run as written against the final code.

## Architecture & Approach

**Phase A – safety net.** Record golden outputs from the current code before any change. Port `tests/chemistry` to the current API (`BasicSpecies`/`SimSpecies`, `stoi=`, `from_product`, species objects). Add tests for the include merge (current behaviour) and the timers.

**Phase B – behaviour-preserving refactor** (no numeric change; crash/hang fixes allowed):
1. Hygiene: move mypy/pytest/pandas-stubs to a `dev` extra and tqdm to an optional extra; untrack `dist/`, `.coverage`, `tests/output/` and gitignore them; exclude `radiopyo/data/IBA/` from the wheel (files stay in place); move ruff settings to `[tool.ruff.lint]`; fix `__all__` (duplicates, missing names), export `MacroPulsedBeam`; logging: keep INFO console output but drop `disable_existing_loggers`, one config source; remove dead code (RON parser + `from_ron`, `run_multi_doses`, `e_LET`, `radiopyo/test.py`, `make_beam_func`, `RawReader`, `_average_val_before_100ns`, `reconstruct_jacobian`, `push_result`/`read_result`, commented-out code, stray `print`).
2. Solver core (Jacobian prep): named frozen matrix dataclasses replacing the positional tuples; one pure `reaction_rates()` used by both `derive` and a vectorised `reconstruct_reactions`; µM↔M conversion and fixed-species append in one place; remove the stale `jacobian` body and the 4-D `JRki` tensor, `use_jac=True` raises a clear `NotImplementedError` until the next job.
3. UnitCell run flow: one `RunOptions` object behind the existing kwargs; `_raw_run` split into small steps; deduplicate `prepare_chunked_run`/`prepare_single_pulse`; reproduce today's `run()` behaviour explicitly so goldens still match.
4. SimEnv and simulation utils: table-driven config readers; one `reactions_involving(...)` with old names as aliases; split `simulation/utils.py` into units / collections / persistence with re-exports.
5. Chemistry: shared reactant/product population in `from_generic`; `BiReaction` owns its two k-reactions with their k set; dedupe `AcidBase.get_couple`; `KReactionCollection` index rebuilt after relabel/sort; `None` instead of the `-1` initial-cc sentinel; proper `__eq__`/`__hash__`.
6. Beams: split `physics/beam.py` into timer / beam / beam_collection with re-exports; `BEAM_TYPES` registry; non-mutating `from_dict`; deduplicated timer methods; one generic `make_matrix`.

**Phase C – result-changing fixes**, one ticket each with before/after numbers for approval:
- `run()`: keyword argument passing (burst applied), no in-place `y0` mutation, no burst when `use_beam=False`, raise the missing `NoBeamDefinedError`.
- Temperature threaded to fixed species and initial values.
- Include merge rewritten as one generic merge on an order-insensitive reaction key (no parser→simulation import): fixes enzymatic k doubling, acid/base override, constant-reaction override; later includes override earlier ones, the including file overrides all; relative paths resolve against the including file (fallback: CWD); only the exact `radiopyo/` prefix maps to package data.
- Macro beam / timers: walrus precedence (pulses per train), hang without `max_dose`, `stop` ignored, macro start time, `tot_pulses` with `start_time`, `copy` keeping `start_time`.
- Deterministic species order; `ODEResult.__add__` returns a new object; Arrhenius with TOML keys `A` and `Ea` (kJ/mol), k = A·exp(−Ea/(R·T)); capitalised enzyme names parsed correctly.

**Phase D:** README fixes (version badge, `from_toml_file`, missing `basic_water_radiolysis.toml`, `max_dose=` keyword in examples, `macro_pulsed` docs, acid_base keys, deps, GitHub URL), CHANGELOG (breaking removals + fixes), ruff and mypy at zero.

## Assumptions & Trade-offs

- n·k convention stays scalar-only and is documented (ADR 0001); bundled configs only use scalar k on nA reactions.
- G(H3O⁺) = G(e_aq) is intended and kept (ADR 0002); the CSV's H3O_plus column is documented as unused.
- `state_at` keeps its `<=` boundary (a time exactly at pulse end is ON); only the docstring changes.
- Delivered as one session, one branch, one PR; one ticket per reviewable change, solver core first among the refactor tickets.
- Not doing: the Jacobian rewrite; physics re-validation of fits.

## Open Questions

- None blocking. Seeds for later: valid dose-rate range of the `GValuesBoscolo._scale_factor` fit (no clamp, inf at 0 Gy/s); sign check of the H2O initial-concentration polynomial (`species.py:363`, rises with T).

## Notes

- Review done by three parallel reviewers (simulation; chemistry/parser; physics/utils/packaging), key claims spot-checked by hand.
- Plausible, to verify during Phase B/C: `species_balance_evolution` counts only k-reactions; `Ez_iy` uses `=` instead of `+=`; `to_hdf` storing object lists as attributes; `parse_dataset_name` uses `.strip(".sim")` (careful with HDF5 compat).
- Crash-only fixes allowed in Phase B (Arrhenius is handled in Phase C because it fixes units too): `update_beam_param(period=)` KeyError, `from_dict` label KeyError, `get_radiolytic_reaction_from` UnboundLocalError, `SimCollection.__len__`, `ResCollection.sort`, empty-list unbound `idx`, `MichaelisMenten.__repr__` missing `)` (if hash-safe).
