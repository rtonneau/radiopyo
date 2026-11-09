# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`radiopyo` is a pure-Python package (no C++, despite the global Geant4-DNA rules — those don't apply here) that simulates water radiolysis and radiobiology chemistry under UHDR/CONV irradiation by solving a system of ODEs with `scipy.integrate.solve_ivp` (LSODA by default). Originally a port of R. Labarbe's MatLab model. The math behind the ODE system is documented in `doc/20230425_RPM_extended_doc.pdf`.

## Commands

Development uses the Anaconda Python on PATH (`radiobio` conda env in the hatch scripts). Install editable with `pip install -e ./`.

```bash
python -m pytest                                   # all tests (pytest uses --import-mode=importlib)
python -m pytest tests/chemistry/test_species.py   # one file
python -m pytest tests/chemistry/test_species.py::test_name
ruff check radiopyo                                # lint (tests/ excluded)
mypy radiopyo                                      # type check (disallow_untyped_defs; tests/ excluded)
python -m build                                    # sdist + wheel into dist/ (hatchling backend; `hatch` CLI not installed)
```

Current state (2026-11): the test suite fails at collection — `tests/chemistry/conftest.py` imports `RawSpecies`/`ChemSpecies`, which no longer exist in `radiopyo/chemistry/species.py` (now `BasicSpecies`, `ChemicalSpecies`, `SimSpecies`, …). `ruff` reports ~33 errors and `mypy` ~44 on the package. Don't assume a green baseline; check before/after your change.

## Architecture

The pipeline is **config file → `SimEnv` → matrices → `derive()` → `ODEResult`**.

1. **Parsing** (`radiopyo/parser/`): `FileParser` subclasses (`toml.py`) produce a `ConfigDict` (TypedDicts in `radiopyo/utils/sim_types.py`). `resolve_includes` recursively loads the `[includes]` files first; later definitions override included ones (that is how a user tweaks a rate constant from a bundled file). An include path starting with `radiopyo/` refers to `radiopyo/data/` (`PACKAGE_DATA_PATH`). Reaction strings (`"A + B -> 2 C"`, `"2H2O2 -- catalase >> O2 + 2 H2O"`) are parsed by Lark grammars in `parser/grammar/` via the shared `RPARSER` singleton into `GenericReaction` dicts (`chemistry/generics.py`).

2. **Chemistry objects** (`radiopyo/chemistry/`): all reactions derive from `BasicReaction`. `KReaction` (plus `Arrhenius`, `CstReaction`) for rate-constant reactions; `BiReaction` subclasses `AcidBase` and `EReaction` are forward/backward pairs that expand into two k-reactions; `MichaelisMenten` for enzymatic reactions; `RadiolyticReaction` for G-value production. Species live in `ListSpecies` and are split into **dynamic** (solved) and **constant** (`[concentrations.fixed]`) species.

3. **Environment** (`simulation/environment_build.py`): `SimEnv.from_dict` builds a `ReactionsList` whose `__post_init__` assigns **global reaction indices in a fixed order**: acid/base (2 each), equilibrium (2 each), k-reactions, enzymatic, then radiolytic. Results analysis relies on these indices (`get_reaction_by_index`), so keep that ordering when adding reaction types. `build_all_matrices()` returns `SimMatrices`: `Ki_mat` (k-reaction rates/stoichiometry tensors), `Gi_mat` (radiolytic vector), `Ez_mat` (enzymatic), `cst_sp` (constant-species vector).

4. **Solver** (`simulation/ode_functions_solver.py`): `derive()` is a fully vectorized right-hand side (einsum/matmul over the matrices, no per-reaction loops). The state vector `y` is in **µmol/L**; `derive` converts to mol/L, appends constant species at the end of the vector, and strips them before returning. A reaction whose product is marked "negative" is zeroed when that species is depleted. `jacobian()` exists but is stale (its signature doesn't match the `derive` args) and is effectively unused.

5. **Driver** (`simulation/unit_cell.py`): `UnitCell` holds a `SimEnv`, the matrices, temperature, a `BeamCollection`, and `GValuesBoscolo`. Setting a beam with an LET rewrites the G-values of e_aq, OH_r, H_r, H2, H2O2 from the Boscolo 2020 fit (`utils/g_values.py`, valid LET 0.13–150 keV/µm) and rebuilds the matrices. On beam ON, the first `T_ZERO` = 1e-7 s is skipped and the radiolytic yield for that interval is added analytically to `y0`. `prepare_chunked_run()` splits pulsed beams into ON/OFF `SimChunk`s with different max step sizes; `run_till_max_dose` builds on `run`.

6. **Beams** (`physics/beam.py`): `ConstantBeam`, `SinglePulseBeam`, `PulsedBeam`, `MacroPulsedBeam`, driven by `Timer`/`MacroTimer` state machines. `beam.at(t).dose_rate()` is what `derive` calls on every evaluation.

7. **Results** (`simulation/ode_result.py`, `simulation/utils.py`): `ODEResult` wraps the scipy solution with the env and matrices so it can recompute per-reaction rates (production/consumption integrals, equilibrium yields, `integrate_species`, `to_pandas`). `simulation/utils.py` also holds unit converters (`O2Convert`, …), `ResCollection`, and HDF5 persistence: `H5Manager`, `SimIdentifier` (a dataclass key serialized into HDF group names, rebuildable from its `repr` via `ast`), and `SimCollection` (lazy, on-demand loading of stored `ODEResult`s). This uses `h5py`.

Logging: importing `radiopyo` configures the `radiopyo` logger from `[tool.logging]` in the repo's `pyproject.toml` (falls back to `utils/default_log_config.toml` when installed). Use `logging.getLogger("radiopyo")` and avoid `print` (ruff `T20`).

## Conventions

- Every module follows the same banner layout: shebang, docstring with `Author:` line, then `External imports` / `Internal imports` / `Type checking ONLY imports` / `Expose importable stuffs` (`__all__`) / `Logging` / `FUNCTION DEFINITIONS` / `CLASS DEFINITIONS` sections. Keep it when adding files.
- `import typing as tp`, `from __future__ import annotations`; fully typed signatures (mypy `disallow_untyped_defs`). Code uses `tp.Self`/`tp.Unpack`, hence `requires-python >= 3.11`.
- Species labels use `_minus` / `_plus` suffixes for charge and `_r` for radicals (e.g. `O2_r_minus`, `Fe_plus_plus`, `e_aq`); config concentrations are mol/L.
- Bundled configs: `radiopyo/data/advanced_water_radiolysis.toml` (`PURE_WATER_CONFIG_FILE`), `reactions_biology.toml`, `config_Labarbe.toml` (`LABARBE_CONFIG_FILE`). The README documents every TOML section and beam type.
- `.work/` is the gps (grill-plan-ship) session state and is committed. `origin` is GitHub (`rtonneau/radiopyo`); the original GitLab remote is kept as `gitlab`.
