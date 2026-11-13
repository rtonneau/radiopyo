# Session: ship rust radkin backend

**Date:** 2026-11-13T10:37:45.991Z
**Status:** Grill phase complete

## Problem Statement

radkin (`${DEV_DIR}/RUST/radkin`, GitHub `rtonneau/radkin`) is a Rust port of radiopyo that reproduces its results (golden tests at 1e-12 on derive/Jacobian, rtol 1e-6 on runs) with the LSODA port rsoda. Today it is only reachable as a Rust library or a CLI writing CSV files. Ship it with radiopyo as an optional, in-process second backend, callable from a JupyterLab notebook, that returns a regular radiopyo `ODEResult`, plus a parallel batch API for parameter sweeps.

## Context & Constraints

- **Current behavior:** radiopyo solves with `scipy.integrate.solve_ivp` (LSODA by default) through `UnitCell.run` / `run_till_max_dose` / `prepare_chunked_run(...).run()` / `prepare_single_pulse(...).run()`. radiopyo is a pure-Python hatchling package installed with `pip install git+https://github.com/rtonneau/radiopyo.git`. `UnitCell` does not keep its `ConfigDict` after construction. radkin offers `UnitCell::from_config/from_toml`, `run`, `run_till_max_dose`, `chunked_run`, `single_pulse_run`, and an `OdeResult` (t, y, burst dose rates, step dose rates, short ON spans, thinned, nfev/njev/nst, exec time) whose fields mirror `ODEResult.append_sol`.
- **Pain point:** no way to use the Rust solver from Python without going through the CLI and CSV files; parameter sweeps run serially in Python.
- **Dependencies:** radiopyo must stay pure Python and installable without Rust; the Rust model lives only in radkin (no vendored copy). radkin pins rsoda by git rev and that pattern is reused: the binding pins radkin by git rev. Changes needed in radkin go through a branch + PR in `rtonneau/radkin` (its CI: fmt, clippy -D warnings, test), merged by the user, then the rev is bumped here (no committed `path` dependency; a temporary `[patch]` for local testing is fine). Old HDF5 files and the golden data must keep loading/passing; ruff and mypy stay at zero errors.
- **Tech stack:** Python 3.11+ (conda env `radiobio`), numpy/scipy/pandas, hatchling; Rust 2024 stable (cargo 1.99 installed), PyO3 + maturin (to install), rayon, numpy (rust-numpy), pyo3-log; GitHub Actions (maturin-action for wheels).

## Success Metrics

- `uc.run(t_span, ..., backend="rust")` and the same option on `run_till_max_dose`, `prepare_chunked_run(...).run()` and `prepare_single_pulse(...).run()` (including point thinning and `t_eval` in run mode) return an `ODEResult` whose species agree with the Python backend at rtol 1e-6 on the bundled configs (pure water, Labarbe, paper config) × beam types × modes; `to_pandas`, `reconstruct_reactions`, integrals and HDF5 persistence work on it.
- Parity guard: before every Rust run, species and reaction labels equal (exact) and k / G vectors equal (rtol 1e-12) between radkin and Python's `SimMatrices`, else a clear error is raised; a test mutating a k in the Python env triggers it. A `UnitCell` built directly from a `SimEnv` (no config) raises a clear error with `backend="rust"`; unsupported options (`method` other than LSODA, etc.) raise.
- `radiopyo.rust.run_batch([uc, ...], t_span, mode=..., n_threads=None, return_exceptions=False, **options)` runs the cells in parallel (GIL released, rayon), returns results in input order, runs every variant to the end then raises the first error naming its index (or puts the exception in its slot with `return_exceptions=True`); a batch is faster than the serial loop on a multi-core machine.
- Without the extension installed, radiopyo imports and all non-Rust tests pass; `backend="rust"` raises an `ImportError` with install instructions; the Rust tests are skipped. CI builds the extension and runs the Rust tests; release CI builds Windows and Linux wheels attached to the GitHub release.
- `examples/notebook_rust_backend.ipynb` runs end to end (single run, batch LET sweep, timing comparison); README and `docs/configuration.md` document install and use; ruff and mypy report zero errors.

## Architecture & Approach

- **Packaging:** new `rust/` directory in radiopyo holding a PyO3 + maturin crate (`rust/Cargo.toml`, `rust/pyproject.toml`, `rust/src/lib.rs`), built as the optional wheel `radiopyo-rs` (abi3, Python >= 3.11), extension module imported as e.g. `radiopyo_rs._core`, depending on radkin via a pinned git rev. radiopyo gets an extra `radiopyo[rust]`; the wheel version is lockstep with radiopyo's (release `versionFiles` gains the rust manifests). Install: `pip install "radiopyo-rs @ git+https://github.com/rtonneau/radiopyo.git#subdirectory=rust"` (needs cargo) or prebuilt wheels from GitHub releases.
- **Boundary (ADR 0010):** Python sends the merged config dict (includes already resolved) plus the current temperature, beam (`beam.as_dict()` + its type name, read by radkin's `Beam::from_table`) and LET; radkin builds its own `UnitCell` and runs the whole simulation (bursts, chunking, short ON spans). It returns t, y, burst dose rates, step dose rates, short ON spans, thinned flag, solver counters, exec time, and its labels and k/G vectors for the parity guard. `UnitCell` stores the merged `ConfigDict` when built via `from_toml` / `from_config_dict` / `from_parser`.
- **Python side:** a module `radiopyo/rust.py` (or `radiopyo/simulation/rust_backend.py` re-exported as `radiopyo.rust`) does the lazy import, the parity guard, the conversion of options, and wraps the Rust output into a regular `ODEResult` (with the Python env/mat/beam), mapping radkin errors to radiopyo's exceptions (`SolverFailedError`, `NoBeamDefinedError`, `ShortBeamOnError`, ...). `UnitCell.run`, `run_till_max_dose` and `SimChunksIter.run` (chunked / single-pulse) take `backend: Literal["python", "rust"] = "python"`; LSODA only.
- **Batch:** `radiopyo.rust.run_batch(cells, t_span, mode="run"|"chunked"|"single-pulse", n_threads=None, return_exceptions=False, **options)`: one mode and one option set for all cells, parity guard per cell, Rust runs them on a rayon pool with the GIL released, results in input order; raise-after-all semantics.
- **Logging:** radkin `log` records forwarded to Python logging (pyo3-log) under the `radiopyo` logger.
- **radkin changes** (expected small: build a `Config` from an already-merged table, expose k/G vectors and labels, needed `pub` items): branch + PR in `rtonneau/radkin`, user merges, rev bumped here.
- **Tests/CI:** pytest comparisons Rust vs Python at rtol 1e-6 (configs × beams × modes), parity guard, error paths, batch order and failure semantics, all skipped without the extension; a CI job installs Rust + maturin, builds the extension and runs them; a release workflow builds abi3 wheels for Windows and Linux and attaches them to the GitHub release.
- **Docs:** `examples/notebook_rust_backend.ipynb`, README install/usage section, `docs/configuration.md`, CHANGELOG fragment.

## Assumptions & Trade-offs

- Two builds of the same system (Python and radkin) instead of one shared matrix set: accepted for reuse of radkin's golden-tested drivers, guarded by the parity check (ADR 0010). Env edits after construction other than temperature / beam / LET are not forwarded; the guard catches those that change labels, k or G.
- Results agree to solver tolerance (rtol 1e-6), not bit for bit (scipy's ODEPACK LSODA vs rsoda).
- Rust backend is LSODA only; solver options radkin does not support raise instead of being ignored.
- Not in this session: batch helper expanding override dicts from a base cell, macOS wheels, PyPI publishing, Rust-side events / dense output.
- Building from source needs network access to GitHub (git deps radkin → rsoda) and a Rust toolchain.

## Open Questions

- Exact radkin additions needed are discovered while implementing the binding; they stay small and go into one radkin PR (more if needed), merged by the user before the rev bump.
- Whether radkin's `Config` can be built from a JSON/TOML table of the merged `ConfigDict` as is, or needs a small adapter (reaction strings, beam section shapes): checked in the first ticket.

## Notes

- radkin's AGENTS.md lists radiopyo at `${DEV_DIR}/PYTHON/radiopyo`; the checkout is `C:\DEV\python\radiopyo` (same path on a case-insensitive FS).
- maturin is not installed yet; install it in the `radiobio` env.
- Licenses compatible: radiopyo LGPL-3.0-only, radkin LGPL-3.0-or-later.
- Glossary terms added: Rust backend, Parity guard, Batch run. ADR 0010 records the config-rebuild boundary.
