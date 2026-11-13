# Session: clone docs/issues/fenton/Fenton_Fixed_Showcase.ipynb and make it run simulation with rust backend

**Date:** 2026-11-13T14:24:28.347Z
**Status:** Grill phase complete

## Problem Statement

`docs/issues/fenton/Fenton_Fixed_Showcase.ipynb` shows the Fenton runaway case after the overshoot-relaxation fix with the Python backend only. We want the same showcase run entirely with the Rust backend (radkin, `backend="rust"`), as a separate notebook.

## Context & Constraints

- **Current behavior:** the showcase runs every simulation with the Python backend: three 1200 s runs (`UHDR`, `conv`, `UHDR_Freeze`), a frequency scan of 0.10 to 0.16 Hz with and without Fe2+, and a 12x12 ROOH map computed by `fenton_map.compute_map` in a `ProcessPoolExecutor` and cached in `.map_cache/`.
- **Pain point:** there is no Rust version of this case. The Python map takes minutes without its cache.
- **Dependencies:** the original notebook stays untouched (its uncommitted kernel-version diff is the user's). Both TOML files are used unchanged. The Python path of `fenton_map` keeps its behavior and its cache key, so the existing `.map_cache` file stays valid. The Rust backend needs a `UnitCell` built from a config: `from_config_dict(load_config(fe))` works (checked: the UHDR 1200 s run takes 0.05 s and gives ROOH 2.40 µmol/L, as in `ANALYSIS.md`).
- **Tech stack:** radiopyo, the `radiopyo-rs` extension (`rp.rust.run_batch`, `backend="rust"`), Jupyter/nbconvert in the `radiobio_env` conda env, numpy/pandas/matplotlib.

## Success Metrics

- `docs/issues/fenton/Fenton_Fixed_Showcase_Rust.ipynb` exists as a pure clone of the original: same sections, cells and figures, with every simulation run with `backend="rust"` and no Python-backend runs.
- `fenton_map.compute_map(backend="rust")` returns the same dict layout as the Python path, computed with `rp.rust.run_batch`, with the same feasibility (NaN) rules and the same "ROOH at max dose" extraction. Its values agree with the Python map within solver tolerance.
- The notebook runs top to bottom in `radiobio_env` without errors (nbconvert --execute) and is committed with its outputs. The `ruff`/`mypy` state of `radiopyo/` is unchanged, and `pytest` still passes.

## Architecture & Approach

1. **Notebook:** copy `Fenton_Fixed_Showcase.ipynb` to `Fenton_Fixed_Showcase_Rust.ipynb`. Add `backend="rust"` to every `.run()` call (sections 2 and 3). Section 4 calls `compute_map(n=12, period=PERIOD, backend="rust", cache=None)`, so the map is computed live and its time is printed. The markdown says that the Rust backend is used, and the claims about caching and process pools are updated. There is no Python-vs-Rust comparison (the user chose a pure clone).
2. **`fenton_map.py`:** `compute_map` gets `backend: "python" | "rust" = "python"`. The feasibility checks (the pulse period must be longer than the pulse, the max dose must be reached within `period`) and the extraction of ROOH at max dose / at the end are factored into helpers shared by `compute_point` and the Rust path. The Rust path builds one cell per feasible point (a deepcopy of the per-Fe cell with `set_beam(beam, LET=LET)`) and calls `rp.rust.run_batch(cells, [1e-9, period], mode="chunked", max_step_size_on=..., max_step_size_off=..., return_exceptions=True)`. A failed point is NaN. The cache key adds `"backend"` only when it is not `"python"`.
3. Execute the notebook with nbconvert in `radiobio_env` and commit it.

## Assumptions & Trade-offs

- Pure clone: no parity or timing cell against Python (the user's choice). The values are expected to match the Python notebook within solver tolerance, but the notebook does not show that.
- The frequency scan stays a plain loop of `.run(backend="rust")` calls, mirroring the original cell. Only the map uses `run_batch`.
- `fenton_map.py` is an issue-folder helper, not package code: no tests are added under `tests/`. The Rust map path is checked against the Python map (cached) by hand.
- The notebook is not wired into CI.

## Open Questions

None.

## Notes

- The size is Standard (light): one notebook and one helper change.
- The kernel is `radiobio_env`. CLAUDE.md's `E:\miniconda3` path is another machine's.
