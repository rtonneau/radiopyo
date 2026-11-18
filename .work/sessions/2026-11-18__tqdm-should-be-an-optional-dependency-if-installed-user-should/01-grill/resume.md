# Session: tqdm should be an optional dependency. If installed, user should be able to silent/disable it.

**Date:** 2026-11-18T16:30:24.278Z
**Status:** Grill phase complete

## Problem Statement

tqdm must be an optional dependency, and when it is installed the user must be able to silence the progress bar of `SimChunksIter.run`.

## Context & Constraints

- **Current behavior:** tqdm is already optional: it is the `tqdm` extra (`pyproject.toml`) and `SimChunksIter.run` (`radiopyo/simulation/unit_cell.py`) falls back to the plain iterator on `ImportError`; `tqdm.notebook` is used under Jupyter. It is the only tqdm use in `radiopyo`.
- **Pain point:** when tqdm is installed the bar is always shown; the only workaround is to consume the iterator by hand (`examples/tonneau_2025_pmb/paper_runs.py`).
- **Dependencies:** keep tqdm out of the core dependencies; default behavior unchanged (bar shown when tqdm is installed); the Rust backend path has no bar and stays as is.
- **Tech stack:** Python 3.11, tqdm (optional), pytest, ruff, mypy.

## Success Metrics

- `SimChunksIter.run(progress=False)` never imports or shows tqdm, even when it is installed; `progress=True` shows it when installed.
- `radiopyo.set_progress_bar(False)` disables the bar for every later `run()` called with `progress=None` (the default); an explicit `progress=` on the call overrides the global setting.
- Without tqdm installed, `run()` works with any setting (no error); tests cover the three cases; ruff, mypy and the full pytest suite stay green.

## Architecture & Approach

New module `radiopyo/utils/progress.py`: a module-level flag, `set_progress_bar(enabled: bool)`, `progress_bar_enabled() -> bool`, and a helper that wraps an iterable in `tqdm` (or `tqdm.notebook.tqdm` under a Jupyter kernel) when enabled and importable, else returns it unchanged. `set_progress_bar` and `progress_bar_enabled` are exported from `radiopyo`. `SimChunksIter.run(backend="python", progress: bool | None = None)`: `None` reads the global setting. README / CHANGELOG document the switch; `paper_runs.py` can use `run(progress=False)`.

## Assumptions & Trade-offs

Both a per-call keyword and a global switch (user's choice): the global flag is process-wide state, the keyword overrides it. No environment variable of our own (tqdm's own `TQDM_DISABLE` may still apply). No progress bar is added to the Rust backend or `rust.run_batch`.

## Open Questions

N/A

## Notes

Sized Small: one change in one method plus a small module, tests and docs.
