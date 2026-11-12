# Implementation Plan

**Session:** paper-reproduction-notebook
**Date:** 2026-11-12T14:49:12.677Z
**Estimated effort:** 1 day

## Strategy

Build the shared run/analysis module first (everything else consumes it), then record the golden data and wire the fast golden test on top of it, then write the notebook that plots the golden data (or recomputes), and close with packaging/docs. The full-grid recording (~10 min) runs once, in ticket 02.

## Tickets Overview

- **Ticket 1:** `examples/tonneau_2025_pmb/paper_runs.py`: paper scenario, grid, one-run analysis, net yield fit, full-grid driver, npz save/load.
- **Ticket 2:** Golden data: `record_golden.py --only paper` writes `tests/golden/data/paper_tonneau_2025.npz`; new fast test re-runs a subset and matches it.
- **Ticket 3:** Notebook `examples/tonneau_2025_pmb/tonneau_2025_pmb.ipynb`: figs 4/5/6/7/8/9/10, comparison table, `RECOMPUTE` flag; executed once with nbconvert, committed with outputs stripped.
- **Ticket 4:** Optional `examples` extra (matplotlib), README examples section, changelog fragment.

## Sequencing Rationale

The notebook and the golden test both import `paper_runs.py`, so it lands first; the notebook loads the golden npz by default, so the golden data must exist before it; docs describe the finished files.

## Risks & Mitigation

- **Risk:** some grid point fails (solver failure, `ShortBeamOnError`) → **Mitigation:** doses 10/20/30 Gy keep pulses ≥ 1 µs at 1e7 Gy/s; ticket 1 runs the grid corners before ticket 2 records the full grid.
- **Risk:** golden test not reproducible across machines (BLAS, solver steps) → **Mitigation:** compare end-of-run scalars and resampled series at rtol ≈ 1e-6 (as the existing golden runs do), loosen only with evidence.
- **Risk:** the npz grows large → **Mitigation:** time series resampled on a fixed log grid (a few hundred points), maps and bars are small arrays.
- **Risk:** test importing a module from `examples/` → **Mitigation:** load it by path with importlib, as `record_golden.py` already does for `test_golden.py`.

## Assumptions

- The radiobio_env on this machine (Python 3.11, matplotlib, nbconvert) is used, with `PYTHONPATH=.` since radiopyo is not installed in it.
- The code's GF fit is used unchanged (ADR 0009).
