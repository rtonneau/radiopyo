# Ticket 01: paper-runs-module

**Status:** ✅ Done

## Local Test Result

Env: `radiobio_env` (Python 3.11, pandas 3.1.1), `PYTHONPATH=.`; `pandas-stubs` and `ruff` were missing from it and were installed (`pip install "pandas-stubs>=1.5.3" ruff`, the dev extra).

- Verification step: `run_point(1e7, 5.0, 30.0, integrals=True)` takes about 1.0 s and gives
  H2O2 final 1.0459 µM (expected about 1.05), ΔO2 0.4417 µM, pH min 5.859; gross H2O2 production 1.051 µM, consumption 0.0054 µM; gross O2 production 0.852 µM, consumption 1.294 µM.
- Corner runs at 5 % pO2 and 30 Gy match the session probes: 0.1 Gy/s gives 3.0635 µM H2O2, ΔO2 1.064, pH min 6.82; 1e3 Gy/s gives 1.7171 µM, ΔO2 0.641, pH min 6.43. 1e7 Gy/s, 0 % pO2, 10 Gy (1 µs pulse) runs: 0.185 µM, pH min 6.28.
- Net yield at 1e3 Gy/s, 5 %: H2O2 0.666 / 1.203 / 1.717 µM at 10/20/30 Gy, slope 0.507 #/100 eV.
- `.scratch/tests/2026-11-12__paper-reproduction-notebook/01-grid_smoke.py`: a 2×2×2 `run_grid`, serial (3.7 s) and `workers=2`, saved and reloaded through npz: every GRID_KEYS array equals the serial one (rtol 1e-10); time-series shape (2, 2, 300).
- `ruff check radiopyo examples/tonneau_2025_pmb/paper_runs.py`: clean. `mypy radiopyo`: no issues (42 files). `mypy examples/tonneau_2025_pmb/paper_runs.py`: no issues.
- `python -m pytest -q`: 793 passed.

## Review Notes

Base: `f23895a69a2b8c920b1a618e6983563ad02ca09e`. Pre-existing changes at start: only the gps-owned `.session-config.json`.

Inline review against the acceptance criteria:
- `examples/tonneau_2025_pmb/paper_runs.py` follows the banner layout, is fully typed, and imports no matplotlib.
- Grid: 19 dose rates (10^(k/2), k = -2..14, plus 50 and 5e6) × 19 pO2 values (bar values, 1 and 5 %, and fill values) × 10/20/30 Gy. Subsets: `PAPER_DOSE_RATES` (table 5), `BAR_PO2`, `TS_PO2`.
- `run_point` returns final O2/H2O2, ΔO2, pH min and O2/H2O2/pH on `TIME_GRID` (300 log points, 1e-7 to 500 s). The gross integrals are optional (`integrals=True`) and reported as positive amounts.
- `net_yield` fits `np.polyfit` over the doses and converts with `kr_to_ge`. `run_grid` (optional `workers`, `progress` callback) computes the integrals at the last dose only. `save_grid` / `load_grid` use the fixed `GRID_KEYS`.
- Logging: the INFO logs of radiopyo are silenced per run by a context manager that restores the logger level, so the global logger state (and pytest caplog) is untouched. One template `UnitCell` is cached and deep-copied per run.

Library fix found while probing the integrals:
- `ODEResult.detail_consumption_of` (`radiopyo/simulation/ode_result.py`) dropped equilibrium columns while iterating over them, which raised an IndexError under pandas 3 as soon as a species has an equilibrium with net production (H2O2, OH_minus).
- It now collects the labels first, then drops them. The outputs do not change where the old code ran.
- New regression tests in `tests/simulation/test_species_integrals.py` (7 tests):
  - consumption is finite and ≤ 0;
  - an equilibrium is never counted in both production and consumption;
  - production + consumption = net change + initial burst.

## Blockers / Challenges

- The dev env on this machine (`radiobio_env`, not `radiobio` as in CLAUDE.md) lacked `pandas-stubs` and `ruff`. Without the stubs, mypy reported 5 pre-existing import-untyped errors. Installing the dev extra fixed it.
- The pandas 3 bug in `detail_consumption_of` (above) had to be fixed before the gross consumption integrals could be computed.

## Commits

- 8c66628 feat: paper reproduction runs module, fix equilibrium consumption under pandas 3 (ticket 01)

## Time Spent

7m (ticket-start.js to ticket-complete.js)
