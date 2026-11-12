---
bump: minor
floor: minor
---
### Added

- `examples/tonneau_2025_pmb/tonneau_2025_pmb.ipynb`: reproduction of the main findings of Tonneau et al.
  2025 (Phys. Med. Biol. 70 235021). It shows the O2, H2O2 and pH time evolution (figs 4, 7, 10), the ΔO2 and net
  H2O2 yield maps over dose rate and pO2 (figs 5, 8), and the gross production / consumption bars (figs 6, 9).
  A table compares the results with the values quoted in the paper. The results use the GF(DR) fit that ships
  with radiopyo, which differs from the paper's by up to 28 % at 1e7 Gy/s, so most values are 3-8 % below the
  paper's. By default the notebook loads the recorded grid; `RECOMPUTE = True` re-runs it (about 10 min).
- `examples/tonneau_2025_pmb/paper_runs.py`: the runs behind the notebook. `run_point` handles one dose rate,
  pO2 and dose; `run_grid` runs the whole grid in parallel; `net_yield` and `save_grid` / `load_grid` are also
  provided.
- Golden data of the paper grid: `tests/golden/data/paper_tonneau_2025.npz`, recorded with
  `record_golden.py --only paper` (new `--workers` option), and a fast test that re-runs 8 points against it.
- Optional dependency extra `examples` (matplotlib).

### Fixed

- `ODEResult.detail_consumption_of` (used by the consumption integrals) raised an `IndexError` under pandas 3
  for a species with an equilibrium of net production (H2O2, OH_minus). Values computed before the fix do not
  change.
