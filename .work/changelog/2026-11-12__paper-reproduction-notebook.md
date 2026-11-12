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
- `examples/tonneau_2025_pmb/gf_fit_comparison.ipynb`: the paper grid re-run with the GF(DR) fit printed in the
  paper (0.982 DR^-0.09, `paper_runs.run_grid(gf_fit="paper")`) and compared with the shipped fit. The paper's fit
  overshoots the quoted UHDR values by up to 35 %; the CDR gaps do not depend on the fit. Golden data:
  `tests/golden/data/paper_tonneau_2025_paper_gf.npz` (`record_golden.py --only paper-gf`).
- `examples/tonneau_2025_pmb/pre_refactor_comparison.ipynb`: the paper grid run with the pre-refactoring code
  (commit 2261b92, `record_pre_refactor.py` from a separate checkout, `data/paper_tonneau_2025_2261b92.npz`)
  against the current code. Final concentrations agree to 1e-8 relative, the minimum pH to 1e-7, the integrals
  to 1e-4; the gaps to the paper are the same before and after.
- Scenario `"original"` of the paper reproduction (`paper_runs.SCENARIOS`, `run_point` / `run_grid(scenario=)`,
  `SCENARIO` switch of `tonneau_2025_pmb.ipynb`): the inputs of the original runs whose values are quoted in the
  paper's text, the author's config `examples/tonneau_2025_pmb/data/CONFIG_PAPER_01.toml` (e_aq + O2, H_r + O2,
  H_r + HO2 rate constants differ from the bundled config) at LET 0.56 keV/µm. It reproduces those runs to about
  0.03 % and the quoted values to about 1 %; the paper's methods (LET 1.17, GF fit of fig. 3B) describe a later
  rerun. The default scenario `"methods"` (bundled config, LET 1.17) is unchanged. Golden data:
  `tests/golden/data/paper_tonneau_2025_original.npz` (`record_golden.py --only paper-original`).
- Optional dependency extra `examples` (matplotlib).

### Fixed

- `ODEResult.detail_consumption_of` (used by the consumption integrals) raised an `IndexError` under pandas 3
  for a species with an equilibrium of net production (H2O2, OH_minus). Values computed before the fix do not
  change.
