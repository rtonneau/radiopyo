# Ticket 01: paper-runs-module

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `examples/tonneau_2025_pmb/paper_runs.py` follows the module banner layout, is fully typed and has no matplotlib import.
- [ ] It defines the scenario of paper tables 3/5: pure-water config, LET 1.17, `variable_GValues` on, `SinglePulseBeam` (dose rate, dose), `prepare_single_pulse` over (0, 500 s), max ON step = pulse / 500, max final step 0.1 s, [O2] = pO2 % × 1300 µM.
- [ ] Grid constants: 19 dose rates (0.1–1e7 at two per decade, plus 50 and 5e6), 19 pO2 values (incl. 0, 0.1, 0.5, 1, 2, 4, 5, 6, 12, 16, 21 %), doses 10/20/30 Gy; bar subsets (table-5 dose rates × 9 bar pO2) and time-series subsets (table-5 dose rates × 1 %, 5 %).
- [ ] One-run analysis returns final O2/H2O2, ΔO2, gross production/consumption of O2 and H2O2 (radiolytic included), and O2/H2O2/pH on a fixed log time grid.
- [ ] Net yield = slope of a linear fit of final [H2O2] over the doses, in #/100 eV; full-grid driver and npz save/load with fixed keys.
- [ ] The grid corners (0.1 and 1e7 Gy/s × 0 and 21 % × 10 and 30 Gy) run without error.

**Files to Touch:**
- `examples/tonneau_2025_pmb/paper_runs.py`

**Verification Step:**

Run:
```bash
PYTHONPATH=. python -c "import importlib.util as u;s=u.spec_from_file_location('p','examples/tonneau_2025_pmb/paper_runs.py');m=u.module_from_spec(s);s.loader.exec_module(m);print(m.run_point(1e7,5.0,30.0).final_h2o2)"
```

Expected:
About 1.05 (µM), no exception; corner runs succeed.

**Notes:**

Use `integrate_production_of` / `integrate_consumption_of` of `ODEResult` (they include the radiolytic yield and the initial burst). pH = -log10([H3O_plus] mol/L). Logging through `logging.getLogger("radiopyo")`, no print in library-like code (the module is an example, but keep the style).
