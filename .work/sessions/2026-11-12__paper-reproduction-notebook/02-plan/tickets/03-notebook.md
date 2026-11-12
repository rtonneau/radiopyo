# Ticket 03: notebook

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `examples/tonneau_2025_pmb/tonneau_2025_pmb.ipynb` explains the scenario (paper tables 3/5, GF fit gap per ADR 0009) and has a `RECOMPUTE = False` flag: False loads the golden npz, True re-runs the grid via `paper_runs.py` and reports the max relative difference to the golden file.
- [ ] It draws fig 4 (O2), fig 7 (H2O2), fig 10 (pH) as 2×2 panels (CDR/UHDR × 1 %/5 %), fig 5 (ΔO2 map, %pO2 and µM), fig 8 (net H2O2 yield map, #/100 eV), fig 6 and fig 9 (gross production/consumption bars with net marker).
- [ ] It shows an informational table: values quoted in the paper text vs radiopyo, with relative difference.
- [ ] It executes top to bottom with nbconvert; committed with outputs stripped.

**Files to Touch:**
- `examples/tonneau_2025_pmb/tonneau_2025_pmb.ipynb`

**Verification Step:**

Run:
```bash
python -m jupyter nbconvert --to notebook --execute examples/tonneau_2025_pmb/tonneau_2025_pmb.ipynb --output-dir .scratch/tests/2026-11-12__paper-reproduction-notebook
```

Expected:
Exit 0; the executed copy holds the 7 figures and the comparison table.

**Notes:**

Inspect the rendered figures against the paper pages in `.scratch/paper/` before stripping outputs. Use `nbconvert --clear-output --inplace` to strip.
