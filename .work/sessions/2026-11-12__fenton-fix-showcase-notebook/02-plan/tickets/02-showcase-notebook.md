# Ticket 02: showcase-notebook

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] Notebook `docs/issues/fenton/Fenton_Fixed_Showcase.ipynb` with sections: recap of old symptoms; UHDR/CONV/UHDR_Freeze 1200 s runs (ROOH, H2O2, Fe3+, solver stats); frequency scan 0.10–0.16 Hz Fe=0 vs 2.2e-8; map of normalized ROOH at max dose styled like the reference image; map of ROOH at end time.
- [ ] Executed top to bottom with outputs saved, no solver failure.
- [ ] Uses the folder's TOML files unchanged; SOD label bug mentioned in markdown.

**Files to Touch:**
- `docs/issues/fenton/Fenton_Fixed_Showcase.ipynb`
- `docs/issues/fenton/fenton_map_cache.npz`

**Verification Step:**

Run:
```bash
cd docs/issues/fenton && python -m jupyter nbconvert --to notebook --execute Fenton_Fixed_Showcase.ipynb --output .scratch_check.ipynb
```

Expected:
Exits 0 (map loaded from cache).

**Notes:**

Build/edit via Jupyter MCP if available, else nbformat. Fe=0 variant: set the fixed Fe_plus_plus to 0 through the env dict / a temp TOML copy in scratch.
