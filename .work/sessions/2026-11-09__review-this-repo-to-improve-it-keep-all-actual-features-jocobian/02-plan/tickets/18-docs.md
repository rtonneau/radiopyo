# Ticket 18: docs

**Model:** fast
**Effort:** medium

**Acceptance Criteria:**
- [ ] README: version badge, `from_toml` (not `from_toml_file`), no reference to missing `basic_water_radiolysis.toml`, `max_dose=` keyword in beam examples, `[[beam.macro_pulsed]]` documented, acid_base keys match the code, G-value list incl. H3O_plus = e_aq note (ADR 0002) and `variable_GValues`, dependency list incl. h5py/tqdm, GitHub install URL, n·k convention.
- [ ] Every code snippet in README runs against the final code (scripted check in the scratch dir).
- [ ] CHANGELOG complete for the session; CLAUDE.md architecture section updated (matrix objects, module splits, removed features, current test state).
- [ ] Seeds recorded for: Jacobian rewrite, `_scale_factor` dose-rate range, H2O initial-cc polynomial sign check, and any finding deferred by earlier tickets.

**Files to Touch:**
- `README.md`
- `CHANGELOG.md`
- `CLAUDE.md`
- `examples/notebook_basic_usage.ipynb` (only if a call changed)

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass; README snippet check script exits 0.

**Notes:**

Edit the notebook only through the Jupyter MCP tools if it is open; otherwise leave it untouched unless a call is broken.
