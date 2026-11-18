# Ticket 03: chemistry-check

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] `SimEnv.check_chemistry() -> list[ChemistryFinding]` (frozen dataclass: `kind`, `species`, `reactions`, `message`) and `UnitCell.check_chemistry()` delegating to it; exported.
- [ ] Species aliases: labels equal once their `_r` tokens are removed, plus table e_aq/e_minus and H_plus/H3O_plus; one finding per group.
- [ ] Charge balance per k-reaction, enzymatic and pair reaction from `_plus`/`_minus` tokens (e_aq = -1; negative products counted on the consumed side; enzyme, `trash`, radiolytic skipped).
- [ ] Dead reactants: dynamic species consumed but with no initial concentration, no radiolytic yield, no producing reaction (pH species H3O_plus/OH_minus have an initial value, so never dead).
- [ ] Unused species: given in `[concentrations.*]` or radiolytic yields but in no reaction.
- [ ] Run at the end of `SimEnv.from_dict`, logging one warning per finding only for a chemistry fingerprint not yet seen in the process; a manual call always logs and returns.
- [ ] All bundled configs (`PURE_WATER_CONFIG_FILE`, `LABARBE_CONFIG_FILE`, a config including `reactions_biology.toml`) give zero findings; full pytest, ruff, mypy pass.

**Files to Touch:**
- `radiopyo/simulation/chemistry_check.py` (new, banner layout)
- `radiopyo/simulation/environment_build.py`
- `radiopyo/simulation/unit_cell.py`
- `radiopyo/simulation/__init__.py` / `radiopyo/__init__.py` (exports)
- `tests/simulation/test_chemistry_check.py`

**Verification Step:**

Run:
```bash
$PY -m pytest tests/simulation/test_chemistry_check.py -q && $PY -m pytest -q && $PY -m ruff check radiopyo && $PY -m mypy radiopyo
```

Expected:
All tests pass; ruff "All checks passed!"; mypy "Success: no issues found".

**Notes:**

If a bundled config gives a finding, first check the rule (false positive → fix the rule); if the config itself is wrong and fixing it changes results, stop and ask. Fingerprint: reaction labels in index order + dynamic/fixed species labels + labels with initial concentrations + radiolytic labels (not k values or temperature). Tests reset the module-level seen-set (fixture).
