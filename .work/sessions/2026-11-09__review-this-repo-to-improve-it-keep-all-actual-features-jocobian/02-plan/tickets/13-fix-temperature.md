# Ticket 13: fix-temperature

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] Temperature has one source of truth (UnitCell vs `env.parameters["temperature"]` reconciled).
- [ ] The fixed-species vector, `initial_values()`, `ODEResult.to_pandas` and `__getitem__` use the simulation temperature.
- [ ] Before/after report at 25 °C (expected: unchanged) and at a non-default temperature, shown to the user before golden regeneration.

**Files to Touch:**
- `radiopyo/simulation/environment_build.py`
- `radiopyo/simulation/unit_cell.py`
- `radiopyo/simulation/ode_result.py`
- `tests/simulation/test_temperature.py`
- `CHANGELOG.md`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass; goldens at default temperature unchanged.

**Notes:**

Mind the mixed units: k polynomials use 1/K, pKa and initial concentration use °C. Document which unit each API takes.
