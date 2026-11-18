# Ticket 01: build-warnings

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] An `[[reactions.acid_base]]` with `acid`/`base` keys and no `reaction` logs a warning showing the converted entry (`reaction = "<acid> <-> H3O_plus + <base>"`, `pKa`, and that `k_forward`/`k_backward` are required), then raises the same `AttributeError` as today with a message pointing to the warning.
- [ ] `H2O` in `[concentrations.fixed]` logs a warning that H2O cannot be replaced and builds normally with the built-in H2O (no `KeyError`).
- [ ] `H2O` in `[concentrations.initial]` logs the same warning; any other label of `[concentrations.initial]` matching no species logs a warning naming it.
- [ ] Tests cover the four cases; full pytest, ruff, mypy pass.

**Files to Touch:**
- `radiopyo/simulation/environment_build.py`
- `tests/simulation/test_config_warnings.py`

**Verification Step:**

Run:
```bash
$PY -m pytest tests/simulation/test_config_warnings.py -q && $PY -m pytest -q && $PY -m ruff check radiopyo && $PY -m mypy radiopyo
```

Expected:
All tests pass; ruff "All checks passed!"; mypy "Success: no issues found".

**Notes:**

`constant_species_from_dict` builds the fixed species before H2O is added: filter H2O out of the fixed dict there (or before) with the warning. The initial-cc loop currently uses `suppress(KeyError)`: replace with an explicit lookup + warning. A fixed species also listed in initial is not an "unknown label" (only BasicSpecies, not SimSpecies): keep it silent unless it is H2O.
