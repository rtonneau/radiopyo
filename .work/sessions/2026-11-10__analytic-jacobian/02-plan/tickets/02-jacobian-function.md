# Ticket 02: jacobian-function

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] Public pure `jacobian(t, y, beam, mat, force_positive=False)` in `ode_functions_solver.py`, in `__all__`, returning a dense (n_dyn, n_dyn) array for y in µmol/L.
- [ ] Mass-action derivatives via `reactant_index` (product rule, repeated reactants), Michaelis-Menten derivatives in enzyme and substrate; radiolytic and fixed species contribute nothing.
- [ ] Piecewise (ADR 0003): blocked reactions (k and enzymatic) contribute 0; rows zeroed by `force_positive` are 0.
- [ ] Test vs central finite-difference Jacobian of `derive` on bundled configs: beam on/off, enzymes, fixed species, an active depletion block, `force_positive`.

**Files to Touch:**
- `radiopyo/simulation/ode_functions_solver.py`
- `tests/simulation/test_jacobian.py`

**Verification Step:**

Run:
```bash
python -m pytest -q tests/simulation/test_jacobian.py && python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
Jacobian tests pass; full suite, ruff, mypy clean.

**Notes:**

Vectorise (no per-reaction Python loops if avoidable). Keep the switch states away from the threshold in the FD test.
