# Ticket 03: vectorise-helpers

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] `depletion_switch`, `depletion_switch_derivatives` and `mass_action_derivatives` have no `np.delete` and no per-reaction / per-column Python loop (a loop over the fixed small number of reactant slots is fine)
- [ ] Same results as the reference within rtol 1e-12, including at c = 0 and for repeated reactants
- [ ] Full pytest, ruff, mypy pass; final before/after bench numbers in the commit log and a CHANGELOG fragment

**Files to Touch:**
- `radiopyo/simulation/ode_functions_solver.py`
- `CHANGELOG.md` or `.work/changelog/`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo && python examples/bench_overshoot_pulsed.py
```

Expected:
All tests pass, zero errors, wall time below ticket 2's.

**Notes:**

Product over "others" without division: prefix/suffix cumulative products (exact zeros safe). For the depletion switch, a masked product (`np.where(negative, per_species[:, None], 1).prod(axis=0)`).
