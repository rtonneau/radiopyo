# Ticket 02: clipped-concentrations

**Effort:** high

**Acceptance Criteria:**
- [ ] `mass_action_rates`, `michaelis_menten_rates`, `mass_action_derivatives`, `michaelis_menten_derivatives` read clipped concentrations max(c, 0) (enzyme and substrate included); `reaction_rates` keeps serving both `derive` and `ODEResult.reconstruct_reactions`.
- [ ] Jacobian: a species with c < 0 contributes 0 to every derivative it enters; at exactly c = 0 the derivative of the c > 0 side is used.
- [ ] New test: `2 A -> B` evaluated at A < 0 gives dA/dt = 0 (no further consumption); a short run started with a slightly negative A stays bounded.
- [ ] New test: mass balance, derive equals S·r for the full rate vector, and a conserved linear combination of species on a small config has zero derivative, including at states with negative and near-depleted species.
- [ ] Golden comparison run before re-recording; before/after numbers (max abs/rel deviation per failing part) recorded in the ticket log; only failing parts re-recorded with `--only`.
- [ ] pytest, ruff and mypy at zero.

**Files to Touch:**
- `radiopyo/simulation/ode_functions_solver.py`
- `tests/simulation/test_reaction_rates.py`
- `tests/simulation/test_jacobian.py`
- `tests/simulation/test_clipped_concentrations.py`
- `tests/golden/data/`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All tests pass, including the new overshoot and mass-balance tests; ruff and mypy report no errors.

**Notes:**

ADR 0004; glossary term "Clipped concentration". Result-changing fix: record before/after numbers before any re-record (`PYTHONHASHSEED=0 python tests/golden/record_golden.py --force --only <part>`). Never rewrite `tests/golden/data/golden_results.h5`. If no golden part fails, say so and re-record nothing.
