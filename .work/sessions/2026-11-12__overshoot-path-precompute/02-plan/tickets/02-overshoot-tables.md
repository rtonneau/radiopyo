# Ticket 02: overshoot-tables

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] A frozen read-only `OvershootTables` in `radiopyo/utils/sim_types.py`, an `init=False` field of `SimMatrices` built in `__post_init__` (via `object.__setattr__`), with per-dynamic-species k and enzymatic column subsets, sub-arrays, slot masks / held-once flags and negative-products-without-j
- [ ] `_k_relaxation_terms`, `_ez_relaxation_terms` and `_relaxation_sum` no longer call `dataclasses.replace`; `SimMatrices` constructor unchanged
- [ ] Equivalence test, golden, Jacobian, overshoot and Fenton stall tests pass; ruff and mypy zero errors
- [ ] Bench numbers after the change in the commit log

**Files to Touch:**
- `radiopyo/utils/sim_types.py`
- `radiopyo/simulation/ode_functions_solver.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo && python examples/bench_overshoot_pulsed.py
```

Expected:
All tests pass, zero lint/type errors, the overshoot share well below the baseline.

**Notes:**

Keep the banner layout and docstrings. The helpers (`mass_action_rates`, `depletion_switch`, ...) may take plain arrays instead of sub-dataclasses if that avoids building objects; keep the public signatures.
