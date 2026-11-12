# Ticket 02: derive-overhead

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] Per-call allocations removed or hoisted in `derive`, `jacobian`, `StateLayout.to_full_molar` and the mass-action / depletion helpers (`np.ones`, `np.broadcast_to` + concatenate), e.g. a preallocated/precomputed path for the 1-D state.
- [ ] Golden tests and `tests/simulation/test_overshoot_equivalence.py` pass unchanged; no golden data re-recorded.
- [ ] `examples/bench_overshoot_pulsed.py` wall time reported before/after in the commit log (expected lower).
- [ ] Full pytest, ruff and mypy pass.

**Files to Touch:**
- `radiopyo/simulation/ode_functions_solver.py`
- `radiopyo/utils/sim_types.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo && python examples/bench_overshoot_pulsed.py --repeat 1
```

Expected:
All tests pass; benchmark prints a wall time lower than before the change.

**Notes:**

Matrices are frozen/read-only (`SimMatrices`); any precomputed buffer must be an `init=False` field derived once, never a mutable per-call buffer shared across calls in a way that leaks into returned arrays.
