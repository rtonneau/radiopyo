# Ticket 01: enzymatic-depletion-block

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `EzMatrices` has a read-only `negative_products` (n_species, n_ez) bool array built in `make_enzymatic_matrix` from each product's `ope`.
- [ ] `reaction_rates` zeroes an enzymatic rate while one of its negative products is at or below `DEPLETION_LIMIT`, through a helper shared with the k-reactions.
- [ ] A test with a toy config shows the block on and off; goldens unchanged; pytest, ruff, mypy clean.

**Files to Touch:**
- `radiopyo/utils/sim_types.py`
- `radiopyo/simulation/environment_build.py`
- `radiopyo/simulation/ode_functions_solver.py`
- `tests/simulation/test_enzymatic_depletion.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All tests pass (golden included), ruff and mypy report no errors.

**Notes:**

Bundled catalase has no negative product: no result change expected. Keep EzMatrices zero-reaction case valid.
