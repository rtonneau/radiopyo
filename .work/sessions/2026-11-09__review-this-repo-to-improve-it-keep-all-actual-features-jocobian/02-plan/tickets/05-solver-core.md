# Ticket 05: solver-core

**Model:** heavy
**Effort:** xhigh

**Acceptance Criteria:**
- [ ] `Ki_mat`/`Gi_mat`/`Ez_mat` positional tuples replaced by named frozen dataclasses (e.g. `KMatrices`, `GVector`, `EzMatrices`) inside `SimMatrices`; no `mat[0]`-style indexing remains.
- [ ] One pure function computes per-reaction rates from (t, full molar state, matrices, beam) and is used by both `derive()` and `ODEResult.reconstruct_reactions()`; the latter is vectorised over time points (no per-row `iloc` loop) and the hardcoded `1e-12` is replaced by the same rule `derive` uses.
- [ ] µM↔M conversion and fixed-species append/strip happen in one helper; works with zero fixed species.
- [ ] The stale `jacobian()` body and the 4-D `JRki` tensor are removed; `run(use_jac=True)` raises `NotImplementedError` with a message saying the Jacobian is being rewritten.
- [ ] Golden tests pass (derive at rtol 1e-12; reconstruct_reactions and runs within the documented tolerances).

**Files to Touch:**
- `radiopyo/simulation/ode_functions_solver.py`
- `radiopyo/simulation/environment_build.py` (matrix builders, `SimMatrices`)
- `radiopyo/simulation/ode_result.py`
- `radiopyo/simulation/unit_cell.py` (call sites)
- `tests/simulation/test_reaction_rates.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass, zero lint/type errors.

**Notes:**

This is the Jacobian prerequisite: keep the rate function differentiable-friendly (stoichiometry maps reusable, the "negative product" zeroing as an explicit separate step). Preserve column order of `reconstruct_reactions` (k → G → Ez) as today even though global indices are k → Ez → G; note the discrepancy in a docstring. Check `Ez_iy` `=` vs `+=` (plausible bug): if real and it changes numbers, leave it and record it for a Phase C ticket.
