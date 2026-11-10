# Ticket 12: fix-run-burst

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] `run()` forwards arguments by keyword: the initial burst is applied on beam ON and `force_positive` reaches the solver.
- [ ] A user-supplied `y0` is never mutated; `SimChunksIter` no longer overwrites the previous chunk's last stored point.
- [ ] `run(use_beam=False)` adds no initial burst; `ODEResult.integrate_production_of` is consistent with whether a burst was applied.
- [ ] Missing beam raises `NoBeamDefinedError`.
- [ ] A before/after report (key species final concentrations and integrals for the golden cases) is written to the scratch dir and shown to the user; goldens regenerated only after approval.

**Files to Touch:**
- `radiopyo/simulation/unit_cell.py`
- `radiopyo/simulation/ode_result.py`
- `tests/golden/data/*`
- `tests/simulation/test_run_options.py`
- `CHANGELOG.md`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass after approved golden regeneration.

**Notes:**

Stop for user approval before regenerating goldens (result-changing fix).
