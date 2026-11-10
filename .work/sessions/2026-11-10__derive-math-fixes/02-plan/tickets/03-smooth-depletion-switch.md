# Ticket 03: smooth-depletion-switch

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] The hard depletion block is replaced by the depletion switch: each k and enzymatic rate is multiplied, for each of its negative products, by c⁺/(c⁺ + ε), ε = DEPLETION_LIMIT = 1e-12 mol/L (c⁺ = clipped concentration).
- [ ] `jacobian` includes the switch derivative ε/(c⁺ + ε)² by the product rule (0 for c < 0); the Jacobian is no longer piecewise except at the clipping kink.
- [ ] New test: analytic Jacobian matches a finite-difference Jacobian on bundled configs at random positive states, at states with a negative product near ε, and with a clipped species.
- [ ] Existing depletion tests (test_reaction_rates.py, test_enzymatic_depletion.py, test_jacobian.py) are rewritten from the step to the switch (e.g. rate halved at c = ε, ~unchanged at 1e-6).
- [ ] Before/after golden numbers and run statistics (nfev, njev) on a bundled config recorded in the ticket log; only failing parts re-recorded.
- [ ] pytest, ruff and mypy at zero.

**Files to Touch:**
- `radiopyo/simulation/ode_functions_solver.py`
- `tests/simulation/test_reaction_rates.py`
- `tests/simulation/test_enzymatic_depletion.py`
- `tests/simulation/test_jacobian.py`
- `tests/golden/data/`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All tests pass, including the Jacobian vs finite-difference test; ruff and mypy report no errors.

**Notes:**

ADR 0005 (supersedes the piecewise part of ADR 0003); glossary term "Depletion switch". Keep the module docstring in sync (it describes the block). Function names are an implementation choice; keep `DEPLETION_LIMIT` exported. Result-changing fix: same golden procedure as ticket 02.
