# Ticket 05: poly-k-warning

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] Building a reaction with a single reactant of stoichiometry n > 1 and a non-scalar k (polynomial, Arrhenius, equilibrium) logs one warning on the `radiopyo` logger saying the k is used as given (ADR 0001).
- [ ] A scalar k on the same reaction logs no such warning; bundled configs build without new warnings unless they hit the case (report it if they do).
- [ ] caplog tests for both cases.
- [ ] pytest, ruff and mypy at zero.

**Files to Touch:**
- `radiopyo/chemistry/k_reactions.py`
- `tests/chemistry/test_nk_warning.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All tests pass; ruff and mypy report no errors.

**Notes:**

Put the check where the n·k convention is applied (`_apply_nk_convention` or the Arrhenius / equilibrium builders); use `logging.getLogger("radiopyo")`, no print.
