# Ticket 02: overshoot-relaxation

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] `derive`: for each dynamic species with y_j < 0, dy_j/dt += J_jj(0⁺)·y_j, where J_jj(0⁺) is the diagonal Jacobian entry at the clipped state (y_j read at 0⁺); other rows unchanged (they read the species as 0).
- [ ] `jacobian`: for negative columns j, entry (j, j) = J_jj(0⁺) and entries (j, k) include y_j·∂J_jj(0⁺)/∂y_k; matches central finite differences of `derive` on states with negative species (away from kinks).
- [ ] No full extra Jacobian per derive call: only the negative columns are computed; states with no negative species have identical output to before.
- [ ] `ODEResult.reconstruct_reactions` unchanged (rates are still reaction rates); docstrings of the module explain the term and cite ADR 0007.
- [ ] Unit tests: a single stiff first-order decay crossing 0 is restored; a species with no consumption stays put; mass balance of other species unaffected.

**Files to Touch:**
- `radiopyo/simulation/ode_functions_solver.py`
- `tests/simulation/test_overshoot_relaxation.py`
- `tests/simulation/test_jacobian.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
all tests pass (golden tests may differ only within their tolerance; if a golden test fails, report the drift, do not re-record), zero ruff/mypy errors.

**Notes:**

Prototype: `.scratch/tests/2026-11-10__lsoda-stalls-fenton-high-concentration/proto.py` (`derive_p1`, `jac_p1`; the prototype Jacobian omits the y_j·∂J_jj/∂y_k term). J_jj(0⁺) for k-reactions: Σ_i S_ji·(∂r_i/∂y_j at 0⁺)·switch_i + Σ_i S_ji·r_i·∂switch_i/∂y_j; same for enzymatic. A species appearing twice in a reactant (2A) has derivative 0 at 0⁺.
