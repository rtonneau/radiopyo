# Ticket 02: jacobian-function

**Status:** ✅ Done

## Local Test Result

Command: `python -m pytest -q tests/simulation/test_jacobian.py && python -m pytest -q && ruff check radiopyo && mypy radiopyo`

- `tests/simulation/test_jacobian.py`: 22 passed
- full suite: 580 passed (was 558 before the 22 new tests)
- `ruff check radiopyo`: All checks passed; `mypy radiopyo`: no issues found in 42 source files (the new test file is also ruff-clean)

## Review Notes

- AC1: `jacobian(t, y, beam, mat, force_positive=False)` added to `radiopyo/simulation/ode_functions_solver.py` and to `__all__`; pure, returns a dense (n_dyn, n_dyn) float64 array (checked by shape/purity tests). Units cancel (y and dy/dt both scaled by 1e6), noted in the docstring.
- AC2: `mass_action_derivatives` (product rule over the reactant slots, one loop over the few slots, `np.add.at` so a repeated reactant such as 2 A sums both slots; unused slots land in a dropped extra column) and `michaelis_menten_derivatives` (enzyme and substrate; `np.add.at`). Radiolytic and fixed species contribute nothing (columns of fixed species are sliced off). No per-reaction Python loop. Hand-computed checks for 2 A + F -> E (2kAF) and the MM partial derivatives.
- AC3: blocked k and enzymatic reactions have their derivative rows zeroed with `negative_product_block`; with `force_positive` the rows where `y < DEPLETION_LIMIT and dydt < 0` are zeroed (same condition as `zero_depleted_consumption`, dydt from `reaction_rates`). Tests: active k block, active enzymatic block, force_positive rows zero (toy env and all-depleted real configs).
- AC4: central finite-difference Jacobian of `derive` on the labarbe (enzyme + fixed species) and water configs, beam on/off x force_positive, plus toy envs. The FD test tolerance is rtol*column max plus an estimate of the FD rounding noise (eps*|f|/h), step 1e-3 relative; a tighter relative step failed only because of rounding noise on tiny entries, not on the analytic values.
- Bundled configs have no negative products, so the active-block cases use toy environments.
- Changed after review: none beyond fixing my own first test drafts (tolerance, and the force_positive scenario: in a random state every species has a source, so I use a fully depleted state where the constant species still consume).
- Not done (out of scope): wiring into `UnitCell`/`use_jac`, CLAUDE.md note about the missing `jacobian` (still says there is none).

Parent: base 7b9f95891388fb3221ab186d4181f8809307f31f, no pre-existing changes. Independent review (opus, high; fresh context, base 7b9f958): no findings; its gaps (mypy, full pytest) closed by the parent run: 580 passed, ruff and mypy clean. Accepted (implementer: opus, high).

## Blockers / Challenges

- FD noise on small entries with a 1e-6 relative step (resolved by a noise-aware tolerance and a 1e-3 step; the central difference is exact for the quadratic/cubic mass-action terms).
- A bash heredoc creating the test file failed once (shell quoting); recreated with the Write tool.
- No blockers.

## Commits

- 694616c feat: analytic Jacobian of derive (ticket 02)

## Time Spent

12m (ticket-start.js to ticket-complete.js)
