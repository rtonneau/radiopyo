# Ticket 03: smooth-depletion-switch

**Status:** ✅ Done

## Local Test Result

```
python -m pytest -q && ruff check radiopyo && mypy radiopyo
615 passed in 84.61s
All checks passed!
Success: no issues found in 42 source files
```
Saved: .scratch/tests/2026-11-10__derive-math-fixes/evidence/03-verification.txt.
Mutation check: with `depletion_switch_derivatives` forced to 0, 15 of the 34 tests in test_jacobian.py fail (the near-eps finite-difference tests, product-rule tests), so the new tests do exercise the switch derivative.

## Review Notes

Implementation (radiopyo/simulation/ode_functions_solver.py):
- `negative_product_block` (bool step) replaced by `depletion_switch(y_full, mat)` (product over the negative products of c+/(c+ eps), eps = DEPLETION_LIMIT, c+ = max(c, 0); shape (n_r,) or (n_r, n_t)) and `depletion_switch_derivatives(y_full, mat)` ((n_r, n_species), product rule over the negative products, 0 for c < 0, 1/eps at c = 0). Both exported; `DEPLETION_LIMIT` still exported; `negative_product_block` removed from `__all__` (no other user in radiopyo/).
- `reaction_rates`: k and enzymatic rates are now law * switch. `jacobian`: product rule via `_product_rule(law, d_law, switch, d_switch)`; no longer piecewise except at c = 0. Module docstring updated (steps list, switch description, ADR 0005); `jacobian` docstring rewritten (awkward line break removed).
- README.md line 319 (enzymatic depletion sentence) rewritten for the switch.

Acceptance Criteria:
- switch c+/(c+ + eps) per negative product on k and enzymatic rates: done (unit-tested: 0 at c <= 0, 1/2 at eps, 2/3 at 2 eps, > 0.99 at 1e-10, > 1 - 1e-5 at 1e-6; scalar and vectorised).
- Jacobian includes eps/(c+ + eps)^2 by the product rule: done; closed-form product-rule values asserted (k and enzymatic toys).
- New FD test: tests/simulation/test_jacobian.py: `negative_cell` = Labarbe config + a k-reaction with two negative products (H2O2, OH_minus) + an enzymatic reaction with one (O2_r_minus), built from a temp TOML (the bundled configs have NO negative product, so they cannot exercise the switch). Tests: random positive states (3 seeds), negative products at 0.1/1/3/100 x eps (rel_step 1e-4), and states with clipped (negative) species incl. a negative product (columns of clipped species are exactly 0). Toy envs also at 0.1/1/3/100 x eps. Bundled labarbe/water FD tests unchanged and green.
- Existing depletion tests rewritten: test_reaction_rates.py (`test_negative_product_depletion_switch`), test_enzymatic_depletion.py (switch parametrised, "halved at eps / unchanged at 1e-6", vectorised), test_jacobian.py (exact formulas now carry the switch factor, which is 1 - 1e-6 at the random D ~ 1 umol/L used by the toys).
- pytest, ruff, mypy: zero (see above). `ruff format --check` reports ode_functions_solver.py, but only for a blank line after the module docstring that was already there at HEAD (not touched).

Golden procedure (result-changing fix, but NOT result-changing on the recorded scenarios): none of the golden scenarios (labarbe, pure_water, pulsed chunked, hdf) has a negative product (0 in both bundled configs), so every switch factor is exactly 1.0. Deviation of the current code vs the recorded data, run unchanged with .scratch/tests/2026-11-10__derive-math-fixes/02-golden_dev.py, is BIT-IDENTICAL before and after the change (diff of evidence/03-golden-deviation-before.txt and -after.txt is empty):
  derive/labarbe: max_abs 2.384e-07, max_rel 7.774e-16 (31/207 differ); derive/pure_water: 0 / 0
  beams: 0 / 0
  runs/labarbe_run species: max_abs 9.604e-11 (max_rel 1.294e+19 on values recorded as ~0); rate: 3.714e-13 / 1.000
  runs/labarbe_curve species: 1.065e-10 / 1.294e+19; rate: 8.752e-10 / 1.450e+01
  runs/chunked_run species: 1.314e-10 / 2.563e-07; rate: 3.711e-12 / 2.563e-07
  hdf/hdf_labarbe rebuilt_rate: 0 / 0 (hdf_expected.npz was already patched in ticket 02)
These deviations are the ones already present at the parent base (tickets 01/02), all within the golden tolerances (golden tests pass). No golden part fails, so nothing was re-recorded: tests/golden/data/ is untouched, golden_results.h5 sha256 unchanged (cb2d1bfbdaf73d0ab10985e2161226518a20ac52dc7d57fa01736643b580cc9e). Consequently the derive.json stale `*_beam_on_force_positive` entries are still there and `test_derive_matches_golden` keeps its ticket 01 subset check on dydt cases (the exact-set restoration was conditional on a derive re-record). The parent may decide to re-record derive purely to drop them.

Solver statistics (Labarbe config, LSODA, use_jac=True; evidence/03-solver-stats-before.txt / -after.txt, 3 repeats each):
  pulsed chunked (3 pulses of 2 Gy, 2 us on, period 10 ms, t_span 1e-7..0.05, LET 0.201): before nfev=8574 njev=274 nlu=274 (4787 points), wall 0.75-1.21 s; after nfev=8574 njev=274 nlu=274, wall 0.74-0.79 s
  plain run 1e-7..10 s: before nfev=4114 njev=144 nlu=144 (2381 points), wall 0.34-0.40 s; after identical counts, wall 0.34-0.45 s
  Identical counts (no negative product in the bundled configs, switch is 1); wall-time differences are noise. `depletion_switch` returns early (np.ones) when a reaction family has no negative product.

Not touched (out of ticket scope, for the parent): CLAUDE.md (item 4 still describes the zeroing, `force_positive` in the signatures), CHANGELOG.md / .work/changelog entries, tests/simulation/test_clipped_concentrations.py (unchanged, passes; its DEPLETED states now slow rather than block).

Parent: base 889480465f458343adb5550ec2cddc570b3358bf; pre-existing changes: only gps session state (.session-config.json, 03-implement/), owned by the parent.
- Parent review round 1: fresh reviewer (opus, base 8894804) REVIEWED, 0 correctness defects; math (a)-(f) verified incl. in-memory FD on reactant-and-negative-product, stoichiometry-2 negative product, fixed species, vectorized inputs. Findings: P3-1 stale 'block' wording in sim_types.py (KMatrices/EzMatrices docstrings) and environment_build.py comment; P3-2 no regression test for 'A + D -> C - 2 D'. Disposition: both accepted; one correction pass, inline by the parent (small, mechanical).
- Correction pass (parent): docstrings in sim_types.py (KMatrices rate line, KMatrices/EzMatrices negative_products) and the environment_build.py comment reworded for the switch (ADR 0005); new test `test_switch_reactant_and_repeated_negative_product` in test_jacobian.py ("A + D -> C - D", "A -> C - 2 D - B", D at 0.1/1/3/100 eps: one switch factor per negative species + FD Jacobian match). Parent re-ran verification: 619 passed, ruff clean, mypy clean, golden_results.h5 sha256 unchanged.
- Parent review round 2: same reviewer re-reviewed the corrected diff: REVIEWED, no findings; P3-1 and P3-2 resolved. Open items (CLAUDE.md item 4, CHANGELOG, stale derive.json entries) are ticket 06 / final-report items. Accepted.

## Blockers / Challenges

None. Notes: (1) the bundled configs have no negative product, hence the toy and `negative_cell` environments for the switch tests and the unchanged golden numbers; (2) the Bash tool choked on heredocs containing triple quotes, so the larger edits were applied through small scripts kept in .scratch/tests/2026-11-10__derive-math-fixes/03-edit_solver.py and 03-edit_jac_tests.py (plus 03-solver_stats.py, 03-neg.toml, 03-conftest_mut.py).

## Commits

- 85168c7 fix: smooth depletion switch replaces the hard depletion block (ticket 03)

## Time Spent

16m (ticket-start.js to ticket-complete.js)
