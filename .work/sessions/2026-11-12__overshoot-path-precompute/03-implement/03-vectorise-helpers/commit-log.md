# Ticket 03: vectorise-helpers

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q` -> 760 passed; `ruff check radiopyo` clean; `mypy radiopyo` no issues.

`python examples/bench_overshoot_pulsed.py`:

```
wall (best of 3): 14.85 s  all: 15.24 14.85 15.24
profiled total: 22.97 s
  derive                             165422 calls   17.23 s (75.0%)
  _overshoot                          27035 calls    2.75 s (12.0%)
  _k_relaxation_terms                 38792 calls    2.40 s (10.4%)
  jacobian                             5904 calls    1.98 s ( 8.6%)
  depletion_switch                   342652 calls    1.54 s ( 6.7%)
  mass_action_derivatives              5904 calls    0.48 s ( 2.1%)
  depletion_switch_derivatives        11808 calls    0.04 s ( 0.2%)
  _ez_relaxation_terms                38792 calls    0.01 s ( 0.0%)
```

Before / after over the session (wall, best of 3): 22.26 s (baseline) -> 18.05 s (ticket 02) -> 14.85 s (-33%). _overshoot 9.14 s -> 3.68 s -> 2.75 s; depletion_switch 6.54 s -> 5.98 s -> 1.54 s.

## Review Notes

Base: cf4f9db86b391989b98c6b533ab48d088b7a40b1, clean tree.
- No `np.delete` and no per-reaction/per-column loop left in the solver (grep). `depletion_switch` is a gather + product over a padded index of the negative products (`pad_slots`, derived `negative_index` init=False field on KMatrices/EzMatrices; overshoot tables store the padded index without j). The derivatives use prefix/suffix cumulative products (`_others_product`) and a bincount scatter (`_scatter_slots`); same multiplication and summation order as the reference for up to 3 slots.
- Equivalence test (rtol 1e-12) covers 2 negative products in one reaction, 2 A, 3 reactants, exact zeros, and a 2-state batch for depletion_switch.
- CHANGELOG: the fragment is written by /gps finish (changelog-prepare.js), so no CHANGELOG file in this ticket; CLAUDE.md architecture note updated instead.

## Blockers / Challenges

mypy: bincount returns an untyped float array (astype float64, copy=False) and `_relaxation_sum` took the old bool type; fixed.

## Commits

- 0cb3ab9 perf: vectorise depletion switch and mass-action derivatives (ticket 03)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
