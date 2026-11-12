# Ticket 02: overshoot-tables

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q` -> 760 passed (equivalence test at rtol 1e-12 included); `ruff check radiopyo` clean; `mypy radiopyo` no issues.

`python examples/bench_overshoot_pulsed.py` (baseline in ticket 01: 22.26 s, _overshoot 9.14 s / 25.9%):

```
wall (best of 3): 18.05 s  all: 18.05 18.28 18.38
profiled total: 29.99 s
  derive                             165422 calls   22.89 s (76.3%)
  depletion_switch                   342652 calls    5.98 s (19.9%)
  _overshoot                          27035 calls    3.68 s (12.3%)
  _k_relaxation_terms                 38792 calls    3.31 s (11.0%)
  jacobian                             5904 calls    3.11 s (10.4%)
  mass_action_derivatives              5904 calls    0.88 s ( 2.9%)
  depletion_switch_derivatives        11808 calls    0.18 s ( 0.6%)
  _ez_relaxation_terms                38792 calls    0.01 s ( 0.0%)
```

Table build (config_Labarbe, 30 species): 0.9 ms, inside build_all_matrices 2.6 ms.

## Review Notes

Base: 85779c4b0984141817fcabb8d1c13fcf861a6151, clean tree.
- `KOvershootTable`, `EzOvershootTable`, `OvershootTables` (frozen, arrays locked) in sim_types; `SimMatrices.overshoot` is `field(init=False, repr=False, compare=False)` set in `__post_init__` via `object.__setattr__`; constructor unchanged (the test building SimMatrices by hand passes untouched).
- `_k/_ez_relaxation_terms` take the species table; `_relaxation_sum` takes weight/flag/negative_without_j. No `dataclasses.replace` left (import removed). Public helpers delegate to private raw-array kernels (`_mass_action`, `_mass_action_derivatives`, `_switch`, `_switch_derivatives`), signatures unchanged.
- The remaining overshoot cost is the per-column Python loop of `_switch`; ticket 03.

## Blockers / Challenges

A bash heredoc with nested triple quotes failed to parse; the edit script went to the scratch dir instead. Otherwise none.

## Commits

- 806b2b2 perf: precompute per-species overshoot tables on SimMatrices (ticket 02)

## Time Spent

5m (ticket-start.js to ticket-complete.js)
