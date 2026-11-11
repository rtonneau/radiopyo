# Ticket 01: bench-and-reference

**Status:** ✅ Done

## Local Test Result

`python -m pytest tests/simulation/test_overshoot_equivalence.py -q` -> 8 passed.
Full `python -m pytest -q` -> 760 passed; `ruff check radiopyo` clean; `mypy radiopyo` no issues.

Baseline `python examples/bench_overshoot_pulsed.py` (config_Labarbe, 200 chunks):

```
wall (best of 3): 22.26 s  all: 22.26 22.72 23.76
profiled total: 35.27 s
  derive                             165422 calls   27.81 s (78.8%)
  _overshoot                          27035 calls    9.14 s (25.9%)
  _k_relaxation_terms                 38792 calls    8.35 s (23.7%)
  depletion_switch                   381444 calls    6.54 s (18.6%)
  jacobian                             5904 calls    3.56 s (10.1%)
  mass_action_derivatives             11314 calls    1.57 s ( 4.4%)
  _ez_relaxation_terms                38792 calls    0.32 s ( 0.9%)
  depletion_switch_derivatives        14513 calls    0.21 s ( 0.6%)
```

## Review Notes

Base: a25fd9261993ec4051fdddf067f231e611a8eb59, clean tree (only this ticket's files are new).
- Bench follows the banner layout; prints wall, chunks, cProfile cumulative times and shares.
- The reference copy in the test is generated verbatim from ode_functions_solver.py at a25fd92 (DEPLETION_LIMIT .. derive, plus ReactionRates). Fixtures: a toy config with negative products, 2 A, 3-reactant, fixed reactant, two enzymatic reactions (one with a negative product), and config_Labarbe; states with negatives, exact zeros and depletion-scale values; also a 2-state batch for depletion_switch.

## Blockers / Challenges

The pure-water config has no LET, so `set_beam` without LET fails (KeyError ''); the bench uses config_Labarbe, which reproduces the scout scenario (200 chunks). Wall time on this machine is 22 s (scout reported 7.5 s, other machine/config); overshoot share 26%, depletion_switch 19%.

## Commits

- 4729c77 test: overshoot equivalence reference and pulsed bench (ticket 01)

## Time Spent

6m (ticket-start.js to ticket-complete.js)
