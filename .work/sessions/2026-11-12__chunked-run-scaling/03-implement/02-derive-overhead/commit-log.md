# Ticket 02: derive-overhead

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`: 768 passed in 37.02s; ruff "All checks passed!"; mypy "Success: no issues found in 42 source files".
`python examples/bench_overshoot_pulsed.py --repeat 2` (200 chunks, config_Labarbe), run wall per tqdm: before (stash, base 8554ac4) 14.7 s / 15.0 s, after 10.4 s / 10.8 s (about -29 %).
Micro benchmark (`.scratch/.../02-micro.py`, config_Labarbe state): derive 51.8 -> 33.2 µs per call, derive with a species below 0 106 -> 72-77 µs, jacobian 144 -> 131 µs. derive and jacobian outputs bit-identical to the base (np.array_equal on a positive and a negative state).

## Review Notes

Mode inline, base 8554ac4501df8a2d2f0291c3ac9fe2ca6ce2bf37.
- `StateLayout.to_full_molar`: 1-D state concatenated directly (no broadcast_to).
- `_mass_action`: 1-D path writes max(y, 0) into a preallocated (n+1) buffer with the padding 1.0 (no np.ones + concatenate), `np.multiply.reduce` (same reduction as np.prod).
- `reaction_rates`: the depletion switch is skipped when no reaction of the family has a negative product (`_switched`; x * 1.0 == x, so values unchanged); `_switch` 1-D path uses np.append instead of ones + concatenate.
- `derive`: `flatnonzero` only when a species is below 0.
- Golden tests and test_overshoot_equivalence pass unchanged; no golden data re-recorded. Self-review (localized hot-path change, bit-identical outputs checked).

## Blockers / Challenges

None.

## Commits

- f0d13fe perf: trim per-call allocations in derive (ticket 02)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
