# Ticket 02: port-tests

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q` -> `261 passed in 18.53s` (no collection errors). Also run with PYTHONHASHSEED=1 and 7: 261 passed each. `tests/golden` still passes.

## Review Notes

- Parent: base 1d430f924f5c4c6774ac99b27ea5be26ca05491e; pre-existing changes: gps session state only.

- AC1: `tests/chemistry` ported to BasicSpecies/SimSpecies, `stoi=`, `RadiolyticReaction.from_product`, `ListSpecies` taking raw species. 192 tests pass. Adaptations to the current contract: adding a species twice to a reaction raises KeyError (was increment); `SimSpecies` has no `__eq__` so comparisons use labels or identity; `ListSpecies.copy()` makes new SimSpecies objects; `BasicSpecies` has no `in` support (replaced by an eq-unsupported-type test); `BasicReaction.__contains__` rejects SimSpecies (pinned); `reactants_as_string` is now `reactants.as_string()`. Removed an empty syntax-broken stub at the end of test_k_reactions.py and added a small `raw_k_value` test.
- AC2: new `tests/parser/test_includes.py` (10 tests: k-reaction override, 2k override, nested includes, built-in include, missing file, unknown extension) and `tests/physics/test_timers.py` (Timer and MacroTimer with hand-computed states and transition times).
- AC3: `TODO(ticket 14)` on first-include-wins; `TODO(ticket 15)` on the walrus precedence in `micro_pulses_per_train` and `duty_cycle` and on `copy()` dropping `start_time` (both timers). Quirk with no ticket: `time_till_next_OFF` gives NaN for infinite periods (`inf*0`), pinned and commented. Not covered: enzymatic doubling, acid-base/constant override (ticket 14), `tot_pulses` start_time (ticket 15, beam-level) and species order (ticket 16), since I could not pin them with confidence from the timer/include scope.
- AC4: full suite passes. No file under `radiopyo/` changed.

- Parent acceptance: reviewer skipped — test-only, no package code touched; parent re-ran full suite (261 passed). Model sonnet (requested medium).

## Blockers / Challenges

None

## Commits

- 9424c85 test: port chemistry tests to current API, pin include and timer behaviour (ticket 02)

## Time Spent

7m (ticket-start.js to ticket-complete.js)
