# Ticket 01: memoize-reaction-parser

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q` → `710 passed in 43.73s` (the two golden time-grid tests listed as known failures in CLAUDE.md pass too); `ruff check radiopyo` → `All checks passed!`; `mypy radiopyo` → `Success: no issues found in 42 source files`. Output in `.scratch/tests/2026-11-11__fast-config-loading/01-pytest.txt`.

Timing (`measure_load.py` in the scratch dir, LABARBE): first load 2.35 s → 0.65 s, second load in the same process 2.37 s → 0.04 s; still 295 `parse_reaction` calls / 77 distinct strings, now 77 Lark parses on the first load and 0 on the second.

## Review Notes

Base: b5b399ea96053e5488549feb1528054344be1759 (no pre-existing changes besides gps state). Inline, localized change: no independent reviewer (one private method plus a test file; diff and checks assessed directly).
- Cache lives in `ReactionParser._parse` (key = stripped string), so `parse_reaction` and `parse_species` both use it; every call returns `copy.deepcopy`. A Lark error raises before the store, so failures are not cached. `clear_cache()` added.
- Tests (`tests/parser/test_reaction_parser_cache.py`): LABARBE load = one Lark parse per distinct string and none on a second load; spacing variants share a parse; mutating constants/reactants/products of a result does not leak; parse_species returns new objects; failures re-raise and are re-parsed; clear_cache forces a new parse.
All acceptance criteria met.

## Blockers / Challenges

A first test version mutated `species.label` (read-only property) then `species.mass` (getter raises: mass lookup is disabled); replaced by mutating the reactants/products dicts and an identity check. Note: CLAUDE.md still says two golden tests fail; they pass on main now.

## Commits

- a163a82 perf: memoize reaction-string parses in ReactionParser (ticket 01)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
