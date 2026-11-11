# Session: fast-config-loading

**Date:** 2026-11-11T19:53:00.446Z
**Status:** Grill phase complete

## Problem Statement

`UnitCell.from_toml(LABARBE_CONFIG_FILE)` is slow (9.5 s cold / 3.3 s warm reported by the scout; about 2.3 s per load measured on 2026-11-11), and every new `UnitCell` in a notebook pays it again. Almost all of it is the Lark Earley reaction-string parser (about 24 ms per string): one load makes 295 `parse_reaction` calls for only 77 distinct strings.

## Context & Constraints

- **Current behavior:** `RPARSER` (the shared `ReactionParser` singleton in `radiopyo/parser/`) parses every string from scratch. `reaction_key` (`parser/utils.py`, used by `_merge_reaction_list`) parses every reaction at each include merge level, then the builders in `simulation/environment_build.py` (`_with_k_constants`, `_build_acid_base`, `_build_equilibrium`) parse them again; `parse_species` is called for H2O, pH species, radiolytic products and fixed species.
- **Pain point:** each string is parsed about 4 times per load and again on every load in the same kernel.
- **Dependencies:** results must be unchanged (golden tests in `tests/golden/`). Callers mutate the returned `GenericReaction` (`gen["constants"][...] = ...`) and attach the species objects to a `ListSpecies`, so cached results must never be shared.
- **Tech stack:** Python 3.11+, lark (Earley; the grammar is not LALR-compatible as written: reduce/reduce collisions `radical : molecule` vs `species : multiple molecule`).

## Success Metrics

- Loading LABARBE runs the real Lark parse exactly once per distinct stripped string (parse-count pytest); a second load in the same process runs no real parse.
- Returned objects are independent: mutating a result (constants, species) does not affect a later call for the same string.
- A string that fails to parse still raises on every call (failures are not cached).
- Full pytest (golden tests included, apart from the two known time-grid failures), `ruff check radiopyo` and `mypy radiopyo` stay at zero errors.

## Architecture & Approach

In-process memoization inside `ReactionParser`: a per-string cache (key = stripped string) of the parsed result for `parse_reaction` and `parse_species`, each call returning `copy.deepcopy` of the cached value (deepcopy about 0.2 ms vs about 24 ms per parse). Unbounded dict, plus a `clear_cache()` method (tests). Failed parses are not cached. The grammar and all callers (merge, builders) are left unchanged: the cache makes the repeated parses free, so the parsed reaction is not threaded from the merge to the builders.

## Assumptions & Trade-offs

- The parse is deterministic for a given string, so a cached result is identical to a fresh one.
- Not doing: an on-disk cache surviving kernel restarts (the first load in a process still pays about 77 parses, roughly 0.6 s per 77 strings), and an LALR grammar rewrite (bigger per-parse speedup, but the grammar is ambiguous for LALR and parse results could change). Both stay possible later.
- No timing-threshold test (machine-dependent, flaky); the guard is the parse count.

## Open Questions

None.

## Notes

Measured 2026-11-11: `import radiopyo` about 2.8 s (not in scope), one LABARBE load about 2.3 s with 295 `parse_reaction` calls / 77 distinct strings. Scout seed: `scout-reports/architecture-review-2026-11-11T13-37-44-572Z.html`.
