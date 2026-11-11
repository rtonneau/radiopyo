# Ticket 01: memoize-reaction-parser

**Model:** inherit
**Effort:** inherit

**Acceptance Criteria:**
- [ ] `ReactionParser` caches the result of the Lark parse per stripped string; `parse_reaction` and `parse_species` return deep copies.
- [ ] A string that fails to parse is not cached and raises on every call.
- [ ] `ReactionParser.clear_cache()` empties the cache.
- [ ] Test: loading LABARBE runs the real parse once per distinct stripped string, and a second load in the same process runs none.
- [ ] Test: mutating a returned reaction (constants, a species) does not affect the next result for the same string.
- [ ] `ruff check radiopyo` and `mypy radiopyo` report zero errors; full pytest passes apart from the two known golden time-grid failures.

**Files to Touch:**
- `radiopyo/parser/reaction_parser.py`
- `tests/parser/test_reaction_parser_cache.py`
- `CHANGELOG.md` (or the gps changelog fragment)

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All tests pass except `test_beam_matrices_match_golden` and `test_chunked_pulsed_run_matches_golden` (known time-grid failures); ruff and mypy report no errors.

**Notes:**

Callers write into `gen["constants"]` and attach species to `ListSpecies`, so the cached value must never be handed out. Keep the module banner layout. Measured before: one LABARBE load = 295 `parse_reaction` calls for 77 distinct strings (about 2.3 s).
