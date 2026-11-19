# Ticket 03: build-error-origins

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q` → 989 passed (golden and tests/rust included, none skipped); `ruff check radiopyo` → All checks passed!; `mypy radiopyo` → Success: no issues found in 48 source files.

## Review Notes

Base: f5bad54c6abcabc29f7434d9286620a084702ca2. Pre-existing user changes (fenton notebooks) untouched.

New `radiopyo/parser/origins.py`: `ConfigOrigins` (entry -> `EntryOrigin(defined, fields, kind)`; `of(config, origin)` for one file, `combined(other)` mirrors combine_config: field by field for reactions and beams, full replace for name = value entries, equilibrium with 2+ constants drops the old constants' origins, a beam of another type or root's 'default' beam is replaced). `config_entry(what, origin)` wraps any non-ConfigError into `ConfigError("<entry>: <message>", origin=...)` with the original as `__cause__`; the origin is computed lazily (on error only). `resolve_includes` delegates to `resolve_includes_with_origins` (signatures kept), also `load_config_with_origins`; `combine_config` gets a keyword-only `beam_origin` (beams are built during the merge, so beam errors are located there). Wrapped per entry: reactions (`reactions_from_dict(..., origins=)`), fixed species, initial concentrations, pH, radiolytic yields, beams (`BeamCollection.from_config_dict/merge/add_beam_from_dict(origin_of=)`, lazy import of origins to avoid the parser <-> physics import cycle). `SimEnv.from_dict(..., origins=None)` and `UnitCell.from_config_dict(..., origins=None)`; without origins the origin is `<dict>`; direct calls of the reader helpers without origins raise the original error as before. The missing-section error is now a ConfigError (still a ValueError). ConfigDict and `UnitCell.config` unchanged (tested). Inline review: checked import cycles from each entry module, no double wrapping, Rust tests pass.

## Blockers / Challenges

Type change (approved in the grill): building errors from SimEnv.from_dict / UnitCell are now ConfigError with the original as __cause__; 3 tests expecting AttributeError/KeyError updated (test_config_warnings, test_env_queries). To list in the CHANGELOG (ticket 04).

## Commits

- b48c794 feat: build errors name the config file of the faulty entry (ticket 03)

## Time Spent

12m (ticket-start.js to ticket-complete.js)
