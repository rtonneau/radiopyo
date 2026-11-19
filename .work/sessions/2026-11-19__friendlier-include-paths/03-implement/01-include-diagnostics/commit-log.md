# Ticket 01: include-diagnostics

**Status:** ✅ Done

## Local Test Result

`python -m pytest tests/parser -q` → 89 passed; `ruff check radiopyo` → All checks passed!; `mypy radiopyo` → Success: no issues found in 47 source files. Full suite: 962 passed in 57.5 s.

## Review Notes

Base: eb34761c8e6ed245519d4547f56d962556926588. Pre-existing changes not owned by this ticket: docs/issues/fenton/Example_Emballement.ipynb, docs/issues/fenton/Fenton_Fixed_Showcase_Rust.ipynb (user's, untouched).

Checked against the criteria: `radiopyo/parser/exceptions.py` (`ConfigError(RadioPyoError, ValueError)` with `message`/`origin`, `IncludeNotFoundError` also a `FileNotFoundError`, `ConfigSyntaxError`), exported from `radiopyo.parser` and `radiopyo`; list / single-string `includes` (`normalize_includes`); unknown table key or list suffix → `ConfigError` naming the file and the accepted forms (raised by the parser, `anchor_includes`, and by `resolve_includes` for hand-built dicts); missing include → places tried, include chain, suggestions (`include_suggestions`: exact bundled file name first, else difflib over bundled entries and local TOMLs); CWD fallback warns; circular include → `ConfigError`; `rp.bundled_configs()` sorted by entry. `resolve_includes` gained `source=` (used by `UnitCell.from_parser`) and `load_config(path)` was added. Inline-only review (no independent reviewer: localized parser change, diff assessed directly). Adjusted after seeing a real message: an exact bundled-name match now suggests only that file.

## Blockers / Challenges

`test_unknown_extension_is_skipped` asserted the old silent skip; renamed to `test_unknown_extension_is_an_error` (approved behavior change). Long inline heredocs broke the Bash tool; edits went through scratch scripts / Edit.

## Commits

- 9db6f04 feat: include diagnostics, list-form includes and bundled_configs (ticket 01)

## Time Spent

6m (ticket-start.js to ticket-complete.js)
