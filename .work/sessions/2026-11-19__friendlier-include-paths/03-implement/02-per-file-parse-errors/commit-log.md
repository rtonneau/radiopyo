# Ticket 02: per-file-parse-errors

**Status:** ✅ Done

## Local Test Result

`python -m pytest tests/parser -q` → 98 passed; `ruff check radiopyo` → All checks passed!; `mypy radiopyo` → Success: no issues found in 47 source files. Full suite before the new tests: 962 passed.

## Review Notes

Base: f55fd61858275518384bfca12509ce40106750d8. Pre-existing user changes (fenton notebooks) untouched.

Checked against the criteria: `TOMLFileParser.parse` reads bytes and wraps `TOMLDecodeError` (and non-UTF-8 input) in `ConfigSyntaxError` naming the file, the faulty line with a caret, the original error chained; the backslash hint (`BACKSLASH_HINT`) is added for escape errors on a line holding a backslash. `FileParser.check_reactions` parses every reaction string of the file at read time (sections `REACTION_LISTS`, single-table form too) and raises `ConfigError` with section and index; lark's token list is cut to its first line plus a syntax example. `resolve_includes` re-chains a parser error of an included file (`ConfigError.with_chain`) so it names `main.toml -> a.toml`. Extra: a not-found include whose entry contains a tab/new line (valid TOML escapes from "C:\temp\new.toml") gets the same hint. Inline review only (localized parser change).

## Blockers / Challenges

The reaction grammar accepts an empty side (`A ->`, `-> B`), so those are not errors; the tests use strings without an arrow. Bash heredocs collapse `\`, so edits with escapes went through Edit.

## Commits

- 804c558 feat: TOML syntax and reaction string errors name their file (ticket 02)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
