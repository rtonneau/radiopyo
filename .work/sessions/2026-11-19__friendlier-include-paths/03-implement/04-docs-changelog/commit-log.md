# Ticket 04: docs-changelog

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q` → 989 passed (docs-only ticket, suite unchanged).

## Review Notes

docs/configuration.md [includes]: list/string/table forms, unknown key/extension is an error, CWD fallback warning, `rp.bundled_configs()` (output pasted from a real run, includes the IBA data file), Windows-path tip, "Errors name their file" with the real missing-include message and an origin example. CHANGELOG: new Unreleased section (Added: list form, bundled_configs, ConfigError hierarchy and origins, *_with_origins API; Changed: build errors are ConfigError with __cause__, unknown include key/extension, CWD warning), ADR 0012 linked. CLAUDE.md parsing paragraph mentions normalize_includes, CWD warning, bundled_configs, ConfigError, ConfigOrigins and config_entry. The user's pre-existing docs/configuration.md edits were left in place (file had no unrelated diff before ticket 1; fenton notebooks untouched).

## Blockers / Challenges

None.

## Commits

- d076781 docs: include forms, bundled_configs and config error origins (ticket 04)

## Time Spent

unknown (no ticket_started event)
