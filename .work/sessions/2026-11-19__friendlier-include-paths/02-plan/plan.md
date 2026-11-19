# Implementation Plan

**Session:** friendlier-include-paths
**Date:** 2026-11-19T13:15:20.342Z
**Estimated effort:** 1 day

## Strategy

Build the error vocabulary and the include diagnostics first (the most visible user win), then make every per-file parse error name its file, then carry field-level provenance through the include merge so SimEnv build errors name their file too, and finish with docs and changelog. Each ticket leaves pytest, ruff and mypy green.

## Tickets Overview

- **Ticket 1:** `ConfigError` hierarchy, list-form `[includes]`, strict keys/suffixes, `IncludeNotFoundError` diagnostics (tried locations, include chain, suggestions), CWD-fallback warning, `rp.bundled_configs()`.
- **Ticket 2:** per-file parse errors name their file: TOML syntax errors (`ConfigSyntaxError`, Windows-path hint) and reaction strings checked at read time.
- **Ticket 3:** field-level origin tracking through the include merge, and SimEnv/beam build errors re-raised as `ConfigError` naming the origin of the faulty field.
- **Ticket 4:** docs (`docs/configuration.md`), CHANGELOG, CLAUDE.md parsing paragraph.

## Sequencing Rationale

Ticket 1 creates the exceptions module and include-chain plumbing that tickets 2 and 3 reuse. Ticket 2 is independent of the merge and small. Ticket 3 depends on the origin of each file being known at read time (tickets 1–2). Docs last, describing the final behavior.

## Risks & Mitigation

- **Risk:** provenance bookkeeping changes the merged ConfigDict or the Rust payload → **Mitigation:** origins live in a separate object, never inside the ConfigDict; a test asserts `radkin_config` and golden outputs are unchanged.
- **Risk:** wrapping build errors as `ConfigError` breaks tests or callers expecting `KeyError`/`TypeError` → **Mitigation:** update the affected tests deliberately, list the change in the CHANGELOG (ADR 0012).
- **Risk:** reaction-string checks at read time double the parse cost → **Mitigation:** `RPARSER` caches parses (test_reaction_parser_cache); measure nothing more than the existing suite time.

## Assumptions

- Only TOML parser exists; the list form infers the parser from the file suffix via `list_available_parsers()`.
- `combine_config` / `resolve_includes` keep their current signatures and results; provenance is exposed through new entry points (e.g. a function returning the config and its origins) used by `UnitCell.from_parser`.
- A config dict passed directly to `from_config_dict` has origin `<dict>`.
