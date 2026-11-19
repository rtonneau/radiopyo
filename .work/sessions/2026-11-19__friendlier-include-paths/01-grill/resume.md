# Session: friendlier-include-paths

**Date:** 2026-11-19T13:04:01.254Z
**Status:** Grill phase complete

## Problem Statement

Users get the `[includes]` paths of their TOML configs wrong (relative vs absolute vs bundled `radiopyo/` entries) and the errors do not help them fix it. More generally, a config error never says which file it comes from, which matters with nested includes and field-by-field overrides. Make include resolution friendlier and make every config error (parse time and SimEnv build time) name its file of origin.

## Context & Constraints

- **Current behavior:** `[includes]` is a table keyed by extension (`toml = [...]`). `radiopyo/<file>` maps to `radiopyo/data/<file>` (`include_path`); relative entries are anchored to the including file (`FileParser.anchor_includes`), else silently fall back to the CWD; absolute paths used as written. A missing file raises a bare `FileNotFoundError(path)`. An unknown extension key only logs a warning and the include is dropped. A top-level `includes = [...]` list crashes with `AttributeError`. TOML decode errors (e.g. `"C:\Users\..."` escapes) are raw tomllib errors with no file name. Reaction strings are only parsed when the merged config is turned into a `SimEnv`; build errors (KeyError, ValueError, TypeError, lark errors) carry no file information.
- **Pain point:** users cannot tell where an include was looked for, nor which file holds a bad entry; some mistakes are silent (dropped include, CWD fallback).
- **Dependencies:** keep the `radiopyo/` prefix and all existing resolution rules; keep precedence and field-by-field merge semantics (`combine_config`); the Rust backend receives `UnitCell.config` (merged ConfigDict) and must not see any provenance data; old configs and the HDF5 golden file must keep loading; ruff/mypy at zero; golden tests unchanged.
- **Tech stack:** Python 3.11 (`tomllib`, `difflib`), Lark reaction grammar (`RPARSER`), pytest.

## Success Metrics

- A missing include raises `IncludeNotFoundError` (subclass of `ConfigError` and `FileNotFoundError`) naming the including file (and include chain), the entry as written, every location tried, and close matches; `radiopyo/data/x.toml` or a bare bundled name gets "did you mean `radiopyo/x.toml`?".
- `includes = ["a.toml", ...]` works (parser chosen from the suffix); an unknown extension key or suffix raises a `ConfigError`; a relative include found only via the CWD logs a WARNING.
- A TOML syntax error raises `ConfigSyntaxError` naming the file, with a single-quote / forward-slash hint on invalid escapes; a bad reaction string raises a `ConfigError` naming its file at read time.
- An error raised while building `SimEnv` from a merged config is a `ConfigError` (original as `__cause__`) naming the file of origin of the faulty field, field by field (e.g. "reaction 'A + B -> C' (defined in a.toml; k_value from main.toml)"); a dict-built config reports `<dict>`.
- `rp.bundled_configs()` returns `{include entry: absolute path}` for every TOML under `radiopyo/data/`.
- Full pytest suite, ruff and mypy pass; golden outputs unchanged; Rust config payload unchanged.

## Architecture & Approach

- New `radiopyo/parser/exceptions.py`: `ConfigError(RadioPyoError, ValueError)` with an `origin` attribute, `IncludeNotFoundError(ConfigError, FileNotFoundError)`, `ConfigSyntaxError(ConfigError)`.
- Include resolution (`parser/utils.py`, `file_parser.py`): normalize the section (table or list), validate keys/suffixes, resolve each entry recording the locations tried, warn on the CWD fallback, build the diagnostic with difflib suggestions over `bundled_configs()` and the files next to the including file; the include chain is passed down for nested errors.
- TOML parser wraps `tomllib.TOMLDecodeError`; each file's reaction strings are checked with `RPARSER` at read time.
- Provenance: an origin table, kept outside the `ConfigDict` (returned alongside it / held by the parse result and `UnitCell`, never sent to radkin), maps each entry (section, reaction key or name) and field to the file that wrote it; `combine_config` updates it field by field following the merge precedence.
- `SimEnv` building wraps errors per entry (reaction readers, species, beams) into a `ConfigError` with the origin from the table.
- `rp.bundled_configs()` exposed at package level next to `PURE_WATER_CONFIG_FILE`.
- Docs: `docs/configuration.md` (list form, where files are looked for, error examples, Windows path tip, `bundled_configs`), CHANGELOG (exception types change), ADR 0012, glossary terms *Include entry* and *Bundled config*.

## Assumptions & Trade-offs

- Only one bundled-config spelling (`radiopyo/`); aliases are not resolved silently, only suggested in errors.
- CWD fallback kept (warning only), not deprecated.
- Build-time errors change type to `ConfigError` (a `ValueError`): callers catching `KeyError`/`TypeError` from config building may break; accepted for one uniform error (ADR 0012), noted in the CHANGELOG.
- Field-level provenance adds bookkeeping to the merge; the merge result itself (ConfigDict) is unchanged.

## Open Questions

None.

## Notes

- Glossary terms added: *Include entry*, *Bundled config*. ADR 0012 records the error-origin decision.
- `radiopyo/data/IBA/20240311/20240311_DATA_AERIAL.toml` is also a bundled TOML; `bundled_configs()` lists data TOMLs recursively.
- The user asked mid-grill that every parsing error contain the file it originates from; the scope was then extended to build-time errors.
