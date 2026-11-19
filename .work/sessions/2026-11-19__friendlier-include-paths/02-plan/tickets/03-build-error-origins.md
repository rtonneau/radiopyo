# Ticket 03: build-error-origins

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] An origin table (e.g. `ConfigOrigins` in `radiopyo/parser/`) records, for each entry (parameters/concentrations/radiolytic by name, reactions by section + `reaction_key`, beams by type/index) and each of its fields, the file that wrote it; it is filled while resolving includes and updated following `combine_config` precedence (field-by-field override, equilibrium full redefinition drops the old constants' origins).
- [ ] The ConfigDict produced is unchanged (same as today's `resolve_includes` output); `UnitCell.config` and `rust.radkin_config` never contain origins.
- [ ] `UnitCell.from_parser` / `from_toml` use the origins; `from_config_dict(config, init, origins=None)` and `SimEnv.from_dict(config, ..., origins=None)` accept them; without origins the origin is `<dict>`.
- [ ] An error raised while building a reaction, species/concentration, radiolytic yield or beam from the merged config is re-raised as `ConfigError` (original as `__cause__`) whose message names the entry and its origins, e.g. "[[reactions.k_reaction]] 'A + B -> C' (defined in a.toml; k_value from main.toml): <original message>".
- [ ] Tests: field overridden in main file with a bad value names main.toml; bad field in an included file names it; dict config names `<dict>`; full suite incl. golden tests passes.

**Files to Touch:**
- `radiopyo/parser/origins.py`
- `radiopyo/parser/utils.py`
- `radiopyo/parser/file_parser.py`
- `radiopyo/simulation/environment_build.py`
- `radiopyo/simulation/unit_cell.py`
- `radiopyo/physics/beam_collection.py`
- `tests/parser/test_config_origins.py`

**Verification Step:**

Run:
```bash
C:/Users/rtonneau/AppData/Local/miniconda3/envs/radiobio_env/python.exe -m pytest -q && C:/Users/rtonneau/AppData/Local/miniconda3/envs/radiobio_env/python.exe -m ruff check radiopyo && C:/Users/rtonneau/AppData/Local/miniconda3/envs/radiobio_env/python.exe -m mypy radiopyo
```

Expected:
Full suite passes (about 750+ tests, golden included); ruff and mypy zero errors.

**Notes:**

Keep `combine_config` and `resolve_includes` signatures; add a sibling that returns `(ConfigDict, ConfigOrigins)` and have `resolve_includes` delegate to it. Wrap errors at the per-entry boundary of `reactions_from_dict` / `_REACTION_BUILDERS` and the species/beam readers, not around the whole build. Do not wrap `ConfigError` twice. ADR 0012.
