---
bump: minor
floor: minor
---
### Added

- `includes` accepts a list of files (`includes = ["radiopyo/advanced_water_radiolysis.toml", "my.toml"]`) or a single string; the parser is chosen from the file extension. The `[includes]` table keyed by extension still works.
- `rp.bundled_configs()` lists the configs shipped with the package, keyed by the `radiopyo/<file>` entry to write in `includes`.
- `rp.ConfigError` (a `ValueError`, under `RadioPyoError`) for every configuration error, with `rp.IncludeNotFoundError` (also a `FileNotFoundError`) and `rp.ConfigSyntaxError`. Its message starts with the file the error comes from (ADR 0012, `.work/adr/0012-config-error-origin.md`):
  - a missing include lists the paths tried, close matches among the neighbouring and bundled files, and a hint about backslashes in double-quoted Windows paths;
  - invalid TOML and invalid reaction strings name the file and the include chain (`main.toml -> a.toml`);
  - an entry that cannot be built (reaction, concentration, parameter, radiolytic yield, beam) names the entry, the file defining it and the file of each field redefined by another file (`a.toml (k_value from main.toml)`); `<dict>` for a config given as a Python dict.
- `load_config_with_origins` / `resolve_includes_with_origins` (`radiopyo.parser.file_parser`) return the merged config and its `ConfigOrigins`; `SimEnv.from_dict` and `UnitCell.from_config_dict` take an optional `origins=`. The `ConfigDict` itself is unchanged (it is still what the Rust backend receives).

### Changed

- Errors raised while building a simulation from a config (`SimEnv.from_dict`, `UnitCell.from_toml` / `from_config_dict`) are now `ConfigError` instead of the raw `KeyError`, `TypeError`, `AttributeError`, ...; the original error is kept as `__cause__`. The missing-section error is a `ConfigError` (still a `ValueError`). The reader helpers called directly (`kreaction_from_dict`, ...) still raise the original error.
- An unknown key of the `[includes]` table, an entry without extension or with an unknown extension is an error naming the file; it used to be skipped silently.
- An include found only relative to the current working directory (not next to the including file) logs a warning.
