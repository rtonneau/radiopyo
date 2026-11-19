# Config errors name their file of origin

_ADR 0012, recorded 2026-11-19 (session 2026-11-19__friendlier-include-paths)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

With nested includes and field-by-field overrides, a bad config entry was reported without saying which file to edit. Every config error (TOML syntax, `[includes]`, include resolution, reaction strings, and the errors raised while `SimEnv` is built from the merged config) is now raised as a `ConfigError` (`RadioPyoError`, `ValueError`) naming the file of origin of the bad field, tracked field by field through the include merge (the include chain for nested files), with the original exception chained as `__cause__`. This changes the exception types callers see (a `KeyError` for a missing `k_value` becomes a `ConfigError`), which we accept for one uniform, catchable error; the origin table is kept outside the `ConfigDict`, so the Rust backend never sees it.
