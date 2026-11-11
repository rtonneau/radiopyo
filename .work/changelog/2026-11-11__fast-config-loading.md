---
bump: patch
floor: patch
---
### Changed

- Faster config loading: the reaction-string parser (`RPARSER`) remembers each string it has parsed and returns a copy, so a string is parsed by Lark only once per process. `UnitCell.from_toml(LABARBE_CONFIG_FILE)` goes from about 2.3 s to 0.65 s on the first load and 0.04 s on the next ones in the same kernel. Results are unchanged. New `ReactionParser.clear_cache()`.
