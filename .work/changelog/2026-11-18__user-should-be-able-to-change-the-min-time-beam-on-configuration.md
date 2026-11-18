---
bump: minor
floor: minor
---
### Added
- `rp.set_min_beam_on_time(value)` changes the shortest beam ON time (default 0.4e-6 s) for every beam built and run started afterwards; `None` restores the default and `rp.min_beam_on_time()` reads it. `UnitCell.MIN_BEAM_ON_TIME` is now a read-only view of it.
- The Rust backend uses the same value (radkin pinned to rev 5f7bba8), and the parity guard checks it.

### Changed
- A run-time beam ON span not longer than `T_ZERO` (1e-7 s) is always added analytically, whatever the setting; a setting below `T_ZERO` logs a warning.
