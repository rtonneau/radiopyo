---
bump: minor
floor: minor
---
### Changed

- `MIN_BEAM_ON_TIME` lowered from 1e-6 s to 0.4e-6 s (still 4x `T_ZERO` = 1e-7 s) so real UHDR pulses of 0.4 to 1 µs can be simulated: a configured beam ON time in that range no longer raises `ShortBeamOnError`, and a run-time ON span of that length is integrated instead of added analytically. The Rust backend (radkin rev `088f06b`) uses the same threshold.
