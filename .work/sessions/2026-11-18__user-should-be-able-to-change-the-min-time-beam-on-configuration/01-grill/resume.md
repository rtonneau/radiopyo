# Session: User should be able to change the min time BEAM ON configuration

**Date:** 2026-11-18T14:25:43.307Z
**Status:** Grill phase complete

## Problem Statement

The minimum beam ON time (`MIN_BEAM_ON_TIME` = 0.4e-6 s) is a hard-coded constant in radiopyo and in radkin: a user cannot change it, e.g. to model shorter pulses or to be stricter.

## Context & Constraints

- **Current behavior:** `radiopyo/physics/beam.py` `MIN_BEAM_ON_TIME` is read directly by `ParticleBeam._check_on_times` (configured ON time below it raises `ShortBeamOnError` at build time). `UnitCell.MIN_BEAM_ON_TIME` (class attribute, same value) drives the run-time short ON span path (`is_short_on_span`, `_short_on_span_step`: yield added analytically, ADR 0008). radkin hard-codes it in `src/consts.rs` (used in `physics/beam.rs` and `unit_cell.rs`); `radiopyo/rust.py` `_RADKIN_CONSTANTS` refuses a cell with another value.
- **Pain point:** no supported way to change it; overriding the UnitCell attribute only affects the run-time path, not the beam check, and breaks the Rust backend.
- **Dependencies:** the default stays 0.4e-6 s and the `MIN_BEAM_ON_TIME` constant is kept as that default (no result change by default); T_ZERO = 1e-7 s is unchanged; ADR 0008 holds; radkin is pinned by git `rev` in `rust/Cargo.toml` (no path dependency committed); version files stay in lockstep.
- **Tech stack:** Python 3.11, scipy LSODA, pytest/ruff/mypy; Rust radkin + PyO3 binding `radiopyo-rs` (maturin).

## Success Metrics

- `rp.set_min_beam_on_time(x)` (any positive finite x, else `ValueError`), `rp.set_min_beam_on_time(None)` restores 0.4e-6 s, `rp.min_beam_on_time()` reads it; `UnitCell.MIN_BEAM_ON_TIME` reports the same value (read-only view).
- A beam built after the change is checked against the new value (`ShortBeamOnError` threshold moves); a run uses the value current at run time for short ON spans; a value < T_ZERO logs a warning and every span <= T_ZERO still takes the analytic path (no solver failure).
- radkin accepts the threshold as an input (default 0.4e-6 s) for the beam check and the short-span check; `backend="rust"` with a non-default value gives the same result as the Python backend (parity test); full pytest, ruff, mypy green; default results unchanged (golden tests).

## Architecture & Approach

Global setting in `radiopyo/physics/beam.py` (`set_min_beam_on_time`, `min_beam_on_time`, exported from `radiopyo`); `_check_on_times` and `UnitCell` read it (`UnitCell.MIN_BEAM_ON_TIME` becomes a read-only view; no per-cell override). Run-time short span test uses max(setting, T_ZERO) as the integrability floor so spans <= T_ZERO are always analytic. radkin (`C:\DEV\RUST\radkin`): threshold becomes an input of the cell/beam build (default 0.4e-6), pushed to GitHub only after the user's OK; radiopyo pins the new rev, `radiopyo/rust.py` sends the current value, drops it from `_RADKIN_CONSTANTS` and covers it in the parity check. Docs: ADR 0008 amendment, glossary, README, CLAUDE.md, CHANGELOG.

## Assumptions & Trade-offs

Global (process-wide) state, chosen over a per-cell or TOML setting: simple, but cells in one process share it; beams already built are not re-checked when the value changes. Any positive value is accepted (values between T_ZERO and 0.4e-6 s give poorly resolved chemistry, the user's choice). Not done: TOML config entry, per-cell override, persisting the value in HDF5 results.

## Open Questions

N/A

## Notes

Sized Standard (public API, two repos). Pushing radkin to GitHub needs the user's explicit OK during ship. `rust.run_batch` runs cells in worker threads: the value is sent per request, so the global is read on the calling thread.
