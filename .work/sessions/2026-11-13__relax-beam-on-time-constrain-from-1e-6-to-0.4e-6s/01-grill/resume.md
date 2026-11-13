# Session: Relax beam ON time constrain from 1e-6 to 0.4e-6s

**Date:** 2026-11-13T14:41:44.139Z
**Status:** Grill phase complete

## Problem Statement

`MIN_BEAM_ON_TIME` = 1e-6 s rejects real UHDR pulses between 0.4 and 1 µs (e.g. linac pulses). Lower it to 0.4e-6 s in both backends.

## Context & Constraints

- **Current behavior:** a beam whose configured ON time is below 1e-6 s raises `ShortBeamOnError` (ADR 0008); run-time short ON spans below it are added analytically. radkin hard-codes the same 1e-6 (`src/consts.rs`) and `radiopyo/rust.py` `_RADKIN_CONSTANTS` refuses a UnitCell with a different value.
- **Pain point:** sub-µs real pulses cannot be simulated.
- **Dependencies:** keep Python/Rust parity; radkin is a separate repo (C:/DEV/RUST/radkin, GitHub rtonneau/radkin) pinned by git rev in `rust/Cargo.toml`.
- **Tech stack:** Python (radiopyo), Rust (radkin, radiopyo-rs PyO3/maturin).

## Success Metrics

- `MIN_BEAM_ON_TIME == 0.4e-6` in radiopyo (`physics/beam.py`, `UnitCell`), in radkin `consts.rs` and in `_RADKIN_CONSTANTS`.
- A 0.5 µs pulse builds and runs; a 0.3 µs pulse raises `ShortBeamOnError`; boundary tests in `tests/physics/test_min_beam_on_time.py` moved to the new value.
- pytest, ruff and mypy stay green; radkin `cargo test` green; rust tests pass with the rebuilt extension.

## Architecture & Approach

Change the constant in radiopyo and update docstrings (unit_cell.py, exceptions.py, ode_result.py), README and CHANGELOG. Change radkin `consts.rs` (and its tests/docs), commit and push radkin (confirm with the user before pushing), bump the rev in `rust/Cargo.toml`, update `_RADKIN_CONSTANTS`. Amend ADR 0008 with a dated note: threshold now 0.4e-6 s (still 4x T_ZERO = 1e-7 s) because real UHDR machines deliver 0.4–1 µs pulses.

## Assumptions & Trade-offs

Not making the threshold configurable. Golden data not re-recorded unless a bundled case has an ON time between 0.4 and 1 µs. Spans between 0.4 and 1 µs are resolved with fewer solver steps; accepted.

## Open Questions

N/A

## Notes

ADR handling: amend 0008, no new ADR.
