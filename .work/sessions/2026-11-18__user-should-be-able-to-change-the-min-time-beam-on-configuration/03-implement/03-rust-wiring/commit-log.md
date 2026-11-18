# Ticket 03: rust-wiring

**Status:** ✅ Done

## Local Test Result

`maturin develop --release` (radkin from GitHub rev 5f7bba8, `CARGO_NET_GIT_FETCH_WITH_CLI=true`): Installed radiopyo-rs-1.0.0. `pytest -q tests/rust tests/golden`: 70 passed, none skipped. Full `pytest -q`: 872 passed; ruff and mypy clean; `cargo fmt --check` in rust/ clean.

## Review Notes

- Pin: done in ticket 02's commit (`rev = "5f7bba8"`, no path dependency; the local `[patch]` used during development was removed).
- Binding: `Cell(..., min_beam_on_time=None)` sets radkin's process-wide value before the beam is built; `Cell.system()` reports `min_beam_on_time`.
- `rust.py`: `build_cell` passes `uc.MIN_BEAM_ON_TIME`; `MIN_BEAM_ON_TIME` left `_RADKIN_CONSTANTS`; `check_parity` compares it.
- The binding also maps the new `RadkinError::Checkpoint` (from the 145 radkin commits) to kind `"checkpoint"`. No parity drift from those commits: all Rust and golden tests pass.
- Tests `tests/rust/test_rust_min_beam_on_time.py` (6): default / 0.2 us / below T_ZERO last-pulse runs match Python, a 0.3 us configured pulse runs with 0.2 us, the parity guard flags a different value, radkin rejects a beam that the higher value makes too short.

## Blockers / Challenges

radkin main gained docs paths too long for Windows checkouts: fixed in radkin (ticket 02, commit 5f7bba8).

## Commits

- e43d526 feat: Rust backend receives the min beam ON time, parity guard compares it (ticket 03)

## Time Spent

0m (ticket-start.js to ticket-complete.js)
