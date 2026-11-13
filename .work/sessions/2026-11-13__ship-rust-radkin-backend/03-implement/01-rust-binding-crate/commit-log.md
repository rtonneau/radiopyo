# Ticket 01: rust-binding-crate

**Status:** ✅ Done

## Local Test Result

`cd rust && cargo fmt --check && cargo clippy --release -- -D warnings`: clean (built against the local radkin through the git-ignored `rust/.cargo/config.toml` `[patch]`, i.e. radkin branch `feat/python-binding-hooks`).
`radiobio_env python -m maturin develop --release`: built `radiopyo_rs-0.7.0-cp311-abi3-win_amd64.whl` and installed it editable.
`python -I -c "import radiopyo_rs as r; print(r.__version__, r.RadkinRustError)"` prints `0.7.0 <class 'radiopyo_rs._core.RadkinRustError'>`.
Smoke test `.scratch/tests/2026-11-13__ship-rust-radkin-backend/01-smoke.py` (config_Labarbe, pulsed beam 1e6 Gy/s peak, 10 pulses, LET 0.201): `Cell.system()` returns 29 species (22 dynamic) and 93 reactions. `run("chunked", (0, 0.1))` returns t (13699,) and y (22, 13699) with 20 step dose rates. `run("run", (0, 1e-3))` returns y (22, 1304). An unknown option raises `RadkinRustError('value', "unknown radkin option 'bogus'")`.
radkin `feat/python-binding-hooks`: `cargo fmt --check`, `cargo clippy --all-targets --release -- -D warnings` and `cargo test --release` all pass.

## Review Notes

Base: 9ebed0e5b3cf9d0ce4eca486bb1a7faeda8f148e. Pre-existing changes: none besides the session files.
Deviations from the ticket text, chosen during implementation:
- The API is a frozen `Cell(config: dict, gvalues_csv, temperature, beam)` class with `system()` and `run(mode, t_span, options)` rather than a single `_core.run(config_toml, ...)`. The Python dict is walked into a `toml::Table` in Rust (no TOML writer in Python) and passed to the new `Config::from_table`. One cell can be checked by the parity guard once and run many times, and ticket 4's batch reuses the `execute()` helper.
- radiopyo's own Boscolo CSV text is passed in (`GValuesBoscolo::from_csv_str` + `UnitCell::from_config_with_gvalues`). radkin's `data_dir()` falls back to `CARGO_MANIFEST_DIR`, which breaks inside a wheel.
- There is no pyo3-log dependency: a small `log::Log` forwards records to `logging.getLogger("radiopyo.radkin")`.
- The `includes` and `beam` config sections are dropped, and the beam is always sent from Python's current state.
- The error kinds are parse/config/io/beam/short_beam_on/no_beam/solver_failed/gvalue/value (`RadkinRustError.args = (kind, msg)`).
The radkin additions are in PR https://github.com/rtonneau/radkin/pull/5 (additive, no behaviour change).
PR #5 was merged as 5e0d181a338323d39b61c62a0e7c740b5c1739ab. `rust/Cargo.toml` now pins `rev = "5e0d181"`, and `Cargo.lock` points at that commit. With the local `[patch]` disabled, fmt and clippy are clean, `maturin develop --release` succeeds and the smoke test gives identical output.

## Blockers / Challenges

- `cargo fetch` failed without `src/lib.rs`, so a placeholder was added first.
- pyo3 0.30 renamed the GIL calls (`Python::attach` / `py.detach`) and `FromPyObject` now has an associated `Error` (converted with `Into<PyErr>`).
- Waited on the merge of radkin PR #5 (done). cargo's built-in git client fails to authenticate on this machine; `CARGO_NET_GIT_FETCH_WITH_CLI=true` fixes it.

**Blocked (2026-11-13T11:02:21.548Z):** waiting for the user to merge radkin PR https://github.com/rtonneau/radkin/pull/5; then bump rust/Cargo.toml rev to the merge commit, rebuild without the local [patch], ticket-complete

## Commits

- d45a37b feat: radiopyo-rs PyO3 binding of radkin (ticket 01)

## Time Spent

1h 02m (ticket-start.js to ticket-complete.js)
