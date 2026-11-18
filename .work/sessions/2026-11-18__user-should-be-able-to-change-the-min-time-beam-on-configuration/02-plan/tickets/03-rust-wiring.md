# Ticket 03: rust-wiring

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] `rust/Cargo.toml` pins the radkin rev from ticket 2 (no path dependency committed).
- [ ] The binding and `radiopyo/rust.py` send the current minimum beam ON time; `MIN_BEAM_ON_TIME` leaves `_RADKIN_CONSTANTS` and `check_parity` compares it.
- [ ] A Rust test with a non-default value (e.g. 0.2e-6 s) matches the Python backend; all `tests/rust/` and golden tests pass with the rebuilt extension.

**Files to Touch:**
- `rust/Cargo.toml`
- `rust/Cargo.lock`
- `rust/src/lib.rs`
- `radiopyo/rust.py`
- `tests/rust/test_rust_min_beam_on_time.py`

**Verification Step:**

Run:
```bash
cd rust && python -m maturin develop --release && cd .. && python -m pytest -q tests/rust tests/golden
```

Expected:
Extension builds; all Rust and golden tests pass (none skipped for a missing extension).

**Notes:**

`CONDA_PREFIX=<env>` for maturin; `CARGO_NET_GIT_FETCH_WITH_CLI=true` if cargo's git auth fails. The 145 radkin commits since 088f06b come along: any parity drift stops the ticket and goes to the user.
