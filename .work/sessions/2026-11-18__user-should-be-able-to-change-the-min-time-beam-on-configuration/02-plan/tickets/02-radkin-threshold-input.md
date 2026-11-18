# Ticket 02: radkin-threshold-input

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] In radkin (C:\DEV\RUST\radkin, branch from current main), the minimum beam ON time is an input of the cell (default `MIN_BEAM_ON_TIME` = 0.4e-6 s) used by the beam ON time check and the short ON span check, with floor T_ZERO for the short span check as in radiopyo.
- [ ] `cargo test` passes; a test covers a non-default value (build check and short span path).
- [ ] Pushed to GitHub (branch + PR in rtonneau/radkin) only after the user's explicit OK; the pushed rev is recorded for ticket 3.

**Files to Touch:**
- `C:/DEV/RUST/radkin/src/consts.rs`
- `C:/DEV/RUST/radkin/src/physics/beam.rs`
- `C:/DEV/RUST/radkin/src/unit_cell.rs`
- `C:/DEV/RUST/radkin/src/config` (or wherever the cell config is read)
- `C:/DEV/RUST/radkin/tests/`

**Verification Step:**

Run:
```bash
cd C:/DEV/RUST/radkin && cargo test
```

Expected:
`test result: ok` for every test binary, 0 failed.

**Notes:**

radkin has its own gps sessions/conventions: follow its AGENTS/CLAUDE files. Ask before pushing.
