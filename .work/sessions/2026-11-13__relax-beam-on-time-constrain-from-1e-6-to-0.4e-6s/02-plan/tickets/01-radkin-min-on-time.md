# Ticket 01: radkin-min-on-time

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] `src/consts.rs` MIN_BEAM_ON_TIME = 0.4e-6; tests and docs that use the 1e-6 boundary moved to the new value.
- [ ] `cargo test` passes; committed and pushed to origin of C:/DEV/RUST/radkin.

**Files to Touch:**
- `C:/DEV/RUST/radkin/src/consts.rs`
- `C:/DEV/RUST/radkin/src/physics/beam.rs` (tests)

**Verification Step:**

Run:
```bash
cd C:/DEV/RUST/radkin && cargo test
```

Expected:
all tests pass

**Notes:**

radkin is a separate repo; report its new commit SHA.
