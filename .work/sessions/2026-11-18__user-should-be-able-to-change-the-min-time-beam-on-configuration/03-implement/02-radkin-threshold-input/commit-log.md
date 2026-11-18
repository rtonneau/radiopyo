# Ticket 02: radkin-threshold-input

**Status:** ✅ Done

## Local Test Result

`cargo test --release` in C:/DEV/RUST/radkin: 302 passed, 0 failed (34 binaries); `cargo fmt --check` clean. `cargo clippy --release --all-targets -- -D warnings` fails on 3 `needless_range_loop` lints already on radkin main (src/irt/solver.rs:145, src/ssa/rdme_tau.rs:155/236), untouched here.

## Review Notes

- radkin commit 3461dae on branch `feat/configurable-min-beam-on-time` (from main 68613f3), pushed after the user's OK; PR https://github.com/rtonneau/radkin/pull/22. Second commit 5f7bba8 on the same PR (user-approved): the docs/literature paper-title paths (up to 252 chars) broke cargo's libgit2 checkout on Windows (MAX_PATH, 'path too long'), so they were renamed to short keys. Rev pinned in radiopyo (`rust/Cargo.toml`, `rust/Cargo.lock`, committed with this ticket): `5f7bba8`.
- Design choice: the threshold is a process-wide setting (`consts::set_min_beam_on_time` / `min_beam_on_time`, atomic), as radiopyo's global, not a cell field: beams are validated when built, before any cell exists. Both the beam ON time check and the short ON span check read it; spans <= T_ZERO stay analytic (`is_analytic_on_span`).
- Test: `tests/min_beam_on_time.rs` (own binary, one test, since it changes the global): invalid values, default / lower / below-T_ZERO / higher thresholds on the build check and the short span path.

## Blockers / Challenges

Pre-existing clippy failures on radkin main (newer clippy lint), left as is and noted in the radkin PR.

## Commits

- b6717c9 build: pin radkin 5f7bba8, configurable min beam ON time (ticket 02)

## Time Spent

1h 55m (ticket-start.js to ticket-complete.js)
