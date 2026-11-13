# Ticket 05: ci-and-release-wheels

**Status:** ✅ Done

## Local Test Result

`pip install -e .[dev]`: resolves. `cargo fmt --check` and `cargo clippy --release -D warnings`: clean. `maturin build --release`: built `rust/target/wheels/radiopyo_rs-0.7.0-cp311-abi3-win_amd64.whl`. `pytest -q`: 852 passed. Both workflow files parse as YAML. The CI jobs and wheels.yml have not run yet: they run on the PR, or on a tag or dispatch.

## Review Notes

Base 64d4e33.
- tests.yml: the new `test-rust` job (Ubuntu and Windows) runs rust-toolchain, rust-cache, cargo fmt/clippy, a maturin release wheel and pip install, asserts `radiopyo.rust.is_available()`, then runs the full pytest. The pure-Python `test` job is unchanged and still skips tests/rust.
- wheels.yml: on a `v*` tag or manual dispatch, builds abi3 wheels with maturin-action (Linux x86_64 manylinux auto, Windows x64) plus an sdist, and attaches them with action-gh-release on tags only.
- radkin and rsoda are public (anonymous `git ls-remote` works), so no token is needed. `CARGO_NET_GIT_FETCH_WITH_CLI` is set everywhere, and also passed into the manylinux container.
- Extra decision: `rust = ["radiopyo-rs"]` is declared by name with a comment (install from a release wheel or `pip install ./rust` first). pip resolves only the extras requested, so `.[dev]` is unaffected. A git direct reference would be refused by PyPI or hatchling.
- gps-config versionFiles now lists pyproject.toml, rust/pyproject.toml and rust/Cargo.toml.
- Ticket 04 follow-up, included here: the full suite exposed a flaky speedup test. Two changes came out of it:
  - Batch workers now buffer their log records and `run_batch` replays them after the batch, so workers no longer take the GIL per record.
  - mimalloc is the global allocator (small gain on Windows).
  The test now uses longer runs (0.5 s span) and best-of-3 timings, and passed 3 times out of 3. Measured: 8 items run 4× faster on 8 threads than on 1. Each solve is about 1.7× slower under parallel load (hybrid cores), so the speedup is below linear.

## Blockers / Challenges

The wheels and CI workflows can only be checked on GitHub (the PR run, plus a `workflow_dispatch` of wheels.yml after the push).

## Commits

- f565911 ci: Rust backend test job and release wheels workflow; rust extra (ticket 05)

## Time Spent

11m (ticket-start.js to ticket-complete.js)
