# Ticket 05: ci-and-release-wheels

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `.github/workflows/tests.yml` gets a job (Ubuntu and Windows) that installs Rust, builds the extension with maturin (release, rust-cache) and runs the full pytest with the Rust tests active, plus `cargo fmt --check` / `cargo clippy -D warnings` in `rust/`; the existing pure-Python job stays and still skips them.
- [ ] `.github/workflows/wheels.yml`: on a `v*` tag (and manual dispatch), maturin-action builds abi3 wheels for Windows x64 and Linux x86_64 (manylinux) plus an sdist and attaches them to the GitHub release.
- [ ] `pyproject.toml` gets the `rust` extra (`radiopyo-rs`, with a comment on installing from the subdirectory or a release wheel); `.work/gps-config.json` release `versionFiles` gains `rust/pyproject.toml` and `rust/Cargo.toml` so the versions stay in lockstep.

**Files to Touch:**
- `.github/workflows/tests.yml`
- `.github/workflows/wheels.yml`
- `pyproject.toml`
- `.work/gps-config.json`

**Verification Step:**

Run:
```bash
<radiobio_env python> -m pip install -e ".[dev]" && cd rust && <radiobio_env python> -m maturin build --release && cd .. && <radiobio_env python> -m pytest -q
```

Expected:
the editable install resolves, a `radiopyo_rs-*-abi3-*.whl` is built in `rust/target/wheels/`, and pytest passes; the PR's CI shows both test jobs green.

**Notes:**

The extra must not break `pip install .[dev]` in the pure-Python job: if `radiopyo-rs` is not resolvable from an index, document the extra as the release-wheel / git-subdirectory install instead of declaring an unresolvable requirement (decide when implementing, and say which). The wheels workflow can only be fully exercised by a tag; validate it with `workflow_dispatch` on the branch if possible.
