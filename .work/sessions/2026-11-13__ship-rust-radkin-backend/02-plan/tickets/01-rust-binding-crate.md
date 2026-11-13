# Ticket 01: rust-binding-crate

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] `rust/Cargo.toml` (crate `radiopyo_rs`, cdylib, PyO3 abi3-py311, rust-numpy, pyo3-log, radkin git dep pinned by `rev`), `rust/pyproject.toml` (maturin backend, package `radiopyo-rs`, version = radiopyo's), `rust/src/lib.rs` exposing module `radiopyo_rs._core`.
- [ ] `_core.run(config_toml: str, state: dict, mode: str, options: dict) -> dict`: builds a radkin `UnitCell` from the merged config text, applies temperature / beam (type + `as_dict()` table via `Beam::from_table`) / LET, runs `run` | `run_till_max_dose` | `chunked` | `single_pulse` and returns t (n,), y (n_dyn, n), burst dose rates, step dose rates, short ON spans, thinned, nfev/njev/nst, exec time, species labels (dyn + fixed), reaction labels, k values and G yields.
- [ ] Errors carry a kind (solver failed, no beam, short beam ON, config, value) for the Python mapping; radkin `log` records go to Python's `radiopyo` logger.
- [ ] `maturin develop --release` in `radiobio_env` builds and imports; `cargo fmt --check` and `cargo clippy -- -D warnings` clean in `rust/`.
- [ ] If radkin needs an addition, it lands as a merged PR in `rtonneau/radkin` and the pinned rev includes it.

**Files to Touch:**
- `rust/Cargo.toml`
- `rust/pyproject.toml`
- `rust/src/lib.rs`
- `rust/README.md`
- `.gitignore`

**Verification Step:**

Run:
```bash
cd rust && cargo fmt --check && cargo clippy --release -- -D warnings && <radiobio_env python> -m maturin develop --release && <radiobio_env python> -I -c "import radiopyo_rs._core as c; print(c.__version__)"
```

Expected:
clippy clean, the extension builds and prints the radiopyo version.

**Notes:**

Python → Rust config: serialize the merged `ConfigDict` (without `includes`, already resolved) to TOML text in Python (a tiny writer, or JSON then `serde_json` → `toml::Value` in Rust if simpler; no new Python dependency) and call `Config::from_toml_str(text, None)`. Check radkin's `RadkinError` variants for the error kinds. Build with `--release` (debug LSODA is slow). Keep radkin's AGENTS.md rule: no committed `path` dependency (`[patch]` only locally).
