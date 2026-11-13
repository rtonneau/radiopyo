# Ticket 04: run-batch

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] `_core.run_batch(items, mode, options, n_threads)` runs every item on a rayon pool (`n_threads=None`: all cores) with the GIL released and returns per-item results or errors in input order.
- [ ] `radiopyo.rust.run_batch(cells, t_span, mode="run"|"chunked"|"single-pulse", n_threads=None, return_exceptions=False, **options) -> list[ODEResult | BaseException]`: parity guard per cell before any solve, then raise-after-all (the first error re-raised with its index), or the exceptions in their slots with `return_exceptions=True`.
- [ ] Tests: results equal the single-cell Rust runs (exact), order kept, a failing variant in the middle (raise and `return_exceptions` paths), different configs in one batch, and a 4-cell batch faster than 4 serial Rust runs when at least 2 cores exist.

**Files to Touch:**
- `rust/src/lib.rs`
- `rust/Cargo.toml`
- `radiopyo/rust.py`
- `tests/rust/test_rust_batch.py`

**Verification Step:**

Run:
```bash
cd rust && cargo clippy --release -- -D warnings && <radiobio_env python> -m maturin develop --release && cd .. && <radiobio_env python> -m pytest tests/rust -q && ruff check radiopyo && mypy radiopyo
```

Expected:
clippy clean, all Rust tests pass, ruff and mypy zero errors.

**Notes:**

Convert all Python inputs into owned Rust values before `py.allow_threads`; build the numpy outputs after re-acquiring the GIL. radkin's `OdeResult` has a `Cell<bool>` (not `Sync`) but each result stays on its thread (`Send` is enough). Log forwarding from worker threads: pyo3-log needs the GIL, so make sure it does not deadlock (cache/ignore records off the main thread if needed).
