# Ticket 04: run-batch

**Status:** ✅ Done

## Local Test Result

`cargo clippy --release -- -D warnings`: clean. `maturin develop --release`: installed. `pytest tests/rust`: 42 passed (8 new in test_rust_batch.py). `ruff check radiopyo tests/rust`: clean. `mypy radiopyo`: no issues.

## Review Notes

Base ec5922c (after ticket 03's session commit d9edc8d).
- `_core.run_batch(items, n_threads)`: every item is `(cell, mode, t_span, options)`, which is more general than one mode and option set per batch. Options are parsed into owned `Job`s before `py.detach`, solves run on a rayon pool, and dicts or `RadkinRustError` instances are built after the GIL is back, in input order.
- `radiopyo.rust.run_batch(cells, t_span, mode, n_threads, return_exceptions, **options)` accepts 'single-pulse' or 'single_pulse'. It checks unknown options before any solve, and runs the beam check and parity guard per cell before solving (a failing cell is not sent to radkin). It raises after all cells finish: the first error is raised with an `add_note` naming its index. With `return_exceptions=True`, the errors stay in their slots.
- Tests cover: equality with single Rust runs (exact, chunked and single-pulse), plain run mode with t_eval, mixed configs, order, a failing cell (no beam, parity) on both paths, bad inputs, and a speedup over serial runs.

## Blockers / Challenges

Worker threads' log records go through `Python::attach`. The main thread has released the GIL, so there is no deadlock; logging is only serialized.

## Commits

- 3470614 feat: radiopyo.rust.run_batch, parallel radkin runs on a rayon pool (ticket 04)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
