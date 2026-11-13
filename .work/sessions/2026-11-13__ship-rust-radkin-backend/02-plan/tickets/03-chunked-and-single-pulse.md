# Ticket 03: chunked-and-single-pulse

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `prepare_chunked_run(...).run(backend=)` and `prepare_single_pulse(...).run(backend=)` (on `SimChunksIter`) dispatch to radkin with the same chunk options (max step sizes on/off/trans/final, trans_time, points per chunk, y0) and return an `ODEResult` with step dose rates, burst dose rates, short ON spans and `thinned` set.
- [ ] Tests: chunked and single-pulse vs Python at rtol 1e-6 on the configs × pulsed / macro-pulsed / single-pulse beams, point thinning, a run cut by `max_dose` producing a short ON span, `reconstruct_reactions` equal to the Python backend's at rtol 1e-6.

**Files to Touch:**
- `radiopyo/simulation/unit_cell.py`
- `radiopyo/rust.py`
- `tests/rust/test_rust_chunked.py`

**Verification Step:**

Run:
```bash
<radiobio_env python> -m pytest tests/rust -q && ruff check radiopyo && mypy radiopyo
```

Expected:
all Rust tests pass; ruff and mypy zero errors.

**Notes:**

`SimChunksIter` already holds the chunk options; forward them rather than rebuilding chunks in Python, so radkin cuts the chunks itself (it reproduces `prepare_chunked_run`).
