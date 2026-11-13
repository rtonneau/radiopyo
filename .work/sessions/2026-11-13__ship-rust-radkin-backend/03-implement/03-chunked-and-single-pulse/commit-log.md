# Ticket 03: chunked-and-single-pulse

**Status:** ✅ Done

## Local Test Result

`python -m pytest tests/rust`: 34 passed (10 new in test_rust_chunked.py). Full suite: 844 passed. `ruff check radiopyo tests/rust` and `mypy radiopyo`: clean.

## Review Notes

Base 6bddbbb. Checked against the criteria:
- `SimChunksIter.run(backend=)` exists. prepare_chunked_run and prepare_single_pulse attach a `RustRequest` (mode, t_span, chunk options) that is forwarded to radkin's chunked/single_pulse modes.
- Pulsed, macro-pulsed and single-pulse runs match the Python backend on the final state (rtol 1e-6).
- Thinned runs match on the full grid.
- The short ON span from a max_dose cut matches.
- reconstruct_reactions matches on a thinned run.
- A hand-built iterator, or one already partly run, raises RustBackendError.

One change after the first run: full-grid comparison was dropped for unthinned runs, because the two LSODA codes choose their own steps (4007 vs 4012 points).

## Blockers / Challenges

The step grids differ between the backends, so for unthinned runs the tests compare only the final state.

## Commits

- ec5922c feat: backend="rust" for chunked and single-pulse runs (ticket 03)

## Time Spent

8m (ticket-start.js to ticket-complete.js)
