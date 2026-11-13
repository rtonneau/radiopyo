# Ticket 01: fenton-map-rust-backend

**Status:** ✅ Done

## Local Test Result

Ran the Verification Step (in `docs/issues/fenton`, radiobio_env); output in `.scratch/tests/<session>/01-verify.txt`:

```
python (cache) load 0.0031 s        # existing .map_cache/fenton_map_4ace6041a0a8.npz still found
rust 12.0 s                         # compute_map(backend="rust", cache=None), 144 points
mask equal True True nan count 0
max rel diff 3.76e-09 9.96e-09      # rooh_max_dose, rooh_end vs the cached Python map
ValueError: backend must be one of ('python', 'rust'), got 'julia'
```

Python path after the refactor: `compute_point` at grid (5, 7) gives (0.21064797446, 2.63716624873) vs cached (0.21064797444, 2.63716624852).

## Review Notes

Base: ed191db7b32909e4fbe5b5aa76f694a3beae763b. Pre-existing changes, owned by the user and not touched: `docs/issues/fenton/Fenton_Fixed_Showcase.ipynb`, `examples/notebook_basic_usage.ipynb`, `examples/notebook_rust_backend.ipynb` (kernel metadata only).

Inline mode, parent review of the diff against the Acceptance Criteria:
- `compute_map(backend="rust")` returns the same keys and shapes and uses `rp.rust.run_batch(mode="chunked")` with `MAX_STEP_ON` / `MAX_STEP_OFF`.
- The shared helpers `_make_beam` (feasibility, with the original `1 / frequency` arithmetic) and `_rooh_pair` (ROOH at max dose and at the end) are used by `compute_point` and `_map_rust`. `return_exceptions=True`, so a failed point is logged and left NaN.
- The cache key adds `backend` only when it is not "python": the existing Python cache file loaded instantly.
- An unknown backend raises `ValueError` before any work.
- The Rust map vs the Python map: same NaN mask, max relative difference 1e-8.
- With the Rust backend, `max_workers` means the rayon thread count.
- No independent reviewer: a localized change in an issue-folder helper, checked directly against the cached reference.

## Blockers / Challenges

The first draft of `_make_beam` computed the pulse period with a wrong expression. I fixed it to the original `1.0 / (avg / dpp)` before verification. Ruff on this file (outside the `radiopyo/` lint gate) reports only issues that were already there (`__all__` order, RUF100 on `noqa: BLE001`, S324 `sha1`), which I left as they were.

## Commits

- 977e8f9 docs(fenton): Rust backend for fenton_map.compute_map via run_batch (ticket 01)

## Time Spent

2m (ticket-start.js to ticket-complete.js)
