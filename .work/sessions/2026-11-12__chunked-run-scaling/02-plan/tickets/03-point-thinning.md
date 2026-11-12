# Ticket 03: point-thinning

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] `prepare_chunked_run` and `prepare_single_pulse` accept `points_per_chunk_on` / `points_per_chunk_off` (int >= 2 or None; ValueError otherwise); the OFF value covers OFF and transition chunks, the final chunk follows its beam state.
- [ ] A thinned chunk stores exactly N points: chunk start (after the initial-burst shift for ON chunks) and N-1 log-spaced times up to the exact chunk stop; solver steps unchanged.
- [ ] Defaults None leave results unchanged (golden tests pass).
- [ ] The next chunk's y0 equals the chunk end of the full-resolution run to solver tolerance (test).
- [ ] `ODEResult` records in memory that it is thinned; `integrate_production_of`, `integrate_consumption_of`, `detail_*`, `species_balance_evolution`, `integrate_species`, `integrate_single_reaction` log one warning (once per result) on a thinned result.
- [ ] A short ON span emits no t_eval warning for generated times.
- [ ] New test `tests/simulation/test_point_thinning.py`; full pytest, ruff and mypy pass.

**Files to Touch:**
- `radiopyo/simulation/unit_cell.py`
- `radiopyo/simulation/ode_result.py`
- `tests/simulation/test_point_thinning.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All tests pass, ruff and mypy clean.

**Notes:**

Carry the point count on `SimChunk` (new optional field, default None) and build the chunk `t_eval` in `SimChunksIter.__next__` / `_raw_run` after `_apply_initial_burst` shifts t_span (replace `_drop_skipped_t_eval` use for generated times). Log-space lower bound: span·1e-6.
