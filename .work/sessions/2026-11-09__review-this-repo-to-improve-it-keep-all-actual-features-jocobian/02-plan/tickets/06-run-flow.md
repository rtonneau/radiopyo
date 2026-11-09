# Ticket 06: run-flow

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] One `RunOptions` dataclass (atol, rtol, method, max_step, first_step, use_beam, force_positive, use_jac, t_eval) used by `run`, `prepare_chunked_run`, `prepare_single_pulse`, `SimChunksIter`; the public keyword arguments and their defaults are unchanged.
- [ ] `_raw_run` is split into small named steps (solver kwargs, initial-burst application, solve, timing); the `t0/100` handling is named and documented.
- [ ] `prepare_chunked_run`/`prepare_single_pulse` share validation and chunk-building code.
- [ ] Today's `run()` behaviour (burst skipped, as caused by the argument-order bug) is reproduced explicitly with a comment pointing to ticket 12.
- [ ] Docstrings state the real defaults (atol 1e-12, FIRST_STEP 1e-15, ...) and real return types.
- [ ] Golden tests pass.

**Files to Touch:**
- `radiopyo/simulation/unit_cell.py`
- `radiopyo/utils/sim_types.py` (`KwargsRun`)
- `tests/simulation/test_run_options.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass, zero lint/type errors.

**Notes:**

Do not fix the burst or `y0` mutation here; ticket 12 does. `set_beam(beam)` without LET silently reusing the previous LET is current behaviour: keep and document.
