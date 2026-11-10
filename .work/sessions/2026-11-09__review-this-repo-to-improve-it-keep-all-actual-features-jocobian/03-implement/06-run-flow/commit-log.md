# Ticket 06: run-flow

**Status:** ✅ Done

## Local Test Result

Command (repo root, Anaconda Python on PATH):

```
$ python -m pytest -q && ruff check radiopyo && mypy radiopyo
316 passed in 25.28s
All checks passed!
Success: no issues found in 36 source files
exit=0
```

That is 287 baseline tests plus 29 new ones in `tests/simulation/test_run_options.py`. Evidence is in `.scratch/tests/2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian/evidence/`:

- `06-verify.txt`: the run above.
- `06-golden-seeds.txt`: `tests/golden` + `test_run_options.py` pass under PYTHONHASHSEED 0, 3 and 13 (34 passed each time).
- `06-run-flow-equivalence.txt`: compares base 33dc22b (`git archive` into `base06/`) with the working tree on 22 run-flow scenarios, using `run_flow_dump.py` and `run_flow_compare.py` in the scratch dir. The scenarios are:
  - `run`: defaults, `force_positive=True`, custom tolerances and max_step, RK45, `t_eval`, empty span with and without the burst, `use_beam=False`, user `y0` (mutated or not)
  - `FIRST_STEP` clamped by max_step and by the span; too-short pulse; `use_jac=True`
  - `run_till_max_dose` with and without kwargs
  - `prepare_chunked_run`: default, with final max step + `force_positive` + `y0`, and `use_beam=False`
  - `prepare_single_pulse`: default and with every kwarg
  - all three entry points without a beam

  Every array is `bit_identical=True` under PYTHONHASHSEED 0 and 7 (`ALL IDENTICAL`, exit 0). This covers `t`, `y`, `n_points`, chunk tables, `y0` in-place deltas and error messages. Species order was fixed by label through `[parameters.species_order]`, because the default order is not reproducible from one run to the next even with a fixed hash seed. A base-vs-base rerun is also `ALL IDENTICAL`, which validates the comparison method.

## Review Notes

- Parent: base 33dc22bdf16a92a8ff5cebc6edfa773c53437829; pre-existing changes: gps session state only.

Acceptance Criteria check:

- **One `RunOptions`.** It is a frozen, keyword-only dataclass in `unit_cell.py` (exported in `__all__`) with fields atol, rtol, method, max_step, first_step, use_beam, force_positive, use_jac and t_eval.
  - `run`, `prepare_chunked_run` and `prepare_single_pulse` build one from their keyword arguments.
  - `SimChunksIter(sim, chunks, options=None, y0=None)` stores it, and each chunk uses `dataclasses.replace(options, max_step=chunk.max_step_size)`.
  - `first_step=None` means "read `self.FIRST_STEP` at run time", so instance overrides still work. `atol`/`rtol=None` resolve to the new module constants `DEFAULT_ATOL`/`DEFAULT_RTOL`.
  - Public signatures (names, order, defaults) of `run`, `prepare_chunked_run`, `prepare_single_pulse` and `run_till_max_dose` are unchanged, and `test_public_signatures_unchanged` pins them.
- **`_raw_run` split into named steps.** `_raw_run(t_span, options, y0=None, beam_on=True)` now runs these steps in order:
  1. `use_jac` check (still `NotImplementedError`)
  2. no-beam check
  3. `widen_empty_span`: the old `t0/100` trick, now a module function with a docstring that explains the `[t, t] -> [t/100, t]` reason and the `t=0` limit
  4. `_solver_kwargs`: first_step/max_step/span clamps, done before the burst shift as at base
  5. `_apply_initial_burst`: T_ZERO shift, the "Too Short" check, and the in-place `y0` increment
  6. `_timed_solve`: `solve_ivp` call and `last_exec_time`
- **Shared prepare code.** `_require_beam()` does the validation (log + raise `NoBeamDefinedError`). `_make_chunks(boundaries, beam_on, max_steps)` builds the `SimChunk` lists for both methods. `_chunks_iter(...)` does the logging (same messages, same order) and builds the `SimChunksIter`.
  - `prepare_chunked_run` keeps its final-chunk rule: `max_step_size_final`, otherwise on/off. `intervals[-2]` is still indexed, so the degenerate case of fewer than 2 rows still fails with `IndexError`.
  - `prepare_single_pulse` keeps "no fallback" for the final chunk. `trans_time` defaults to `DEFAULT_TRANS_TIME = 100e-6`.
- **Bug reproduced explicitly.** At base, `run()` called `_raw_run` positionally with 10 arguments, while `_raw_run` has 11 parameters, with `beam_ON` before `force_positive`. So run's `force_positive` landed in `beam_ON`, and `_raw_run`'s `force_positive` kept its default `False`. `run()` now builds `RunOptions(force_positive=False, ...)` and calls `_raw_run(..., beam_on=force_positive)`, with a `# TODO(ticket 12)` comment. The run docstring describes the effect: the default skips the burst, `force_positive=True` applies the burst and never forces positivity. This is pinned by `test_run_reproduces_argument_order_bug[False|True]` and `test_run_default_skips_burst`.
- **Docstrings.** They now state the real defaults: atol 1e-12, rtol 1e-10, FIRST_STEP 1e-15 (the class docstring said 1e-12), max step inf, method "LSODA", trans_time 100e-6, `time_start` 1e-9. They also give the real return types: `run` returns `ODEResult` (was documented as `ScipyODEResult`), and `set_beam` returns `Self`.
  - `prepare_single_pulse` now documents `max_step_size_trans` and `trans_time`. It used to list a non-existent `max_step_size_off`.
  - `set_beam` documents that `LET=None` silently reuses the previous beam's LET, as the ticket note asks.
  - `run_till_max_dose`, `SimChunk.beam_on` and the `SimChunksIter` attributes are documented. The "chucks" typo is fixed.
- **`KwargsRun`** (`sim_types.py`) is now `total=False` and gains `use_beam`/`force_positive`. Before, `Unpack` of a total TypedDict made every key mandatory for type checkers. This change is typing only, with no runtime change: a `t_span` passed to `run_till_max_dose` is still silently overwritten (documented).
- **Golden tests** pass unchanged; nothing in `tests/golden/` was touched.
- **New tests.** `tests/simulation/test_run_options.py` has 29 tests. They cover RunOptions defaults and frozen/kw-only behaviour, constants, `widen_empty_span`, the `_solver_kwargs` clamps, the in-place burst and the too-short error, the `_raw_run` step wiring, `run` and `run_till_max_dose` forwarding, the argument-order bug, public signatures, the chunk tables of both prepare methods (including the final-chunk rules), `NoBeamDefinedError`, and `SimChunksIter` per-chunk calls (options, `beam_on`, previous last point as `y0`). The tests reuse a parsed cell through module-scoped templates plus `copy.deepcopy`, so the file runs in about 6 s instead of about 34 s.

Left as is on purpose:

- **Missing `raise`.** `_raw_run` still creates `rp_excep.NoBeamDefinedError("Simulation aborted")` without raising it, with a NOTE pointing to ticket 12. A beamless `run()` still ends in `KeyError: ''`, identical to base.
- **`y0` mutation.** It is kept: the burst increments the caller's `y0`, or in chunked runs the view `self.sol.y[:, -1]` of the previous chunk's last stored point. This is documented in the docstrings.
- **Burst with `use_beam=False`.** It still uses `self.beam` even when `use_beam=False` (ticket 12).
- **Temperature** is untouched (ticket 13).

API notes for the parent:

- `SimChunksIter.__init__` changed from `(sim, chunks, atol, rtol, use_jac, method, y0, use_beam, force_positive)` to `(sim, chunks, options=None, y0=None)`. Its old public attributes `atol`, `rtol`, `use_jac`, `method`, `use_beam` and `force_positive` are now under `.options`. `SimChunksIter` is not in `__all__` and is only built by `prepare_*`; grep finds no other user in `radiopyo/` or `tests/`.
- `_raw_run` (private) has a new signature.
- New public module names in `unit_cell.py`: `RunOptions` (in `__all__`), plus `DEFAULT_ATOL`, `DEFAULT_RTOL`, `DEFAULT_TRANS_TIME` and `widen_empty_span`.

- Parent acceptance: independent reviewer (opus/high): no findings; parent re-ran pytest (316 passed), ruff 0, mypy 0; base-vs-new bit identical on 22 run-flow scenarios (evidence/06-run-flow-equivalence.txt). run() burst bug reproduced explicitly with TODO(ticket 12). SimChunksIter ctor changed (internal, not exported).

## Blockers / Challenges

- My first base-vs-new equivalence attempt failed even when comparing base with itself under the same PYTHONHASHSEED. The default species order comes from set iteration over objects and changes from one process to the next. I fixed it by building the cells in the scenario script with a label-sorted `[parameters.species_order]`. After that, base vs base and base vs new are both bit-identical.
- `ruff format --check radiopyo/utils/sim_types.py` reports a missing blank line after the module docstring. This is pre-existing (HEAD has the same result), not in the Verification Step, and I left it alone. `ruff check` is clean.

## Commits

- d24bc92 refactor: RunOptions and stepwise _raw_run in UnitCell (ticket 06)

## Time Spent

21m (ticket-start.js to ticket-complete.js)
