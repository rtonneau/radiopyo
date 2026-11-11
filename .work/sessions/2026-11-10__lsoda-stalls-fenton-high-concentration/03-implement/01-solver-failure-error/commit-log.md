# Ticket 01: solver-failure-error

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`
- pytest: 631 passed in 90.93s
- ruff: All checks passed!
- mypy: Success: no issues found in 42 source files
Output saved in `.scratch/tests/2026-11-10__lsoda-stalls-fenton-high-concentration/evidence/ticket01-verify.txt`.

## Review Notes

- `SolverFailedError(SimulationError)` (so a `RadioPyoError`) added to `radiopyo/simulation/exceptions.py`, in `__all__`.
- `UnitCell._raw_run` raises it when `sol.success` is False; message: method, t_span, last reached t, scipy message. `run` and chunked/single-pulse runs go through `_raw_run`.
- `tests/simulation/test_solver_failure.py`: hierarchy, `run` raises (message content checked), chunked run raises, successful run unaffected. Failure forced by monkeypatching `unit_cell.solve_ivp` (real solve, then success=False), deterministic and fast.
- CHANGELOG Unreleased/Added bullet.
- Review rounds: 0 (self-check only).

## Blockers / Challenges

Two existing tests in `tests/simulation/test_run_options.py` faked the solver with the string `"sol"`, which has no `.success`; they now return a `types.SimpleNamespace(success=True)` sentinel (same assertions otherwise).

## Commits

- 5804aa1 feat: raise SolverFailedError when solve_ivp fails (ticket 01)

## Time Spent

6m (ticket-start.js to ticket-complete.js)
