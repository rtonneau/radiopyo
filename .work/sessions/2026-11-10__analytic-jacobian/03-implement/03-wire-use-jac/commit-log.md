# Ticket 03: wire-use-jac

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q` -> 585 passed (82 s); `ruff check radiopyo` -> All checks passed; `mypy radiopyo` -> no issues in 42 source files.

## Review Notes

AC1: NotImplementedError removed from `_raw_run`; `_timed_solve` passes `jac=jacobian` (same args as derive via `args`) only when `use_jac` and method in `JAC_METHODS` (LSODA, BDF, Radau). AC2: test_use_jac_rk45_no_warning runs RK45 with use_jac=True under warnings-as-errors. AC3: new tests/simulation/test_use_jac.py (LSODA/BDF/Radau, chunked, default False); old NotImplementedError test removed from test_reaction_rates.py; docstrings updated. Review rounds: 0. Gap: tests do not assert that jac is actually received by solve_ivp (only result agreement).

Parent: base 35152f95d390a0b5b12aa899068047ff77071ee7, no pre-existing changes. Independent review skipped: localized change in one method (_timed_solve) plus docstrings; parent read the full diff. Implementer's gap closed by the parent: solve_ivp patched with a wrapper on LABARBE config shows jac=jacobian passed for LSODA and BDF and absent for RK45. Accepted (implementer: sonnet, medium).

## Blockers / Challenges

None

## Commits

- d37c2f8 feat: wire use_jac to solve_ivp for LSODA, BDF and Radau (ticket 03)

## Time Spent

7m (ticket-start.js to ticket-complete.js)
