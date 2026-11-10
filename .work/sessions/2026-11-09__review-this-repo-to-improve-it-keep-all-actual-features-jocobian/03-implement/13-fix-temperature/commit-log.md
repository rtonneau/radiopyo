# Ticket 13: fix-temperature

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`: 471 passed (460 baseline + 11 new in tests/simulation/test_temperature.py), ruff "All checks passed!", mypy "Success: no issues found in 42 source files". Goldens (tests/golden) untouched and passing at the default temperature. Before/after at 25 and 37 degC: `.scratch/tests/2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian/evidence/13-before-after.md` (25 degC identical, 37 degC H2O +0.78%).

## Review Notes

- Parent: base ed6d27c4549dbefc5731d1cc66da78eab0e46839; pre-existing changes: gps session state only.

- One source of truth: `env.parameters["temperature"]` (degC); `UnitCell.T` is a view with a setter, `_T` kept as a read-only property (used by existing tests), `set_temperature` delegates to `SimEnv.set_temperature` then rebuilds the matrices. Config `temperature` key now honoured.
- Fixed-species vector, `initial_values()`, `ODEResult.to_pandas` and `__getitem__` use the simulation temperature (ODEResult deep-copies the env, so it keeps its own snapshot). `build_all_matrices` / `make_kreaction_matrix` / `make_cst_species_vector` / `initial_values` default to the env temperature; explicit `T` still overrides. Docstrings state units (k: 1/K internal, pKa and initial cc: degC).
- Before/after report written. Left alone as out of scope: `ReactionsList.iter_*_kvalues(T=None)` still defaults to 25 degC when called directly; `convert_to_config_dict` writes constant-polynomial species at default T (temperature independent).

- Parent acceptance: USER APPROVED on 2026-11-10 (evidence/13-before-after.md): 25 degC bit-identical, no golden change. Reviewer skipped: parent reviewed diff directly (localized to temperature plumbing, goldens unchanged); 471 passed, ruff 0, mypy 0. Seed for ticket 18 (user asked to log, not fix): H2O cc_init polynomial has all-positive coefficients, 55.93 M at 25 degC and rising with T (real 55.34 M, falling).

## Blockers / Challenges

None. Note: bundled configs have no temperature-dependent k, so the k-value test builds a config with an extra k polynomial.

## Commits

- bafa16f fix: one temperature source of truth, fixed species at simulation temperature (ticket 13)

## Time Spent

10m (ticket-start.js to ticket-complete.js)
