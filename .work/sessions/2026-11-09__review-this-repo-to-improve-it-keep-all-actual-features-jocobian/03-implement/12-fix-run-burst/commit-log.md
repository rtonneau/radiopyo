# Ticket 12: fix-run-burst

**Status:** ✅ Done

## Local Test Result

Phase 1 (goldens NOT regenerated, awaiting user approval).

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`:
- pytest: 458 passed, 2 failed (460 collected vs 450 at baseline: +10 net new tests in tests/simulation/test_run_options.py). The 2 failures are the expected golden tests, caused only by the burst / y0 fix: `tests/golden/test_golden.py::test_labarbe_run_matches_golden` (run() now starts at 1e-7 + T_ZERO = 2e-7 s instead of 1e-7 s, first point carries the burst, t_eval times before 2e-7 s dropped) and `::test_chunked_pulsed_run_matches_golden` (the burst no longer overwrites the last stored point of the previous OFF chunk: H2, H2O2, H3O_plus, H_r, OH_r, e_aq differ at the end-of-OFF-chunk points only). The derive, beams and HDF golden tests still pass.
- ruff check radiopyo: All checks passed
- mypy radiopyo: Success, no issues found in 42 source files
- Before/after report: .scratch/tests/2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian/evidence/12-before-after.md (raw: metrics_before.json, metrics_after.json, 12-pytest-after.txt)

## Review Notes

- Parent: base 8c8964990989bbd72ff06e5af4165b825bf2d0fc; pre-existing changes: gps session state only.

- AC1: run() passes force_positive to RunOptions (reaches derive via solve_ivp args) and always calls _raw_run(beam_on=True); tests test_force_positive_reaches_derive, test_run_applies_burst_by_default.
- AC2: _raw_run works on a copy of y0 (user y0 and the previous chunk's stored last point untouched); tests test_run_never_mutates_y0, test_sim_chunks_iter_keeps_previous_last_point (fails when the copy is removed), test_sim_chunks_iter_never_mutates_y0.
- AC3: burst only when beam_on and use_beam. ODEResult records the dose rate of every applied burst (`burst_dose_rates`, set via `append_sol(..., bursts=)`, merged in `__add__`); integrate_production_of sums them. Results read from files (unknown, None) keep the former n_pulse * peak assumption. Note: reconstruct_reactions still uses the result's beam (not use_beam) for the radiolytic rate; left as is.
- AC4: _raw_run calls _require_beam() (raises NoBeamDefinedError; also with use_beam=False).
- AC5: report written; goldens untouched pending approval.
- Beyond the spec: with the burst on, a user t_eval starting inside the skipped T_ZERO interval made solve_ivp raise ("Values in t_eval are not within t_span"); such times are now dropped with a warning (UnitCell._drop_skipped_t_eval). Added UnitCell.last_burst_dose_rate. Replaced the tests pinning the old behaviour (argument-order test, default-skips-burst test, in-place/`is y0` assertions). CHANGELOG `### Fixed` entry added.

- Parent acceptance: USER APPROVED before/after (evidence/12-before-after.md) on 2026-11-10; independent reviewer (sonnet/high): 3 P3 — short-span ValueError documented + guard now >= (fixed by parent), burst_dose_rates not persisted (documented in CHANGELOG, seed ticket 18), chunked integrals noisy (noted to user). Goldens: only runs.npz re-recorded (record_golden.py --only runs, PYTHONHASHSEED=1); HDF5 written by original code kept. Parent: 460 passed, golden 5/5 at seeds 0/7/42, ruff 0, mypy 0.

## Blockers / Challenges

- The tool's heredocs choked on triple quotes: edit scripts were written under .scratch/ instead (edit_*.py, untracked/ignored).
- Decision for the parent/user: dropping t_eval times in the skipped interval (vs. raising) is a behaviour choice; also the Labarbe golden `labarbe_curve` will lose its first t_eval point(s) and its time axis starts at 2e-7 s once regenerated.

## Commits

- af5b8ca fix: apply initial burst in run(), never mutate y0, raise NoBeamDefinedError (ticket 12)

## Time Spent

6h 25m (ticket-start.js to ticket-complete.js)
