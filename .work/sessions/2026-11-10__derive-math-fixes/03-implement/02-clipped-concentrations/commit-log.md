# Ticket 02: clipped-concentrations

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`

- pytest: 597 passed (92.6 s); ruff: All checks passed!; mypy: Success: no issues found in 42 source files.
- The new/extended tests (test_clipped_concentrations.py, test_jacobian.py, test_reaction_rates.py) were also run against the pre-change `ode_functions_solver.py`: 7 of them fail there (self-reacting species consumed below 0, reverse reaction run backwards, Jacobian column non-zero for c < 0, rates not equal to clipped-state rates, MM with negative enzyme/substrate), all pass with the fix (evidence: .scratch/tests/2026-11-10__derive-math-fixes/evidence/02-tests-on-old-code.txt).

## Review Notes

Parent: base d843b340a94b270952bb75c710cbf0dfdde332f7; pre-existing changes: only gps session state (.session-config.json, 03-implement/), owned by the parent.

AC1 (clipped rate laws): `mass_action_rates`, `michaelis_menten_rates`, `mass_action_derivatives`, `michaelis_menten_derivatives` read max(c, 0) (enzyme and substrate included); `reaction_rates` unchanged in structure, still shared by `derive` and `ODEResult.reconstruct_reactions`. The negative-product block (c <= DEPLETION_LIMIT) is unchanged and agrees with clipping (a clipped 0 is depleted).
AC2 (Jacobian): a species with c < 0 has its column zeroed in both derivative functions (k and Michaelis-Menten, each of the enzyme and substrate); at exactly c = 0 the c > 0 side is used (tests: A = 0 gives d rate/dA = k B; S = 0 gives k E / Km).
AC3: test_clipped_concentrations.py (2 A -> B at A < 0 gives dA/dt = dB/dt = 0; LSODA run with the analytic Jacobian from A0 = -1e-3 uM under irradiation stays >= A0 and ends > 0; reverse reaction not run backwards).
AC4: derive == S.r for the full stacked rate vector, and the conserved combinations A+C+2E, B+C, S+P have zero derivative, on positive / negative / near-depleted / mixed states; plus rate-equals-clipped-state tests (scalar and (n_species, n_t)) in test_reaction_rates.py and Jacobian tests in test_jacobian.py (also vs finite differences at negative states).
AC5 (golden): ran tests/golden unchanged after the code change. Deviation vs recorded data (max abs / max rel), before the fix (HEAD solver) -> after the fix, PYTHONHASHSEED=0 (scripts and output in .scratch/tests/2026-11-10__derive-math-fixes/, evidence/02-golden-deviation-before.txt and 02-golden-deviation.txt):
- derive: labarbe 2.384e-07 / 7.774e-16 -> identical; pure_water 0 / 0 -> identical. Test passes: not re-recorded.
- beams: 0 / 0 -> 0 / 0. Test passes: not re-recorded.
- runs: labarbe_run species 1.25e-10 / 6.3e-09 -> 9.6e-11 (rel meaningless, golden ~0); labarbe_curve species 3.7e-10 -> 1.1e-10; labarbe_curve rate 8.75e-10 -> 8.75e-10; chunked_run species 1.314e-10 / 2.56e-07 -> identical; chunked_run rate 3.71e-12 / 2.56e-07 -> identical. These deviations exist before the fix too (solver/step noise vs the old recording, within the RTOL_RUN / ATOL tolerances); tests pass: not re-recorded.
- hdf: test_hdf5_written_by_old_code_is_readable FAILS: only the 'rebuilt_rate' expectation of hdf_labarbe (rates recomputed by the current code on the stored time series). Before 0 / 0, after max abs 1.603e-11, max rel 1.0 (191 of 3720 values). Cause: the old stored series holds H_r = -1.60e-13 uM (a real overshoot), the old code gave rate(RH + H_r -> RH_r) = -1.60e-11, the clipped law gives 0. The stored species and stored rates (read from the .h5) still match.
Re-record: the `--only hdf` recorder would rewrite golden_results.h5 and its json store, which is forbidden. So I did NOT run record_golden.py for any part; instead I updated only the key `hdf_labarbe__rebuilt_rate` of tests/golden/data/hdf_expected.npz with the current code's output on the unchanged .h5 (script .scratch/tests/2026-11-10__derive-math-fixes/02-patch_hdf_expected.py). golden_results.h5 sha256 verified unchanged, the .golden_results/ json store untouched. derive.json untouched, so the stale `*_beam_on_force_positive` entries remain and `test_derive_matches_golden` keeps its subset check (restoring the exact-set check only applies if derive is re-recorded).
AC6: pytest, ruff, mypy at zero (see Local Test Result).
- Parent review (inline; no separate reviewer: localized to ode_functions_solver.py, whole diff inspected, covered by FD Jacobian tests): clipping in mass_action_rates / michaelis_menten_rates and the c<0 zero columns in both derivative functions are correct (the chained out[:, :n][:, mask] assignment writes through a view). Golden: only hdf_labarbe__rebuilt_rate changed (stored H_r = -1.6e-13 µM gave a negative rate before, 0 now), patched in hdf_expected.npz; golden_results.h5 untouched (sha256 cb2d1bf…). Accepted; implementer reported 597 passed, ruff and mypy clean.

## Blockers / Challenges

The only golden failure is the hdf `rebuilt_rate` expectation, which cannot be re-recorded with `record_golden.py --only hdf` without rewriting golden_results.h5; handled by patching hdf_expected.npz alone (see Review Notes). The parent may prefer a different handling (e.g. a tolerance) and should confirm. Initial attempt to append tests with a shell heredoc failed on quoting before running anything; used the file tools instead.

## Commits

- ba8d133 fix: rate laws read clipped concentrations max(c, 0) (ticket 02)

## Time Spent

8m (ticket-start.js to ticket-complete.js)
