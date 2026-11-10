# Ticket 04: catalase-k

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q -p no:cacheprovider && ruff check radiopyo && mypy radiopyo`

- pytest: 619 passed (85 s); ruff: All checks passed!; mypy: Success: no issues found in 42 source files (evidence/04-verification.txt).
- golden_results.h5 sha256 after: cb2d1bfbdaf73d0ab10985e2161226518a20ac52dc7d57fa01736643b580cc9e (unchanged).

## Review Notes

AC1: reactions_biology.toml catalase k_value = 6.62e7, comment says per-reaction k, used as given (ADR 0006), assumed Labarbe's Rb10 value and not verified.
AC2 (golden, PYTHONHASHSEED=0, script 02-golden_dev.py; evidence/04-golden-deviation-before.txt / -after.txt; deviation vs the old recording, max abs / max rel):
- derive/labarbe: 2.384e-07 / 7.8e-16 -> 2.554e-03 / 1.44e-08; derive/pure_water 0 -> 0 (pure water config does not use the biology file).
- beams: 0 -> 0.
- labarbe_run species 9.6e-11 -> 3.28e-05; rate 3.7e-13 -> 3.36e-07; labarbe_curve species 1.07e-10 -> 1.57e-01; rate 8.75e-10 -> 8.30e-06.
- chunked_run: deviation vs the old record unchanged (does not use catalase); but recording the `runs` part rewrote it too: this is the first runs.npz re-record of the session, so it folds in the within-tolerance drift from tickets 02/03 (chunked_run species max_abs 1.3e-10, rate 3.7e-12, n_points changed; labarbe rates gained exact zeros from clipping, hence the smaller compressed npz; same 20 keys and shapes). hdf rebuilt_rate: 0 -> 0 (unchanged, test passes).
- Failing tests after the change: test_derive_matches_golden, test_labarbe_run_matches_golden only. Re-recorded with `PYTHONHASHSEED=0 python tests/golden/record_golden.py --force --only derive runs` (derive.json, runs.npz modified; the recorder rewrote the other files with identical content, git shows no change for them, h5 sha256 verified). hdf was not re-recorded; the pinned doubled-k json-store test passes unchanged.
- Because derive was re-recorded, the stale *_beam_on_force_positive entries are gone; test_derive_matches_golden restored to the exact-set check (== instead of <=).
AC3: grep of 1.324e8 / 6.62e7: remaining hits are tests/parser/test_includes.py (lines 311, 338, 342), which write their own base toml with their own catalase reaction, so they do not mean the bundled value and were left unchanged; tests/golden/test_golden.py:486 is the pinned doubled-k test (kept); CHANGELOG line 90 documents the old stored-file doubling (kept). README already says 6.62e7. CHANGELOG has no entry for this fix; not in the ticket scope.
AC4: see Local Test Result.

- Parent: base e54ee0f74f4bc6c0e0d8a2360f1709897a09b376; pre-existing changes: only gps session state, owned by the parent.
- Parent review round 1: fresh reviewer REVIEWED, no correctness defects. It checked: runs.npz keeps the same keys and shapes (the size drop comes from the exact zeros clipping produces); derive.json only lost the force_positive entries, pure_water entries are bit-identical, the labarbe dH2O2 and dO2 changes match a halved catalase k exactly; the exact-set test is correct; the TOML comment is honest. One P3: the chunked_run line in this log was inaccurate. The parent fixed it inline and re-ran verification: 619 passed, ruff clean, mypy clean, h5 sha256 unchanged. Accepted (a log-only change, no re-review needed).

## Blockers / Challenges

None

## Commits

- 8c8eb90 fix: catalase k is the per-reaction value 6.62e7, used as given (ticket 04)

## Time Spent

9m (ticket-start.js to ticket-complete.js)
