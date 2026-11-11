# Ticket 04: benchmark-golden-docs

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo && python examples/bench_fenton_stall.py` run as separate commands (`.scratch/tests/2026-11-10__lsoda-stalls-fenton-high-concentration/evidence/04-verification.txt`, `.scratch/tests/2026-11-10__lsoda-stalls-fenton-high-concentration/evidence/04-bench.txt`).

- pytest: 702 passed, 2 failed (expected, golden not re-recorded: `test_beam_matrices_match_golden`, `test_chunked_pulsed_run_matches_golden`, time grid only; table in `.scratch/tests/2026-11-10__lsoda-stalls-fenton-high-concentration/evidence/04-golden-before-after.md`).
- ruff check radiopyo: All checks passed. mypy radiopyo: Success, no issues in 42 source files.
- Benchmark: 11 grid points, every one ok=True, wall 1.2 s to 27.4 s (slowest: D=0.1 Gy at 1e5 Gy/s 27.4 s, at 1e3 Gy/s 20.7 s), ROOH within 3.4e-8 relative of the LSODA atol-1e-10 references and 2e-9 of the BDF atol-1e-10 references (BDF failed with SolverFailedError at D=2 Gy, 1e7 Gy/s: no reference there). "0 failing point(s)".
- New fast test `tests/simulation/test_fenton_stall.py`: about 8 s; the same point (3 pulses of 0.1 Gy at 1e5 Gy/s) took 81 s on the code of 0cda498 (before the tickets) and 0.9 s now.
- Derive cost (`.scratch/.../derive_cost.py`): 63 us/derive and 210 us/jacobian with all c >= 0; 1075 us/derive and 3131 us/jacobian with 9 species below 0 (about 17x and 15x).

## Review Notes

- AC1: `examples/bench_fenton_stall.py` (outside pytest) runs the grid D {0.1, 1, 2} Gy x DR {1e3, 1e5, 1e7} Gy/s plus Fe 8.9e-7 (D = 1 Gy, 1e3 and 1e5), with the ROOH checks against BDF/LSODA atol 1e-10. Points with on_time <= T_ZERO are skipped (0.1 and 1 Gy at 1e7 Gy/s): 1 Gy at 1e7 has on_time == T_ZERO = 1e-7 and raises "Too Short beam pulse" (the ticket said on_time >= T_ZERO; the check is `>`). Each grid point = 10 Gy in pulses every 0.1 s, then 2 s.
- AC2: `tests/simulation/test_fenton_stall.py` (stalling point, wall limit 30 s, ROOH vs atol 1e-10 run).
- AC3: before/after table done (`04-golden-before-after.md`); NOT re-recorded, waits for approval. Commands: `PYTHONHASHSEED=0 python tests/golden/record_golden.py --force --only beams` then `--only runs`.
- AC4: README (chunks, SolverFailedError, overshoot relaxation, short pulse), CLAUDE.md (derive/jacobian, chunks, baseline line), CHANGELOG (overshoot relaxation bullet, benchmark, derive cost, golden note, follow-up on short pulses in Known issues) updated.

## Blockers / Challenges

None blocking. The reference BDF run fails at D=2 Gy, 1e7 Gy/s (SolverFailedError), so that point has an LSODA reference only. The full benchmark takes about 4.5 min with references (`--no-reference` skips them). Files README/CLAUDE/CHANGELOG are CRLF and were kept CRLF.

## Commits

- c087ebd test: Fenton stall regression, benchmark, golden re-record and docs (ticket 04)

## Time Spent

1h 08m (ticket-start.js to ticket-complete.js)
