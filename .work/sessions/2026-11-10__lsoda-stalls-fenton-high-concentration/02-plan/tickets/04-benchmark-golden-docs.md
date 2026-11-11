# Ticket 04: benchmark-golden-docs

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] A benchmark script under `tests/` or `examples/` (not in the default pytest run, or marked slow) runs the Fenton reproduction grid (fe config, pulse dose {0.1, 1, 2} Gy × {1e3, 1e5, 1e7} Gy/s with on_time ≥ T_ZERO, plus Fe 8.9e-7): every point succeeds under 60 s with defaults; ROOH within 1e-6 relative of BDF/LSODA-atol-1e-10 references where those succeed.
- [ ] A fast regression test (one stalling point from the reproduction, bounded wall time) in the default suite.
- [ ] Before/after table of golden differences per species shown to the user; golden run data re-recorded (`record_golden.py --force --only <part>`, PYTHONHASHSEED fixed) only after approval.
- [ ] README/CLAUDE.md/CHANGELOG updated (overshoot relaxation, chunk boundaries, SolverFailedError); follow-up recorded for pulses shorter than T_ZERO.

**Files to Touch:**
- `tests/simulation/test_fenton_stall.py`
- `examples/bench_fenton_stall.py`
- `tests/golden/data/`
- `README.md`
- `CLAUDE.md`
- `CHANGELOG.md`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo && python examples/bench_fenton_stall.py
```

Expected:
all tests pass, zero ruff/mypy errors, benchmark prints every grid point ok=True under 60 s.

**Notes:**

The golden re-record needs explicit user approval (grill decision): stop and ask with the before/after table. Reproduction configs are in `.scratch/tests/2026-11-10__lsoda-stalls-fenton-high-concentration/` (fe.toml, bundled_fe.toml); copy what is needed into the test tree.
