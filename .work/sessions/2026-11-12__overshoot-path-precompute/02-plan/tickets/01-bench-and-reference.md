# Ticket 01: bench-and-reference

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `examples/bench_overshoot_pulsed.py` runs the scout scenario (PulsedBeam 1e4 Gy/s peak, 10 ms period, 10 µs ON, 10 Gy, chunked run) on the pure-water config, prints wall time, chunk count and the cumulative time of the overshoot path (cProfile), and follows the module banner layout
- [ ] `tests/simulation/test_overshoot_equivalence.py` holds a frozen copy of the current `derive` / `jacobian` math (overshoot path and helpers) and checks the live functions match it within rtol 1e-12 (with a small atol) on random states with species below 0 and exactly 0, for a config with negative products and enzymatic reactions
- [ ] Baseline bench numbers written in the ticket commit log
- [ ] Full pytest, ruff and mypy pass

**Files to Touch:**
- `examples/bench_overshoot_pulsed.py`
- `tests/simulation/test_overshoot_equivalence.py`

**Verification Step:**

Run:
```bash
python -m pytest tests/simulation/test_overshoot_equivalence.py -q && python examples/bench_overshoot_pulsed.py
```

Expected:
Tests pass; the bench prints the wall time (about 7.5 s) and the overshoot share (about 38%).

**Notes:**

Look at `examples/bench_fenton_stall.py` and `tests/simulation/test_overshoot_relaxation.py` / `test_jacobian.py` for how configs and states are built. Use a fixed RNG seed.
