# Ticket 04: scaling-benchmark-docs

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `examples/bench_chunked_scaling.py` (module banner layout, `--pulses` list, `--points-on/--points-off` optional) runs config_Labarbe pulsed runs and prints chunks, stored points, wall time and per-chunk cost per pulse count, up to 10k pulses.
- [ ] Measured numbers (at least 100 / 400 / 1600 pulses, and 10k with thinning) recorded in the commit log, showing a flat per-chunk cost.
- [ ] README documents the thinning options and the integral caveat; CHANGELOG (Unreleased) lists the storage fix, derive trim and thinning; CLAUDE.md architecture notes updated (storage, thinning, benchmark).

**Files to Touch:**
- `examples/bench_chunked_scaling.py`
- `README.md`
- `CHANGELOG.md`
- `CLAUDE.md`

**Verification Step:**

Run:
```bash
python examples/bench_chunked_scaling.py --pulses 100 400 1600 && ruff check radiopyo && mypy radiopyo
```

Expected:
Per-chunk cost roughly constant across pulse counts; ruff and mypy clean.

**Notes:**

Same beam as the grill measurement: PulsedBeam peak 1e4 Gy/s, 10 µs ON, 10 ms period, 0.1 Gy per pulse.
