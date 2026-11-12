# Ticket 01: map-helper

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `compute_map(n=12, ...)` returns dose-per-pulse axis, average-dose-rate axis, ROOH at max dose and ROOH at end time (µM), NaN for infeasible/failed points.
- [ ] Runs points in a `ProcessPoolExecutor`; results cached in an `.npz` keyed by grid parameters; second call loads the cache.
- [ ] Module follows repo banner layout, uses logging not print.

**Files to Touch:**
- `docs/issues/fenton/fenton_map.py`

**Verification Step:**

Run:
```bash
cd docs/issues/fenton && python -c "import fenton_map as m; r=m.compute_map(n=2, period=60, cache=None); print({k: getattr(v,'shape',v) for k,v in r.items()})"
```

Expected:
Prints arrays of shape (2,) axes and (2, 2) maps with finite values for feasible points, no exception.

**Notes:**

Build the UnitCell from `Configuration_Bordet_Emballement.toml` per worker; build a `PulsedBeam` with `max_dose=` by keyword (see CLAUDE.md). ROOH at max dose = value at the time of the last pulse end (dose reaches 10 Gy); ROOH at end = final value.
