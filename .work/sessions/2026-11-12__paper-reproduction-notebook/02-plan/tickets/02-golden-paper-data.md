# Ticket 02: golden-paper-data

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `tests/golden/record_golden.py` gets a `paper` part writing `tests/golden/data/paper_tonneau_2025.npz` from `paper_runs.py` (full grid); the other parts and files are untouched.
- [ ] A new test (in `tests/golden/`) loads `paper_runs.py` by path, re-runs a few grid points (net yield at ≥ 2 points, bar integrals, two time series) and matches the npz at rtol ≈ 1e-6; it runs in a few seconds.
- [ ] The npz is recorded once and committed; `python -m pytest` passes.

**Files to Touch:**
- `tests/golden/record_golden.py`
- `tests/golden/test_paper_tonneau_2025.py`
- `tests/golden/data/paper_tonneau_2025.npz`

**Verification Step:**

Run:
```bash
python -m pytest tests/golden -q
```

Expected:
All golden tests pass, including the new paper test, in seconds.

**Notes:**

Document the new part in the recorder docstring (same re-record rules). Recording the full grid takes ~10 min: run it in the background.
