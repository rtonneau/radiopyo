# Ticket 01: golden-outputs

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] A recorder script (`tests/golden/record_golden.py`) generates reference data from the current code, run once at 2261b92 before any package change.
- [ ] Goldens cover: `derive()` at several fixed states (beam on and off) for both bundled configs; `run(t_span=[1e-7, 10])` on `LABARBE_CONFIG_FILE`; a chunked pulsed run on `PURE_WATER_CONFIG_FILE` with LET set; `reconstruct_reactions()`; beam `make_matrix` output for Constant, SinglePulse, Pulsed and MacroPulsed beams (finite `max_dose`); an HDF5 file written by H5Manager/SimCollection with the old code.
- [ ] All goldens are keyed by species/reaction label, not index.
- [ ] `tests/golden/test_golden.py` compares the current code to the goldens and passes on the unmodified package.
- [ ] Tolerances are explicit constants in the test file, with a comment explaining each.

**Files to Touch:**
- `tests/golden/record_golden.py`
- `tests/golden/test_golden.py`
- `tests/golden/data/*` (npz/csv/h5 reference files)
- `tests/conftest.py` (only if shared fixtures are needed)

**Verification Step:**

Run:
```bash
python -m pytest tests/golden -q
```

Expected:
All golden tests pass (tests/chemistry may still fail at collection until ticket 02; run tests/golden alone).

**Notes:**

Record before editing anything under `radiopyo/`. Keep reference files small (subsample time points if needed). `run()` today skips the initial burst because of the argument-order bug: that is the behaviour to record. Run the recorder twice with different `PYTHONHASHSEED` values to confirm label-keyed comparison is stable.
