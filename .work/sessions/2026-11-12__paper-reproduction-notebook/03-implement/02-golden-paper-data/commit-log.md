# Ticket 02: golden-paper-data

**Status:** ✅ Done

## Local Test Result

- Recording: `PYTHONHASHSEED=0 PYTHONPATH=. python tests/golden/record_golden.py --only paper --workers 10` ran the 1083 runs (19 × 19 × 3) in about 6 min, exit 0 (log: `.scratch/tests/2026-11-12__paper-reproduction-notebook/02-record.log`). It wrote `tests/golden/data/paper_tonneau_2025.npz` (174 kB). Every array is finite.
- Verification step: `python -m pytest tests/golden -q` gives 14 passed in 8.1 s (9 new paper tests).
- `python -m pytest -q` gives 802 passed.
- Sanity check of the recorded grid:
  - ΔO2 maximum 0.091 %pO2 at 3.2e3 Gy/s, 0.25 % pO2;
  - net H2O2 hot spot 0.85 #/100 eV at 1e3–3e3 Gy/s, 0.25 % pO2;
  - maximum 1.05 #/100 eV at 0.1 Gy/s, 21 %; near 0 at 0 % pO2 with CDR;
  - lowest pH minimum 5.80.

  These match the shapes of the paper's figs 5 and 8.

## Review Notes

Base: `5046325` (after ticket 01). There were no pre-existing changes other than gps-owned files.

Checked against the acceptance criteria:
- `record_golden.py` has a new `paper` part (in `PARTS`, so it is included in the default "all") and a `--workers` option (default: half the CPUs). The docstring documents the part, its cost, the same re-record rule, and the new file.
  - `_load_test_module` now takes the file and module name.
  - The derive/runs/beams/hdf code paths are unchanged, and `git status` shows no other golden file modified.
- `tests/golden/test_paper_tonneau_2025.py` holds `PAPER_FILE` and `load_paper_runs()`. The loader loads `examples/tonneau_2025_pmb/paper_runs.py` by path, registered as `paper_runs`, with its folder on sys.path so that worker processes can unpickle it. The recorder reuses both.
- The test re-runs 8 points:
  - (0.1 Gy/s, 1 %) and (1e7 Gy/s, 5 %) at all three doses: net yield, final O2/H2O2/pH min, integrals and time series;
  - (1e3, 21 %) and (1e7, 0 %) at 30 Gy: bar integrals.

  Tolerances are rtol 1e-6, plus atol 1e-9 µM for the near-zero H2O2 consumption at UHDR. A layout test pins the axes, the shapes and the time grid to the module constants, so a grid change without a re-record fails loudly.

## Blockers / Challenges

None.

## Commits

- 3bd6c60 test: golden data of the Tonneau 2025 paper reproduction (ticket 02)

## Time Spent

4m (ticket-start.js to ticket-complete.js)
