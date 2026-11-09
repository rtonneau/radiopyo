# Ticket 01: golden-outputs

**Status:** ✅ Done

## Local Test Result

`python -m pytest tests/golden -q` -> `5 passed in 22.31s` (package unmodified, goldens recorded with PYTHONHASHSEED=1).
Also run under 25 different PYTHONHASHSEED values: 25/25 pass. Recorder run a second time with PYTHONHASHSEED=7 (--force): tests pass at seeds 1, 7 and 99; the committed goldens (seed 1) pass at seed 7.
`python -m pytest tests/chemistry` still fails at collection (RawSpecies import), as expected and unchanged.
Evidence: .scratch/tests/2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian/evidence/ (final_pytest.txt, seeds_final.txt, recorder_two_seeds.txt, spread.txt).

## Review Notes

- Parent: base 0da0d7ac6fe7d3b4ea144e3250077d674f8d4894; pre-existing changes: only gps session state under .work/ (parent-owned). Package code at 2261b92 (unchanged).

Files added (all unstaged, nothing under radiopyo/ touched): tests/golden/test_golden.py, tests/golden/record_golden.py, tests/golden/data/ (derive.json, runs.npz, beams.npz, golden_results.h5, hdf_expected.npz, .golden_results/*.json; ~304 kB total).
AC check: goldens for derive() at 4 states x 4 cases (beam on/off, no beam, force_positive) on both bundled configs; run(t_span=[1e-7,10]) on LABARBE_CONFIG_FILE; chunked pulsed run on PURE_WATER_CONFIG_FILE with LET=2.5 (plus the Boscolo G-values); reconstruct_reactions() rates; make_matrix of Constant, SinglePulse, Pulsed (2 variants) and MacroPulsed beams with finite max_dose; HDF5 written by SimCollection/H5Manager and read back in the test (on a tmp copy). All keyed by species label or str(reaction), never by index. derive() is called only in the _derive_adapter function. Tolerances are named constants with comments (algebra rtol 1e-12; run rtol 1e-6 with absolute floors: 1e-10 uM for species, 1e-7 of the largest rate for rates).
Changed after review: first version failed under other hash seeds (linear resampling of raw solver steps, per-row floors too tight for HO3 rate), replaced by t_eval dense-output comparison, chunk-boundary points only, and a global-scale rate floor (measured spread <= 3.5e-10 of global max rate vs 1e-7 allowed).

- Parent acceptance: independent reviewer skipped — test-only, no package code touched; parent re-ran `python -m pytest tests/golden -q` (5 passed), checked data size (304K) and tolerance constants. Model used: sonnet (requested medium). Noted findings (catalase k doubling on HDF5 json round trip; run point count varies with PYTHONHASHSEED) feed tickets 14 and 16.

## Blockers / Challenges

1. Non-reproducibility of full runs: with the unmodified code, PYTHONHASHSEED changes the species order, hence the summation order, hence the LSODA step sequence. The number of points of the Labarbe run varies (~2294-2448) and raw points cannot be compared. Hence the t_eval/chunk-edge design and the n_points tolerance (15%).
2. K-reaction as_label() is index based (KR_<n>), so it changes with the hash seed; reactions are keyed by str(reaction), checked unique (85 Labarbe, 73 pure water).
3. SimIdentifier.to_string uses field.type(...), so test_golden.py must not use `from __future__ import annotations`.
4. Pre-existing bug pinned, not fixed: the json round trip of H5Manager (read_ODEResult) doubles the catalase k_value (1.324e8 -> 2.648e8), so rates recomputed on a result read back from HDF5 are 2x for `2H2O2 -- catalase >> O2 + 2H2O`. The golden records the rebuilt rates on purpose, so a fix will flip that test and must re-record deliberately.
5. Pre-existing bug recorded as is: UnitCell.run() passes force_positive in the beam_ON slot of _raw_run, so the initial radiolytic burst is skipped in plain runs.
6. H5Manager opens the file in r+ mode, so the HDF test works on a tmp_path copy of tests/golden/data.

## Commits

- dba163a test: record golden outputs of current behaviour (ticket 01)

## Time Spent

37m (ticket-start.js to ticket-complete.js)
