# Session: paper-reproduction-notebook

**Date:** 2026-11-12T14:28:40.875Z
**Status:** Grill phase complete

## Problem Statement

radiopyo is the code behind Tonneau et al. 2025 (Phys. Med. Biol. 70 235021, `docs/literature/Tonneau_2025_Phys._Med._Biol._70_235021.pdf`), but nothing in the repo reproduces the paper's main findings. We need a notebook in `examples/tonneau_2025_pmb/` that reproduces them (the ΔO2 and net H2O2 yield maps, the gross production/consumption bars of O2 and H2O2, the time evolution of O2, H2O2 and pH) and whose data become golden data guarded by a test.

## Context & Constraints

- **Current behavior:** The model runs the paper's scenario: the bundled pure-water config (`PURE_WATER_CONFIG_FILE`), LET 1.17 keV/µm (Boscolo G-values match table 3 to 3 digits), `variable_GValues` applies GF(DR) from the peak dose rate, `SinglePulseBeam` + `prepare_single_pulse`. A 30 Gy, 500 s run takes ≈0.5 s. At 5 % pO2 the current code gives final H2O2 3.06 / 1.72 / 1.05 µM at 0.1 / 1e3 / 1e7 Gy/s (paper ≈3.3 / 1.7–1.8 / 1.1), minimum pH 5.86 at 1e7 Gy/s (paper 5.8), ΔO2 1.06 µM at 0.1 Gy/s (paper 0.8–1.2). The only example notebook is `examples/notebook_basic_usage.ipynb` (outputs stripped).
- **Pain point:** The paper's figures cannot be regenerated from the repo, and no test would catch a code change that moves the paper's results.
- **Dependencies:** Keep the existing golden test rules (`tests/golden/record_golden.py`, re-record only for an approved result-changing fix; existing golden files untouched). `MIN_BEAM_ON_TIME` = 1 µs: at 1e7 Gy/s a dose below 10 Gy raises `ShortBeamOnError`. The GF fit in code is 0.972·DR^-0.1097, the paper prints 0.982·DR^-0.09: the code fit is used as-is (ADR 0009). matplotlib is not a package dependency.
- **Tech stack:** Python 3.11, radiopyo, numpy, scipy (LSODA), pandas, matplotlib (notebook only), Jupyter/nbconvert (radiobio env), pytest.

## Success Metrics

- `examples/tonneau_2025_pmb/tonneau_2025_pmb.ipynb` runs top to bottom (nbconvert) with `RECOMPUTE = False` in seconds from the golden file and draws figs 4, 5, 6, 7, 8, 9 and 10 (paper layout: 2×2 CDR/UHDR × hypoxia/physioxia for the time series, contour maps with pO2 % and µM axes, bars with net marker) plus a paper-vs-radiopyo comparison table; with `RECOMPUTE = True` it re-runs the grid through `paper_runs.py` and reports the difference to the golden file.
- `tests/golden/data/paper_tonneau_2025.npz` is written by `record_golden.py --only paper`, and a new fast golden test re-runs a few grid points and two time series through `paper_runs.py` and matches the file (rtol ≈ 1e-6), in a few seconds.
- `python -m pytest`, `ruff check radiopyo` and `mypy radiopyo` stay green; the notebook is committed with outputs stripped.

## Architecture & Approach

- `examples/tonneau_2025_pmb/paper_runs.py` (no matplotlib, module banner layout): the scenario (paper tables 3 and 5: pure-water config, LET 1.17, `variable_GValues` on, one `SinglePulseBeam` per run, `prepare_single_pulse` over (0, 500 s), max ON step = pulse / 500, max final step 0.1 s, [O2] = pO2 % × 1300 µM, pH 7, 25 °C), the grid, one-run analysis (final concentrations, ΔO2, gross production/consumption of O2 and H2O2 via `integrate_production_of` / `integrate_consumption_of`, time series of O2, H2O2, pH resampled on a fixed log time grid), the net yield (slope of a linear fit of [H2O2] at 500 s over 10/20/30 Gy, in #/100 eV), the full-grid driver and npz save/load.
- Grid: 19 dose rates (0.1 to 1e7 Gy/s at two per decade, plus 50 and 5e6) × 19 pO2 values (0 to 21 %, including the bar values 0, 0.1, 0.5, 2, 4, 6, 12, 16, 21 % and 1 % / 5 %) × doses 10/20/30 Gy (≈1,100 runs, ≈10 min). Bars use the 11 dose rates of table 5 × the 9 bar pO2 at 30 Gy; time series the 11 dose rates × 1 % and 5 %.
- Notebook `tonneau_2025_pmb.ipynb`: `RECOMPUTE = False` loads the golden npz, `True` re-runs the grid; draws figs 4/5/6/7/8/9/10 and an informational comparison table of the values quoted in the paper text (no assertion). Notes the GF fit gap (ADR 0009).
- Golden: `tests/golden/data/paper_tonneau_2025.npz`, new `--only paper` part of `record_golden.py`, new test loading `paper_runs.py` by path (importlib) and re-running a small subset.
- Packaging/docs: optional `examples` extra (matplotlib), README examples section, changelog fragment.

## Assumptions & Trade-offs

- The code's GF fit is used, not the paper's printed fit (ADR 0009): the golden data record what radiopyo produces today; the figures will differ slightly from the published ones.
- The paper's numbers are compared informationally only; the model changed since publication (clipping, overshoot relaxation…), the golden test is what catches drift.
- The net yield fit uses 10/20/30 Gy (not down to 0.1 Gy as in the GF calibration) so every dose rate is feasible under the 1 µs minimum pulse.
- Not done: fig 3 (GF calibration), fig 11 (reaction yields), Blain 2022 experimental overlay, a configurable GF fit, committing the executed notebook, and the user's uncommitted README / docs / logo / `docs/literature` changes.

## Open Questions

N/A

## Notes

- Text and figures of the paper were extracted to `.scratch/paper/` with PyMuPDF installed into `.scratch/pdftools` (not a project dependency).
- On this machine the env is `radiobio_env` (`C:\Users\rtonneau\AppData\Local\miniconda3\envs\radiobio_env`), not the `E:\miniconda3\envs\radiobio` path in CLAUDE.md; radiopyo is not installed in it, so scripts run with `PYTHONPATH=.`.
- Timing probe: `.scratch/tests/2026-11-12__paper-reproduction-notebook/probe.py`.
