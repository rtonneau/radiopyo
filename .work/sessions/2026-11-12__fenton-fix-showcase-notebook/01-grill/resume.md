# Session: fenton-fix-showcase-notebook

**Date:** 2026-11-12T08:46:51.454Z
**Status:** Grill phase complete

## Problem Statement

The Fenton "runaway" case (`docs/issues/fenton/Example_Emballement.ipynb`) used to freeze (UHDR_Freeze) and produced islands in the ROOH map. These were solver artifacts, now fixed. A new notebook next to the old one should show the fix with the same config files, with clear graphs, plus a simple reproduction of the ROOH map (dose per pulse × average dose rate, 50 µM O2, no GSH, Fe2+ = 2.2e-8 M).

## Context & Constraints

- **Current behavior:** Only the old notebook exists; findings are in `docs/issues/fenton/ANALYSIS.md`.
- **Pain point:** No visual demonstration that the current radiopyo fixes the freeze and islands.
- **Dependencies:** `Configuration_Bordet_Emballement.toml` and this folder's `reactions_biology.toml`, used unchanged (the SOD `O2_r_minus` label bug is only mentioned in markdown).
- **Tech stack:** radiopyo, numpy, matplotlib, concurrent.futures, JupyterLab (built/executed via Jupyter MCP).

## Success Metrics

- New notebook `docs/issues/fenton/Fenton_Fixed_Showcase.ipynb` runs top to bottom without solver failure.
- Section 1: recap of old symptoms (from ANALYSIS.md).
- Section 2: UHDR / CONV / UHDR_Freeze 1200 s runs with ROOH, H2O2, Fe3+ curves and solver stats.
- Section 3: frequency scan 0.10–0.16 Hz, Fe = 0 vs 2.2e-8, smooth (no freeze/islands).
- Section 4: 12×12 map of ROOH at max dose (10 Gy), normalized, dose/pulse 0.03–15 Gy × average dose rate 0.2–1000 Gy/s, log axes, jet-like colormap, label box like the image.
- Section 5: second map, same grid, ROOH at the end simulation time (1200 s).
- Map runs in parallel and cached in a `.npz` next to the notebook so re-runs are instant.

## Architecture & Approach

Notebook only. Grid computation in a helper function run with a process pool; each point builds a PulsedBeam from (dose per pulse, average dose rate) with the config's on_time and LET, runs chunked to 1200 s, and records ROOH at the time 10 Gy is reached and at 1200 s. Results cached to `.npz`; plotting with pcolormesh/contourf on log axes.

## Assumptions & Trade-offs

Old version is not executed (narrative only). Grid is coarse (12×12), not the full-resolution image. Points where 10 Gy cannot be delivered within the period are masked.

## Open Questions

N/A

## Notes

Map 2 (ROOH at end time) added on user request after design review.
