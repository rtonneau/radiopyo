# Implementation Plan

**Session:** fenton-fix-showcase-notebook
**Date:** 2026-11-12T08:48:55.934Z
**Estimated effort:** 3 hours

## Strategy

First an importable helper module that computes (and caches) the dose-per-pulse × average-dose-rate grid, since a process pool on Windows cannot pickle functions defined in a notebook. Then the notebook itself, built and executed through the Jupyter MCP, using the helper for both maps.

## Tickets Overview

- **Ticket 1:** `docs/issues/fenton/fenton_map.py`: parallel, cached computation of ROOH at max dose and at end time over a log grid.
- **Ticket 2:** `docs/issues/fenton/Fenton_Fixed_Showcase.ipynb`: recap, time-course runs, frequency scan, two maps; executed with outputs.

## Sequencing Rationale

The notebook's map sections depend on the helper module.

## Risks & Mitigation

- **Risk:** 144 runs × 1200 s too slow → **Mitigation:** process pool, cache `.npz`, allow a smaller grid via argument; mask failed/unreached points.
- **Risk:** Very low dose per pulse with high average dose rate implies frequency > 1/on_time → **Mitigation:** mask infeasible points (duty cycle ≥ 1).

## Assumptions

- Same TOML files used unchanged; beam built from the config's `conv` beam on_time (1.2e-6 s) and LET 0.2116.
- Peak dose rate = dose_per_pulse / on_time, frequency = avg_dose_rate / dose_per_pulse, max_dose 10 Gy, period 1200 s.
