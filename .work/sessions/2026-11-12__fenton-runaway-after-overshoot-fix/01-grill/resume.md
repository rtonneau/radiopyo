# Session: fenton runaway after overshoot fix

**Date:** 2026-11-12T08:07:21.120Z
**Status:** Grill phase complete

## Problem Statement

A user report (notebook and config in `docs/issues/fenton/`): adding Fe2+ = 2.2e-8 M makes ROOH keep growing to the end of a 1200 s run. Some beams show "islands" (a final ROOH much lower than neighbouring beams), and one beam (`UHDR_Freeze`) froze with `lsoda: Repeated convergence failures`. Question: was this already fixed by the PR "Fix LSODA stalls near zero concentration (overshoot relaxation, exact beam chunks, solver failure error)", and is the remaining growth usable?

## Context & Constraints

- **Current behavior:** On current main the reported case runs with no LSODA warning. `UHDR_Freeze` takes about 1 s and gives ROOH = 3.01 µM, in line with its neighbours (2.79 at 0.12 Hz, 3.02 at 0.14 Hz). ROOH still grows linearly, about 1.2e-3 µM/s at 1200 s (3e-6 without Fe).
- **Pain point:** The linear growth comes from the model inputs: catalase = 0 keeps H2O2 at about 1.06 µM, Fe2+ is a fixed species and is never consumed, and the ROO_r + RH -> ROOH + R_r chain multiplies the Fenton OH flux (2.3e-5 µM/s) by about 50. While checking this, a real bug turned up: the bundled `reactions_biology.toml` SOD reaction used `O2_r_minus`, a species nothing produces (water config: `O2_minus`).
- **Dependencies:** The golden derive/runs data change with the SOD fix; the HDF5 golden file must stay untouched.
- **Tech stack:** radiopyo, scipy LSODA, TOML configs.

## Success Metrics

- The reported case reproduces with no solver failure and no island on current main (verified).
- `docs/issues/fenton/ANALYSIS.md` explains the drift mechanism, the config levers that stop it (with measured numbers) and answers the "use the runaway anyway?" question.
- The bundled SOD reaction reads `2 O2_minus -> H2O2 + O2`; derive/runs golden data re-recorded; the test suite has no new failures.

## Architecture & Approach

No solver change. Verification and lever-variant runs as scratch scripts. One data fix in `radiopyo/data/reactions_biology.toml` plus re-recorded golden data and a CHANGELOG Fixed entry. A write-up in `docs/issues/fenton/ANALYSIS.md`.

## Assumptions & Trade-offs

The reproduction uses the current `advanced_water_radiolysis.toml` (the user's run used an older `basic_water_radiolysis.toml`). The Fe3+ + O2_minus -> Fe2+ + O2 recycling k (1.5e8) and GSH = 6.5e-3 M are test choices. The user's copy in `docs/issues/fenton` is not edited. Not done: a steady-state warning in radiopyo.

## Open Questions

The CONV lever variants were still running when the session closed; only the CONV baseline (ROOH 6.61 µM at 1200 s, slope 1.1e-3 µM/s) is in the write-up. Four tests (`test_short_on_span.py` x3, `test_run_options.py` x1) fail on clean main too, probably the pandas of the Python used; not investigated.

## Notes

Recommendation given to the user: do not use the runaway runs as results, because their final ROOH depends on the stop time. Use catalase on (or a fuller antioxidant set) to get a settled value, and scale a tolerable Fe2+ by the healthy/tumour ratio.
