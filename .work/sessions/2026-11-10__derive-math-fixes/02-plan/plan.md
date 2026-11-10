# Implementation Plan

**Session:** derive-math-fixes
**Date:** 2026-11-10T17:46:38.034Z
**Estimated effort:** 1 day

## Strategy

Remove `force_positive` first so `derive` and `jacobian` are simpler before the math changes, then introduce clipped concentrations (the smooth switch reads c⁺), then the smooth depletion switch, then the catalase constant on its own so its effect on results is visible alone, then the polynomial-k warning, then docs. Each result-changing ticket (02, 03, 04) shows before/after numbers of the golden comparison and re-records only the failing golden parts (`tests/golden/record_golden.py --force --only <part>`, `PYTHONHASHSEED` fixed); the old `tests/golden/data/golden_results.h5` must keep loading.

## Tickets Overview

- **Ticket 1:** remove `force_positive` and `zero_depleted_consumption` from the API (RunOptions, UnitCell run/prepare_*, derive, jacobian) and the tests.
- **Ticket 2:** rate laws and their derivatives read clipped concentrations max(c, 0); overshoot and mass-balance tests; golden re-record if needed.
- **Ticket 3:** smooth depletion switch c⁺/(c⁺+ε) replaces the hard depletion block for k and enzymatic reactions; Jacobian includes the switch derivative; Jacobian vs finite-difference test; golden re-record if needed.
- **Ticket 4:** catalase k_value 1.324e8 -> 6.62e7 in reactions_biology.toml; golden re-record.
- **Ticket 5:** warning at build for a non-scalar k on a single-reactant n > 1 reaction; caplog test.
- **Ticket 6:** README, CLAUDE.md, CHANGELOG fragment, ADR 0003 status note.

## Sequencing Rationale

Ticket 1 removes code paths (force_positive rows in the Jacobian) that tickets 2 and 3 would otherwise have to keep consistent. Ticket 3 builds on the clipped concentrations of ticket 2. Ticket 4 is isolated so its golden diff is attributable. Ticket 5 is independent. Docs last so they describe the final behaviour.

## Risks & Mitigation

- **Risk:** golden diffs mix several changes and hide a bug → **Mitigation:** one result-changing ticket each, before/after numbers recorded in the ticket log before any re-record.
- **Risk:** the smooth switch changes solver behaviour (stiffness, step counts) on bundled configs → **Mitigation:** Jacobian vs finite-difference test; compare run statistics (nfev, njev) before/after in ticket 3's log.
- **Risk:** re-recording breaks loading of the old HDF5 golden file → **Mitigation:** never rewrite `golden_results.h5`; only re-record `derive`, `runs`, `beams`, `hdf` expectation files as needed and keep the HDF loading test.
- **Risk:** the catalase value reading (6.62e7 per reaction) is wrong → **Mitigation:** recorded as an open question in the grill resume and ADR 0006.

## Assumptions

- The baseline is green (551 tests, ruff and mypy at zero) on main at the session start.
- Hard removal of `force_positive` (breaking API change) is approved; no deprecation period.
- ε = DEPLETION_LIMIT = 1e-12 mol/L, not configurable; `DEPLETION_LIMIT` keeps its name.
