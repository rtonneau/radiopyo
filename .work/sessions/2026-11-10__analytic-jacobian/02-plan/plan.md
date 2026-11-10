# Implementation Plan

**Session:** analytic-jacobian
**Date:** 2026-11-10T10:55:16.444Z
**Estimated effort:** ~1 day

## Strategy

Block depleted enzymatic reactions first so the Jacobian targets the final `derive`; then write and prove the Jacobian alone; then wire `use_jac` (default still False); then benchmark and, behind a gate, switch the default to True; docs last.

## Tickets Overview

- **Ticket 1:** Enzymatic depletion block (F9): `EzMatrices.negative_products`, shared block in `reaction_rates`.
- **Ticket 2:** Public pure dense `jacobian()` validated against a finite-difference Jacobian of `derive`.
- **Ticket 3:** Wire `use_jac` to `solve_ivp(jac=)` for LSODA/BDF/Radau only; remove NotImplementedError; default still False.
- **Ticket 4:** Benchmark analytic vs FD; gate; default `use_jac=True`; golden tolerance check.
- **Ticket 5:** README, CHANGELOG, CLAUDE.md.

## Sequencing Rationale

Each ticket builds on the previous: the Jacobian must match the final rates (01), the run wiring needs the function (02), the default switch needs the wiring and the benchmark (03, 04), docs describe the final state.

## Risks & Mitigation

- **Risk:** noisy FD comparison (µmol/L scaling, k up to ~1e10) → **Mitigation:** central differences, step scaled to |y|, per-entry relative tolerance plus an absolute floor relative to the matrix norm.
- **Risk:** analytic Jacobian slower or less robust than FD → **Mitigation:** benchmark gate in ticket 04: stop and ask the user before changing the default.
- **Risk:** goldens drift with the default True → **Mitigation:** compare to unchanged data with a tolerance written in the test; re-record only with approval.

## Assumptions

- Dense Jacobian is enough (bundled systems ≤ ~40 species); LSODA requires dense.
- `tests/golden/data` stays unchanged.
