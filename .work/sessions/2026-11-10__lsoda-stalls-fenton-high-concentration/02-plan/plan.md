# Implementation Plan

**Session:** lsoda-stalls-fenton-high-concentration
**Date:** 2026-11-11T07:41:39.801Z
**Estimated effort:** 1 day

## Strategy

First make failures loud (cheap, independent, protects every later measurement). Then the overshoot relaxation in derive/jacobian (the main stall cause), then the beam-side fix (exact pulse-end cut + constant dose rate per chunk). Last, the regression benchmark on the reproduction and the golden before/after table, which pauses for the user's approval before re-recording.

## Tickets Overview

- **Ticket 1:** a failed solve_ivp raises `SolverFailedError` instead of continuing silently.
- **Ticket 2:** overshoot relaxation of dynamic species below 0 in `derive`, exact in `jacobian` (ADR 0007).
- **Ticket 3:** chunks cut exactly at pulse end (no GROWTH_TIME extension) and a constant dose rate per chunk.
- **Ticket 4:** regression benchmark on the Fenton reproduction, golden before/after table, re-record after approval, docs.

## Sequencing Rationale

Ticket 1 makes any stall/failure in tickets 2–4 visible. Ticket 2 and 3 are independent code areas but both needed for the benchmark of ticket 4; ticket 2 first because it is the larger result change. Ticket 4 needs both.

## Risks & Mitigation

- **Risk:** relaxation term makes the analytic Jacobian inexact → **Mitigation:** include c·∂J_jj/∂y_k and test against finite differences on states with negative species.
- **Risk:** relaxation slows derive (extra derivative work) → **Mitigation:** compute only the negative columns; benchmark wall time against the prototype.
- **Risk:** constant per-chunk dose rate wrong for beams whose dose rate varies inside an interval → **Mitigation:** only for chunked runs, whose intervals come from timer transitions (piecewise-constant); test every beam type.
- **Risk:** golden drift larger than solver tolerance → **Mitigation:** before/after table and explicit user approval before re-recording.

## Assumptions

- All beam types are piecewise constant between `make_matrix` transitions.
- The user's real runs use defaults with Fe2+ fixed, matching `fe.toml`.
