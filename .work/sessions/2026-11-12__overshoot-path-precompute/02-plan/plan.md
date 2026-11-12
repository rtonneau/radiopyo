# Implementation Plan

**Session:** overshoot-path-precompute
**Date:** 2026-11-11T23:33:55.122Z
**Estimated effort:** 3-4 hours

## Strategy

First freeze the current behaviour (a reference copy of the solver math in a test, plus a benchmark with baseline numbers), then precompute the per-species overshoot tables on `SimMatrices`, then vectorise the derivative helpers. Each step is checked against the reference at rtol 1e-12 and measured with the bench.

## Tickets Overview

- **Ticket 1:** benchmark script for the scout pulsed scenario and an equivalence test comparing `derive` / `jacobian` to a frozen reference copy of the current code; baseline numbers recorded.
- **Ticket 2:** `OvershootTables` built once in `SimMatrices.__post_init__`; the relaxation terms read it, no more `dataclasses.replace` per call.
- **Ticket 3:** vectorise `depletion_switch`, `depletion_switch_derivatives` and `mass_action_derivatives`; final before/after numbers.

## Sequencing Rationale

The reference test and baseline must exist before any change to prove equivalence and measure the gain. Tables before vectorisation: ticket 2 removes the largest cost (38%) and its new tables feed the helpers vectorised in ticket 3.

## Risks & Mitigation

- **Risk:** a vectorised product silently differs from the loop in an edge case (repeated reactant 2 A, several negative products, species exactly 0) → **Mitigation:** the equivalence test draws random states including exact zeros and negatives on configs with negative products and enzymatic reactions, plus the existing FD-Jacobian tests.
- **Risk:** the tables add build cost to every `SimMatrices` → **Mitigation:** the bench reports the matrix build time; tables are only index arrays and small sub-arrays.

## Assumptions

- `SimMatrices` is only built in `environment_build.py` and by tests through the same constructor; it is not persisted, so an `init=False` field is compatible.
- The reference copy lives under `tests/` only (excluded from ruff/mypy) and is removed by nobody: it documents the pre-refactor math.
