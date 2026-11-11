# Implementation Plan

**Session:** beam-on-time-min-and-single-pulse-off-step
**Date:** 2026-11-11T10:36:30.181Z
**Estimated effort:** 3-5 hours

## Strategy

Smallest change first (the `prepare_single_pulse` max step), then the build-time check with the new exception and the threshold constant, then the run-time handling of short ON spans (burst over the span, warning, `ODEResult.short_on_spans`) together with the docs and the ADR.

## Tickets Overview

- **Ticket 1:** `prepare_single_pulse` OFF chunks before tm use `max_step_size_trans`.
- **Ticket 2:** `MIN_BEAM_ON_TIME` = 1e-6 s and `ShortBeamOnError`, raised when a beam with a too-short configured ON time is built.
- **Ticket 3:** run-time short ON spans (< 1e-6 s, partial / clipped pulses): burst over the span, warning, `ODEResult.short_on_spans`; docs, glossary, ADR 0008, CHANGELOG fragment.

## Sequencing Rationale

Ticket 1 is independent and trivial. Ticket 3 relies on ticket 2's constant and on the guarantee that any short ON span at run time comes from a partial or clipped pulse, not from a bad beam.

## Risks & Mitigation

- **Risk:** existing tests or golden beams use `on_time` < 1e-6 s → **Mitigation:** grep tests/examples first; adjust tests only where they pin the old behaviour, never golden data without approval.
- **Risk:** the one-point synthetic result of a short span breaks `ODEResult` bookkeeping (`n_points`, `step_dose_rates`, `burst_dose_rates`, `reconstruct_reactions`) → **Mitigation:** build it through `append_sol` with a minimal solution object and test reconstruct_reactions / dose integral on a run with a short span.
- **Risk:** floating point makes a full-length ON chunk look slightly shorter than on_time → **Mitigation:** compare with a small relative tolerance (only spans clearly below the threshold are short).

## Assumptions

- `T_ZERO` stays 1e-7 and spans >= 1e-6 s keep today's path, so golden data does not change.
- `short_on_spans` is in-memory only (no HDF5 change).
