# Session: beam-on-time-min-and-single-pulse-off-step

**Date:** 2026-11-11T10:30:36.124Z
**Status:** Grill phase complete

## Problem Statement

1. A beam ON interval shorter than `T_ZERO` (1e-7 s) cannot be computed correctly (the initial burst covers the whole span, `_apply_initial_burst` raises a bare `ValueError("Too Short beam pulse")`). Be more conservative: a configured beam ON time < 1e-6 s must raise a clear error, and short ON spans that come from a partial or clipped pulse must be handled explicitly.
2. In `UnitCell.prepare_single_pulse`, the beam OFF chunks before the time of max dose use `max_step_size_final`; they should use `max_step_size_trans`.

## Context & Constraints

- **Current behavior:** `UnitCell._apply_initial_burst` skips `T_ZERO` and raises `ValueError` when the ON span is not longer than `T_ZERO`; no check exists when a beam is built. `prepare_single_pulse` maps OFF chunks before tm to `max_step_size_final` (documented as such).
- **Pain point:** a short `on_time` fails late (inside a chunked run) with an unclear message; a last pulse cut by `max_dose` or a pulse clipped by `t_span` can fail the same way although the beam is valid; pulses shorter than `T_ZERO` were a recorded follow-up (CHANGELOG Known issue, lsoda-stalls session).
- **Dependencies:** keep `T_ZERO` = 1e-7 and the normal path for ON spans >= 1e-6 s unchanged (golden data must not change); keep `ValueError` catchable (new exception subclasses it); HDF5 files written by older versions must keep loading; ruff and mypy at zero.
- **Tech stack:** Python 3.11+, numpy, scipy `solve_ivp` (LSODA), pytest, h5py.

## Success Metrics

- Building a `PulsedBeam` / `MacroPulsedBeam` with `on_time` or `micro_on_time` < 1e-6 s, or a `ConstantBeam` / `SinglePulseBeam` whose single ON span (`time_at_max_dose - start_time`, finite `max_dose`) is < 1e-6 s, raises `ShortBeamOnError` (subclass of `BeamError` and `ValueError`); on_time == 1e-6 s is accepted.
- A short ON span (< `MIN_BEAM_ON_TIME`: last pulse cut by `max_dose`, micro-pulse cut by its train end, pulse clipped by `t_span`) in `prepare_chunked_run`, `prepare_single_pulse` or `run()` is not integrated: its full radiolytic yield (burst over the span length) is added analytically, one result point is stored at the span end, a `logger.warning` is logged and `(start, stop)` is appended to `ODEResult.short_on_spans`; the delivered dose matches the expected total.
- `prepare_single_pulse` OFF chunks before tm carry `max_step_size_trans`; the final OFF chunk keeps `max_step_size_final`.
- Full pytest suite unchanged except new tests (the two known golden time-grid failures aside); ruff and mypy zero.

## Architecture & Approach

- `UnitCell.MIN_BEAM_ON_TIME = 1e-6` (class constant next to `T_ZERO`), backed by a module constant in `physics` read by the beam checks.
- New `ShortBeamOnError(BeamError, ValueError)` in `physics/exceptions.py`.
- Build-time check in the beam constructors / `update_beam_param` (pulsed: on_time, micro_on_time; constant/single pulse: ON span with finite max_dose).
- Run-time: ON spans shorter than `MIN_BEAM_ON_TIME` are short ON spans: burst over the span length (`compute_initial_cc_before(t0=span)`), no solve_ivp call, one point at the span end, warning, recorded in `ODEResult.short_on_spans` (in-memory only, concatenated by `__add__`, empty for results loaded from HDF5).
- `prepare_single_pulse`: `max_step_size_on if on else max_step_size_trans` before tm; docstring updated.
- Docs: glossary term *Short ON span*, ADR 0008, CHANGELOG fragment, README / CLAUDE.md, the "follow-up pending" note replaced.

## Assumptions & Trade-offs

- A short ON span is approximated (yield delivered instantly at the span scale) rather than refused: the user chose a warning plus result flag over an error for partial/clipped pulses. Chemistry inside that span is not resolved.
- The threshold is a class constant (overridable by assignment), not a run argument.
- `short_on_spans` is not persisted to HDF5 (no file-format change).
- Not done: changing `T_ZERO`, resolving pulses between 1e-7 and 1e-6 s by integration.

## Open Questions

None.

## Notes

Sized Standard (light). User answers: check at both build and run time; partial last pulse -> warning (not error) with burst over span; clipped spans treated the same; logger warning + result flag; new exception subclass; class constant. User chose `/gps auto` after approval.
