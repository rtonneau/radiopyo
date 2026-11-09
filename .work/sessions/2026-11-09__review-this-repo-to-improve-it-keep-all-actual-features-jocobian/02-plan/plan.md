# Implementation Plan

**Session:** review this repo to improve it. Keep all actual features. Jocobian will be one of the first next job. Just improve flow and readability. Don't break how code is used by user. Other things can be adapted if needed
**Date:** 2026-11-09T20:25:07.087Z
**Estimated effort:** 4-6 days

## Strategy

Lock today's behaviour first (golden outputs recorded from 2261b92 + ported tests), then refactor without changing a single number (each ticket proven against the goldens), solver core first because the Jacobian job builds on it. Only then land the result-changing fixes, one per ticket, each stopping for the user's approval of before/after numbers before goldens are regenerated. Docs close the session. ruff + mypy reach zero in ticket 04 and stay at zero for every later ticket.

## Tickets Overview

- **Ticket 1:** Golden outputs recorded from the current code (derive, run, chunked pulsed run, reconstruct_reactions, beam matrices, old-format HDF5 file) + a pytest that compares against them.
- **Ticket 2:** Port tests/chemistry to the current API; add tests pinning current include-merge and timer behaviour.
- **Ticket 3:** Packaging and repo hygiene: dev extras, untrack artefacts, IBA data out of the wheel, ruff lint section, `__all__`, export `MacroPulsedBeam`, logging without side effects.
- **Ticket 4:** Remove dead code and agreed features (RON parser, `from_ron`, `run_multi_doses`, `e_LET`, …); bring ruff and mypy to zero.
- **Ticket 5:** Solver core: named matrix objects, one pure `reaction_rates()` shared by `derive` and a vectorised `reconstruct_reactions`; stale Jacobian removed.
- **Ticket 6:** UnitCell run flow: `RunOptions`, `_raw_run` split into steps, `prepare_*` deduplicated (today's run() behaviour reproduced explicitly).
- **Ticket 7:** SimEnv: table-driven config readers, one `reactions_involving(...)` with aliases, dead flags/fields removed, crash fixes.
- **Ticket 8:** Split `simulation/utils.py` into units / collections / persistence with re-exports; crash fixes in `ResCollection`/`SimCollection`.
- **Ticket 9:** Chemistry cleanup: shared `from_generic` population, `BiReaction` owns its k-reactions, `get_couple` dedupe, collection index, `None` sentinel, `__eq__`/`__hash__`.
- **Ticket 10:** Split `physics/beam.py` into timer / beam / beam_collection; `BEAM_TYPES` registry; non-mutating `from_dict`; timer dedupe; crash fixes.
- **Ticket 11:** One generic `make_matrix` for all beams driven by timer transitions.
- **Ticket 12:** [result-changing] `run()` keyword passing (initial burst applied), no in-place `y0`, no burst with `use_beam=False`, raise `NoBeamDefinedError`.
- **Ticket 13:** [result-changing] Temperature threaded to fixed species and initial values.
- **Ticket 14:** [result-changing] Include merge rewritten: order-insensitive key, enzymatic/acid-base/constant overrides fixed, later-include-wins, file-relative paths, exact `radiopyo/` prefix.
- **Ticket 15:** [result-changing] Macro beam and timers: pulses-per-train precedence bug, hang without `max_dose`, `stop`, start times, `tot_pulses`, `copy`.
- **Ticket 16:** [result-changing] Deterministic species order.
- **Ticket 17:** Arrhenius (`A`, `Ea` in kJ/mol), capitalised enzyme names, `ODEResult.__add__` returns a new object.
- **Ticket 18:** README and CHANGELOG brought in line with the code; seeds recorded.

## Sequencing Rationale

Tickets 1-2 must exist before any code change, otherwise "no numeric change" cannot be proven. Hygiene and dead-code removal (3-4) shrink the surface and set the lint/type baseline at zero so later diffs stay clean. The solver core (5) comes first among refactors because it is the Jacobian prerequisite and the run flow (6) and results code depend on its matrix objects. SimEnv (7) and utils split (8) follow the solver; chemistry (9) is independent but touches the builders refactored in 7. Beams (10-11) are self-contained; the generic `make_matrix` is separated from the split because it is the riskiest refactor. Result-changing fixes (12-17) come after all refactors so each before/after diff isolates exactly one fix; 14 relies on the chemistry cleanup (9) for a reaction key. Docs (18) describe the final code.

## Risks & Mitigation

- **Risk:** species indices depend on `set` iteration order (varies with `PYTHONHASHSEED`), making goldens flaky → **Mitigation:** goldens store and compare by species/reaction label, never by index; ticket 16 later makes the order deterministic.
- **Risk:** ODE solver runs amplify last-bit rounding changes from refactored arithmetic → **Mitigation:** tight tolerance (rtol 1e-12) for derive, matrices and beam matrices; solver runs compared with a documented tolerance tied to the solver's rtol/atol, and any loosening is stated in the ticket's commit log.
- **Risk:** a refactor silently changes results → **Mitigation:** every Phase B ticket runs the golden test; goldens are only regenerated in result-changing tickets after user approval.
- **Risk:** breaking user imports when moving modules → **Mitigation:** re-export every moved public name from its old module path; a test imports the old paths.
- **Risk:** HDF5 files from the current version stop loading → **Mitigation:** ticket 1 stores a file written by the old code; a test reads it with the new code in every later ticket.

## Assumptions

- Development uses the Anaconda Python on PATH; `python -m pytest`, `ruff`, `mypy` are available.
- Bundled configs (`PURE_WATER_CONFIG_FILE`, `LABARBE_CONFIG_FILE`) are representative enough for goldens, complemented by a small pulsed and a macro-pulsed beam case.
- Before/after reports for result-changing tickets go to the session scratch dir `.scratch/tests/<session>/`.
