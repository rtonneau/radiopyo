# Ticket 01: single-pulse-trans-step

**Status:** ✅ Done

## Local Test Result

`python -m pytest tests/simulation -q`: 276 passed in 66.54s.

## Review Notes

Base 54d05c72a997e0c5641a36756cb687f9d9e7bee3, clean tree, inline. AC1: pre-tm OFF chunks (leading OFF, gaps) now use `max_step_size_trans`; the transition keeps trans, the final OFF keeps final. AC2: docstrings of `max_step_size_trans` / `max_step_size_final` updated. AC3: `test_single_pulse_off_chunks_before_tm_use_the_transition_max_step` (pulsed beam starting at 1e-4 s, 3 pulses). Diff inspected (one code line, docstring, test); localized low-risk change, no independent review.

## Blockers / Challenges

None

## Commits

- 6a8670f fix: prepare_single_pulse OFF chunks before tm use max_step_size_trans (ticket 01)

## Time Spent

2m (ticket-start.js to ticket-complete.js)
