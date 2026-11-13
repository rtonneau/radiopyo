# Ticket 01: radkin-min-on-time

**Status:** ✅ Done

## Local Test Result

`cargo test` in C:\DEV\RUSTadkin: all suites ok (18+4+3+6+1+4+4+4+5+2+2+7 passed, 0 failed). Output: .scratch/tests/2026-11-13__relax-beam-on-time-constrain-from-1e-6-to-0.4e-6s/01-cargo-test.txt

## Review Notes

MIN_BEAM_ON_TIME set to 0.4e-6 with doc comment. tests/beams.rs short_pulse_is_rejected: on_time 5e-7 -> 3e-7 (5e-7 now valid). tests/unit_cell.rs last-pulse short-span test: max_dose 4.5 -> 4.3 and expected span 5e-7 -> 3e-7. Two design/plan docs updated to 0.4e-6. THINNING_FIRST_OFFSET untouched.

## Blockers / Challenges

None

## Commits

- 5b3d9b9 chore(gps): plan 2026-11-13__relax-beam-on-time-constrain-from-1e-6-to-0.4e-6s

## Time Spent

10m (ticket-start.js to ticket-complete.js)
