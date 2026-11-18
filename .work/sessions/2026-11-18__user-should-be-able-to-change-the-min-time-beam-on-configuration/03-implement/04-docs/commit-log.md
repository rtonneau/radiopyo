# Ticket 04: docs

**Status:** ✅ Done

## Local Test Result

`grep -c set_min_beam_on_time README.md CLAUDE.md .work/adr/0008-...md`: 1 / 1 / 1.

## Review Notes

ADR 0008 amendment 2026-11-18 (global setting, default 0.4e-6 s, spans <= T_ZERO analytic, radkin rev 5f7bba8 takes it). Glossary 'Short ON span' names the setting, no fixed value. README beam timing and CLAUDE.md (driver + Rust parity) describe rp.set_min_beam_on_time.

## Blockers / Challenges

None

## Commits

- 28ea86d docs: configurable min beam ON time (README, CLAUDE.md) (ticket 04)

## Time Spent

1m (ticket-start.js to ticket-complete.js)
