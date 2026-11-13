# Implementation Plan

**Session:** Relax beam ON time constrain from 1e-6 to 0.4e-6s
**Date:** 2026-11-13T14:43:36.409Z
**Estimated effort:** 2 hours

## Strategy

Change radkin first (it is pinned by rev), push it, then change radiopyo, bump the rev and rebuild the extension.

## Tickets Overview

- **Ticket 1:** radkin: MIN_BEAM_ON_TIME = 0.4e-6, tests updated, committed and pushed.
- **Ticket 2:** radiopyo: constant, docs, tests, ADR 0008 amendment, CHANGELOG, `_RADKIN_CONSTANTS` and the radkin rev bump in `rust/Cargo.toml`.

## Sequencing Rationale

radiopyo's Cargo.toml pins a radkin rev, so radkin must be pushed first.

## Risks & Mitigation

- **Risk:** a golden or bundled case has an ON time between 0.4 and 1 µs and its result changes → **Mitigation:** run the full pytest; report it instead of re-recording.

## Assumptions

- The user has approved pushing radkin to GitHub (asked before ship).
