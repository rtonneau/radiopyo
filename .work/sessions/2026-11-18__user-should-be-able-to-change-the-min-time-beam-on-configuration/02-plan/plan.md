# Implementation Plan

**Session:** User should be able to change the min time BEAM ON configuration
**Date:** 2026-11-18T16:48:06.628Z
**Estimated effort:** 1 day

## Strategy

Python side first (global setting, beam check, run-time short span path, interim Rust guard), so radiopyo is usable on its own. Then radkin takes the threshold as an input on its current main (pushed only after the user's OK), then radiopyo re-pins radkin and sends the value, then docs.

## Tickets Overview

- **Ticket 1:** `rp.set_min_beam_on_time` / `rp.min_beam_on_time`, read by the beam check and the UnitCell run-time short span path (floor T_ZERO, warning below it); Rust backend refuses a non-default value until ticket 3.
- **Ticket 2:** radkin (C:\DEV\RUST\radkin, on main): the minimum beam ON time becomes an input (default 0.4e-6 s) of the beam check and the short span check; pushed after the user's OK.
- **Ticket 3:** radiopyo pins the new radkin rev, the binding and `radiopyo/rust.py` send the current value, parity covers it; Rust tests with a non-default value.
- **Ticket 4:** docs: ADR 0008 amendment, glossary, README, CLAUDE.md.

## Sequencing Rationale

Ticket 1 has no Rust dependency and fixes the Python semantics the radkin change must mirror. Ticket 3 needs ticket 2's pushed rev. Docs last, describing the final behavior.

## Risks & Mitigation

- **Risk:** re-pinning radkin brings 145 commits (IRT stage, 36 touching src/) → **Mitigation:** full `tests/rust/` parity and golden tests must pass with the new rev before ticket 3 is done; any drift stops the ticket.
- **Risk:** global state leaks between tests → **Mitigation:** an autouse fixture restores the default in every test that changes it.
- **Risk:** a span between T_ZERO and the setting, or below T_ZERO, reaching the solver → **Mitigation:** the short span test uses max(setting, T_ZERO); tests with a value below T_ZERO and a span shorter than T_ZERO.

## Assumptions

- The radiobio_env (C:\Users\rtonneau\AppData\Local\miniconda3\envs\radiobio_env) has a Rust toolchain and maturin for ticket 3; radkin pushes go to GitHub `rtonneau/radkin` only after explicit user OK.
