# Implementation Plan

**Session:** config-file-warnings
**Date:** 2026-11-18T16:30:57.511Z
**Estimated effort:** 1 day

## Strategy

Start with the local, low-risk warnings in the SimEnv build (old acid/base syntax, H2O and unknown labels in concentrations), then the config-level removal of disabled reactions (touches the TOML, dict and Rust paths), then the chemistry check, which reads the cleaned env, and finally docs and the changelog.

## Tickets Overview

- **Ticket 1:** Warnings for the old acid/base syntax, H2O in `[concentrations.*]` (no crash) and unknown `[concentrations.initial]` labels.
- **Ticket 2:** Disabled reactions: a zero scalar rate constant removes the reaction (whole pair for acid/base and equilibrium) from the config, with a warning (ADR 0011).
- **Ticket 3:** Chemistry check: `SimEnv.check_chemistry()` / `UnitCell.check_chemistry()` (aliases, charge balance, dead reactants, unused species), auto-run at build, logged once per distinct chemistry.
- **Ticket 4:** Docs (README, docs/configuration.md) and changelog entry.

## Sequencing Rationale

Ticket 1 is independent and small. Ticket 2 must precede ticket 3: the check must not report species that only lived in a disabled reaction, and its fingerprint is computed on the cleaned chemistry. Docs come last, describing the final behaviour.

## Risks & Mitigation

- **Risk:** a bundled config gives chemistry-check findings (e.g. charge read from labels) → **Mitigation:** test asserting zero findings on every bundled config; a fix that changes results stops for approval.
- **Risk:** disabled reactions removed in SimEnv but still in `UnitCell.config` break the Rust parity guard → **Mitigation:** removal is one config-level function used by both `SimEnv.from_dict` and `UnitCell.from_config_dict` (stored config is the cleaned one); idempotent, so no duplicate warnings.
- **Risk:** new warnings change golden outputs or existing tests that assert on logs → **Mitigation:** full pytest including tests/golden after each ticket.

## Assumptions

- Warnings go through `logging.getLogger("radiopyo")`; tests use `caplog`.
- Dev env: `C:/Users/rtonneau/AppData/Local/miniconda3/envs/radiobio_env/python.exe` (written `$PY` below).
