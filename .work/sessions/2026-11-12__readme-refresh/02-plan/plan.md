# Implementation Plan

**Session:** readme-refresh
**Date:** 2026-11-12T08:47:18.301Z
**Estimated effort:** 2 hours

## Strategy

Move the detailed reference out first (nothing is lost), add the logo and CI, then rewrite the README that links to all of them.

## Tickets Overview

- **Ticket 1:** Move the detailed config/beam/run reference to docs/configuration.md and add docs/logo.svg.
- **Ticket 2:** Add a GitHub Actions workflow running ruff, mypy and pytest.
- **Ticket 3:** Rewrite README.md (< 150 lines) with logo, badges, install/update/usage, condensed config, citation.

## Sequencing Rationale

The README links to docs/configuration.md, docs/logo.svg and the workflow badge, so those exist first.

## Risks & Mitigation

- **Risk:** CI fails on GitHub (golden tests are hash-seed sensitive) → **Mitigation:** set PYTHONHASHSEED=0 and run the same commands locally.
- **Risk:** reference content lost in the move → **Mitigation:** verbatim move, heading count check.

## Assumptions

- No PyPI release; install from GitHub only. Version badge static at 0.6.0.
