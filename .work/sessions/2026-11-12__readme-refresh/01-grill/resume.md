# Session: readme-refresh

**Date:** 2026-11-12T08:42:02.935Z
**Status:** Grill phase complete

## Problem Statement

The README is 485 lines, mostly a detailed TOML/beam reference. A new user cannot quickly find how to install, update and use radiopyo, and nothing asks users to cite the paper.

## Context & Constraints

- **Current behavior:** README mixes a pitch, developer notes, install, a full config reference, run modes and acknowledgements; a few static badges; no logo; no citation; no CI.
- **Pain point:** too long, not straight to the point; no citation request (Tonneau et al 2025 Phys. Med. Biol. 70 235021, DOI 10.1088/1361-6560/ae1fcc).
- **Dependencies:** keep all reference content (moved, not deleted); keep install from GitHub (`rtonneau/radiopyo`), Python >= 3.11, LGPL-3.0, version 0.6.0.
- **Tech stack:** Markdown, shields.io badges, SVG, GitHub Actions.

## Success Metrics

- README under ~150 lines with logo, badges, pitch, Install, Update, Usage, condensed config table, dev section, Citation (link + BibTeX), Acknowledgement.
- docs/configuration.md holds the full TOML/beam/run reference moved verbatim, linked from the README.
- docs/logo.svg (molecule + radiation wave + "radiopyo" wordmark) renders.
- .github/workflows/tests.yml runs ruff, mypy and pytest on Python 3.11 / ubuntu; README has its badge.

## Architecture & Approach

New docs/logo.svg; README rewritten; detailed sections (Configuring Simulation, beam definitions, Ways of Running Simulation) moved to docs/configuration.md; badges: Tests CI, version, Python, license, Ruff, mypy, DOI; update = `pip install --upgrade git+https://github.com/rtonneau/radiopyo.git` or `git pull` for editable installs.

## Assumptions & Trade-offs

Not publishing to PyPI. CI badge only goes green once the workflow runs on GitHub. Version badge stays static (manual bump on release).

## Open Questions

N/A

## Notes

User chose: condensed config in README + detailed doc in docs/, pip git upgrade, static + CI badges, molecule + wordmark logo.
