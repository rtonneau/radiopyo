# Ticket 02: radiopyo-min-on-time

**Status:** ✅ Done

## Local Test Result

radiobio_env python -m pytest -q: 854 passed (incl. tests/rust with the rebuilt extension, golden tests unchanged); ruff check radiopyo: All checks passed; mypy radiopyo: no issues in 43 source files.

## Review Notes

MIN_BEAM_ON_TIME = 0.4e-6 in beam.py (UnitCell inherits it), rust.py _RADKIN_CONSTANTS, unit_cell docstring, README, CHANGELOG (new Unreleased section), CLAUDE.md, GLOSSARY, ADR 0008 amendment (2026-11-13). Cargo.toml radkin rev 088f06b (Cargo.lock updated). test_min_beam_on_time.py: SHORT 5e-7 -> 3e-7, new tests (0.5 us pulse runs, 0.3 us raises). test_short_on_span.py: short-span lengths moved below 0.4 us (5e-7 -> 3e-7, max_dose 10.5 -> 10.3, micro train 2.05e-5 -> 2.03e-5, clip at 9.7e-6). No golden data re-recorded.

## Blockers / Challenges

maturin develop failed because rust/python/radiopyo_rs/_core.pyd was locked by a running radiobio_env python process; renamed it to _core.pyd.old (git-ignored) and rebuilt. No process was killed. The old .pyd.old can be deleted once that process exits.

## Commits

- cdc4f7e feat: MIN_BEAM_ON_TIME = 0.4e-6 s, radkin pinned at 088f06b (ticket 02)

## Time Spent

5m (ticket-start.js to ticket-complete.js)
