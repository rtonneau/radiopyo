# Ticket 03: packaging-hygiene

**Status:** ✅ Done

## Local Test Result

Ran `python -m pytest -q && python -m build && <wheel IBA check>`: 262 passed; build OK; printed `False`. Sdist has 0 IBA entries. `ruff check radiopyo` shows no deprecation warning (28 pre-existing lint errors remain).

## Review Notes

- Parent: base ed67c93a4897709415b825bc8587359acc17fc3a; pre-existing changes: gps session state only.

All criteria addressed: dev/tqdm extras; .gitignore entries added (parent must `git rm --cached` dist/, .coverage, tests/output/); IBA excluded via tool.hatch.build exclude; ruff moved to [tool.ruff.lint]; __all__ deduplicated, MacroPulsedBeam exported, logging/beam __all__ extended; logging config now only from default_log_config.toml (always loaded, disable_existing_loggers=false, InfoOnly filter removed, [tool.logging] dropped from pyproject); URL set to GitHub; tests/test_logging.py added (subprocess). INFO output format verified unchanged.

- Parent acceptance: parent reviewed the full diff inline (packaging + logging, localized, no numeric code); reviewer skipped. Parent untracked dist/, .coverage, tests/output/ and committed by hand. Model haiku (requested fast).

## Blockers / Challenges

None

## Commits

- 9a8e758 chore: packaging hygiene, logging without side effects (ticket 03)

## Time Spent

2m (ticket-start.js to ticket-complete.js)
