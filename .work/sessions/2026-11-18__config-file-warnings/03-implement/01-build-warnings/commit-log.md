# Ticket 01: build-warnings

**Status:** ✅ Done

## Local Test Result

`pytest tests/simulation/test_config_warnings.py -q` → 6 passed. Full `pytest -q` → 860 passed (51.5 s). `ruff check radiopyo` → All checks passed!. `mypy radiopyo` → Success: no issues found in 43 source files. Loading `PURE_WATER_CONFIG_FILE` and `LABARBE_CONFIG_FILE` logs no new warning.

## Review Notes

Base: a9d6c7ed1afe1aa56b93380292bb8b3e31bebe04. Pre-existing user changes (not owned by this ticket, untouched): docs/issues/fenton/* (modified and untracked), .work/sessions/2026-11-18__user-should-be-able-to-change-the-min-time-beam-on-configuration/.
Inline review against criteria: old syntax → `_warn_old_acid_base_syntax` logs the converted entry (reaction string, pKa, k_forward/k_backward required), then `AttributeError` "Old acid/base syntax…see the warning above"; an entry without reaction nor acid/base keeps the former message. H2O in fixed/initial → `_without_water` warns and drops it before the species are built; the built-in H2O stays the fixed value (cst_sp equal to the reference). Unknown initial label → explicit warning instead of `suppress(KeyError)`; a fixed species listed in initial stays silent.

## Blockers / Challenges

Heredoc escaping turned `
` into real newlines in the warning string (SyntaxError); fixed with Edit.

## Commits

- aa6698d feat: warn on old acid/base syntax, H2O and unknown labels in concentrations (ticket 01)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
