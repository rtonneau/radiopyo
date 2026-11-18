# Ticket 02: disabled-reactions

**Status:** ✅ Done

## Local Test Result

`pytest -W error::RuntimeWarning tests/parser/test_disabled_reactions.py -q` → 16 passed. Full `pytest -q` → 876 passed (59.9 s, radiopyo_rs installed so tests/rust ran, 0 skipped). `ruff check radiopyo` → All checks passed! (after fixing one SIM102). `mypy radiopyo` → Success: no issues found in 43 source files. Manual: bundled pure-water config included, `e_aq + OH_r -> OH_minus` overridden with k_value = 0, `run(backend="rust")` passes the parity guard (one warning, result shape (15, 706)).

## Review Notes

Base: b7d50453d980de574e7815d9596090cb0d179d43 (pre-existing user changes in docs/issues/fenton/ and the other session dir untouched).
Criteria: `disabled_constant` / `drop_disabled_reactions` in parser/utils.py (`_RATE_CONSTANTS` per section; bool and list constants never zero); one warning per removed reaction, "pair" for acid_base/equilibrium. Called in `SimEnv._from_validated_dict` and `UnitCell.from_config_dict` (stored config is the cleaned one); returns the same object when nothing is dropped, so the second call is silent. Input not mutated (test). Include override with k_value = 0 removes the included reaction (TOML test). Species only used by the removed reaction are absent; no RuntimeWarning (test under `simplefilter("error")`). ADR 0011 matches the implementation, unchanged.

## Blockers / Challenges

None.

## Commits

- e33ef0e feat: a zero rate constant disables the reaction (ticket 02)

## Time Spent

3m (ticket-start.js to ticket-complete.js)
