# Ticket 03: chemistry-check

**Status:** ✅ Done

## Local Test Result

```
$PY -m pytest tests/simulation/test_chemistry_check.py -q   -> 30 passed
$PY -m pytest -q                                            -> 906 passed in 52.85s
$PY -m ruff check radiopyo                                  -> All checks passed!
$PY -m mypy radiopyo                                        -> Success: no issues found in 44 source files
```
Bundled configs (`.scratch/tests/2026-11-18__config-file-warnings/check_bundled.py`):
`advanced_water_radiolysis.toml` 0 finding; `config_Labarbe.toml` 2 `charge` findings (kept, see Review Notes).

## Review Notes

- New module `radiopyo/simulation/chemistry_check.py` (banner layout): frozen `ChemistryFinding(kind, species, reactions, message)`, `check_chemistry(env, log=True)`, `auto_check_chemistry`, `chemistry_fingerprint`, `label_charge`, `electron_count`, `ALIAS_TABLE`. `SimEnv.check_chemistry()` and `UnitCell.check_chemistry()` call it; a manual call always logs.
- Automatic run at the end of `SimEnv.from_dict`, logged once per fingerprint (reactions in index order, dynamic/fixed species, species with an initial value; rate constants and temperature excluded), so a modified chemistry is checked again. A build with `log=False` neither logs nor marks the chemistry as seen.
- Aliases: `ALIAS_TABLE` (e_minus -> e_aq, H_plus -> H3O_plus) and the `_r` rule, refined (my implementation choice): a label without `_r` joins its `_r` sibling only if its formula has an odd electron count (OH, O_minus, O2_minus). The plain rule flagged RH / RH_r in config_Labarbe, a real molecule and its radical; labels that are not formulas (RH, GSH) are never aliases.
- Charge: per k-reaction, pair (once) and enzymatic reaction; negative products counted on the consumed side, `trash` and enzymes skipped, e_aq = -1.
- Dead reactants: fixpoint over both directions of each pair; a species is available if fixed, initial != 0, radiolytic, or produced by a reaction whose needs (reactants, negative products, enzyme) are available. Catches chains (O2 missing -> O3_minus dead too).
- Unused species: given in fixed/initial/radiolytic but in no reaction (H2O, H3O_plus, OH_minus, trash exempt).
- User decision (asked during this ticket): the 2 charge findings of the bundled `reactions_biology.toml` (included by config_Labarbe), `2O2_minus -> H2O2 + O2` (protons left out) and `RH + e_aq -> RH_r` (RH_r is the radical anion), are **kept as warnings**: no config change, so config_Labarbe logs 2 charge warnings once per process. Test `test_labarbe_known_charge_simplifications` pins them.

## Blockers / Challenges

- Bundled config_Labarbe is not finding-free (see Review Notes): plan criterion "zero findings" amended by the user's answer.
- `.work/sessions/.current-session` was switched to `2026-11-18__user-should-be-able-to-change-the-min-time-beam-on-configuration` at 17:40:20 during this ticket, not by this session (another process?); reported to the user before completing.
- ruff S105 false positive on a `token` variable: renamed.

## Commits

- e458755 feat: chemistry check of the defined chemistry (ticket 03)

## Time Spent

2h 18m (ticket-start.js to ticket-complete.js)
