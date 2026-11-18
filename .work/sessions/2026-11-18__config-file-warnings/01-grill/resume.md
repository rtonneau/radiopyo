# Session: config-file-warnings

**Date:** 2026-11-18T16:19:49.648Z
**Status:** Grill phase complete

## Problem Statement

Config-file mistakes are either silent or end in an unhelpful crash: the old acid/base syntax fails with a generic message, a zero rate constant silently builds a reaction with log10(k) = -inf, H2O in `[concentrations.*]` crashes (fixed) or is silently ignored (initial), and mixed spellings of one species (OH / OH_r, O_minus / O_r_minus) become separate species. Users also have no way to check the chemistry they defined.

## Context & Constraints

- **Current behavior:** (probed 2026-11-18) old `acid`/`base` keys → `AttributeError: Need a reaction string AND a pKa…`; `k_value = 0` → accepted, numpy `RuntimeWarning: divide by zero in log10`; `H2O` in `[concentrations.fixed]` → `KeyError: 'H2O is already recorded'`; `H2O` or any unknown label in `[concentrations.initial]` → silently ignored (`suppress(KeyError)` in `SimEnv` build); `OH` + `OH_r` → two distinct dynamic species.
- **Pain point:** an `[includes]` override cannot delete a reaction, so users write `k_value = 0` to disable one; alias mixing (seen in `docs/issues/fenton/basic_water_radiolysis.toml`, `O_r_minus` vs bundled `O_minus`) silently splits a species.
- **Dependencies:** Rust backend rebuilds from the stored config (ADR 0010) and the parity guard compares labels: disabled reactions must be removed at config level. Global reaction index order (acid/base, equilibrium, k, enzymatic, radiolytic) unchanged. Golden data unchanged (no bundled config has a zero k). Warnings through `logging.getLogger("radiopyo")` (no `print`, ruff T20). ruff and mypy stay at zero.
- **Tech stack:** Python 3.11, `radiopyo/parser/` (TOML, include merge `combine_config`), `radiopyo/simulation/environment_build.py` (`SimEnv.from_dict`, `_REACTION_BUILDERS`), `radiopyo/simulation/unit_cell.py`, pytest.

## Success Metrics

- An old-syntax acid_base entry logs a warning showing the converted entry (`reaction = "<acid> <-> H3O_plus + <base>"`, k_forward/k_backward required) and then raises the same error type as today, whose message points to the warning.
- A reaction with a zero scalar rate constant (k_value of k_reaction/constant/enzymatic, Arrhenius A, k_forward/k_backward/k_equilibrium of a pair; one zero side drops the whole pair) is absent from the built env and from `UnitCell.config`, with one warning per reaction and no numpy warning; this holds for TOML files (after include merge) and for `SimEnv.from_dict` on a dict.
- `H2O` in `[concentrations.fixed]` or `[concentrations.initial]` logs a warning that it cannot be replaced, keeps the built-in H2O and does not crash; an unknown label in `[concentrations.initial]` logs a warning.
- `SimEnv.check_chemistry()` and `UnitCell.check_chemistry()` return findings for: species aliases, unbalanced charge per reaction, dead reactants, unused species; the automatic run at build logs each distinct chemistry's findings once per process; a manual call always reports.
- All bundled configs (`advanced_water_radiolysis.toml`, `reactions_biology.toml`, `config_Labarbe.toml`) give zero findings; full pytest, ruff and mypy pass.

## Architecture & Approach

1. **Old acid/base syntax**: in `_build_acid_base`, when `reaction` is missing and `acid`/`base` are present, log the warning with the converted entry, then raise as today.
2. **Disabled reaction** (ADR 0011): one function (config level, e.g. `drop_disabled_reactions(config)`) removes every reaction with a zero scalar rate constant, logging one warning each; it runs on the merged config so both the TOML path and `SimEnv.from_dict` go through it, and `UnitCell.config` (sent to Rust) no longer holds them. Polynomial k are never zero (log10 coefficients) and are untouched.
3. **Concentrations**: in `SimEnv` build, `H2O` in fixed/initial is skipped with a warning before `add_cst_species` (no more KeyError); an initial label not found among species logs a warning instead of the silent `suppress(KeyError)`.
4. **Chemistry check**: `SimEnv.check_chemistry() -> list[ChemistryFinding]` (frozen dataclass: kind, species, reactions, message), delegated by `UnitCell.check_chemistry()`. Checks: species aliases (labels equal after removing `_r` tokens, plus a small table: e_aq/e_minus, H_plus/H3O_plus); charge balance per k/enzymatic/pair reaction from `_plus`/`_minus` tokens (e_aq = -1; negative products counted on the consumed side; enzymes, `trash` and radiolytic reactions skipped); dead reactants (dynamic species with no initial concentration, no radiolytic yield, never produced, but consumed); unused species (in `[concentrations.*]` or radiolytic yields, in no reaction). Run at the end of `SimEnv.from_dict`; a module-level set of chemistry fingerprints (reactions, species, concentrations given) makes the automatic run log only for a chemistry not yet seen in the process.
5. **Docs**: README / `docs/configuration.md` (disable by k = 0, H2O note, `check_chemistry`), CHANGELOG Unreleased entry.

## Assumptions & Trade-offs

- k = 0 means disabled (implicit, no `enabled` field): chosen by the user; existing results do not change since a zero-k reaction contributed nothing, but reaction indices and species lists shrink accordingly.
- One zero side of a pair drops the whole pair (keeps the 2-reactions-per-pair index rule).
- The `_r` rule relies on same formula + charge ⇒ same species; the alias table covers only known non-`_r` cases. Formula reordering (HO2 vs O2H) is not detected.
- Charge balance reads charge from labels only; labels without `_plus`/`_minus` tokens are neutral except table entries (e_aq).
- Not doing: auto-conversion of the old acid/base syntax, mass balance, an option to switch the chemistry check off (log-once covers repeated builds).

## Open Questions

- If a bundled config produces a finding whose fix would change simulation results, stop and ask before changing it (result-changing fix).

## Notes

- Glossary: Disabled reaction, Species alias, Chemistry check ("sanity check" is avoided: it names the parity guard).
- ADR 0011: a zero rate constant disables the reaction.
- Probe script: `.scratch/tests/2026-11-18__config-file-warnings/probe.py`.
