# Ticket 14: fix-include-merge

**Status:** ✅ Done

## Local Test Result

Verification Step: `python -m pytest -q && ruff check radiopyo && mypy radiopyo`

```
501 passed in 45.25s
All checks passed!
Success: no issues found in 42 source files
```

- `tests/parser/test_includes.py`: 40 tests (was 10), all pass; includes the former `TODO(ticket 14)` first-wins test turned into later-wins.
- `tests/golden/` not touched and all golden tests pass unchanged, including the HDF5 one (see Blockers / Challenges).
- Before/after report: `.scratch/tests/2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian/evidence/14-before-after.md`.

## Review Notes

- Parent: base f8c1777ecd4d71ab65151036f83bb51370900fc3; pre-existing changes: gps session state only.

Checked against the Acceptance Criteria:

- One generic `combine_config` keyed on `reaction_key` (frozensets of (label, stoi) for reactants, (label, stoi, operation) for products, enzyme label): order-insensitive, arrow direction kept. It only calls `RPARSER.parse_reaction`, builds no reaction object, and `radiopyo/parser/*.py` has no `radiopyo.simulation` import (grep, and an AST test over the parser sources: `test_parser_package_does_not_import_the_simulation_layer`; `sys.modules` cannot prove it because `radiopyo/__init__` imports UnitCell).
- Enzymatic override keeps k (`5e7` stays `5e7`, was `1e8`); acid/base override with only `k_forward` keeps `k_backward` and `pKa` (and builds a valid `SimEnv`); `constant` override replaces (was duplicated, it was compared with the k-reactions); `B + A -> C` overrides `A + B -> C`; arrhenius merges field by field too. An `equilibrium` override giving 2 or more of k_forward / k_backward / k_equilibrium replaces all three (otherwise a stale `k_backward` wins over the new `k_forward` in `eq_reaction_from_dict`).
- Precedence: later include overrides earlier ones (`[a, b]` gives b), the including file overrides all of them, nested includes work; same rule for parameters, fixed/initial concentrations, radiolytic and reaction lists. Overridden reactions keep their position, new ones are appended.
- Paths: relative entries are anchored to the including file by `FileParser.anchor_includes` when it is parsed (nested relative includes work from any CWD), falling back to the CWD (`include_path`); only the exact `radiopyo/` prefix maps to `PACKAGE_DATA_PATH`, so `radiopyo_mine.toml` is a plain path.
- Bundled configs: `config_Labarbe.toml` resolves to the same entries as before in every section (set equality, checked by `ba14.py` and by `test_bundled_config_is_the_union_of_its_includes`); only the ORDER of the k-reactions changes (includes in list order: advanced_water_radiolysis, then biology; was biology first), so global reaction indices of configs with includes change. A run on Labarbe differs by at most 7.9e-10 relative on species above 1e-6 uM (trace species such as HO3 at 4e-30 uM move by 7%, summation order); goldens are keyed by label and pass unchanged.
- `combine_config` / `resolve_includes` no longer mutate their inputs (tests `test_combine_config_does_not_modify_its_arguments`, `test_resolve_includes_does_not_modify_the_given_config`).
- Docs: README `[includes]` section rewritten (paths, precedence, merge rules); CHANGELOG entry appended under `[Unreleased]` / `### Fixed` with the behaviour changes flagged (last include wins, path resolution, reaction order).
- Scope kept: no beams (15), species order (16) or Arrhenius/enzymes add (17) change; nothing under `tests/golden/` edited; nothing staged or committed.

Changed after review of my own work:

- Equilibrium rule added after reading `eq_reaction_from_dict` (k_backward is preferred over k_forward when both are present).
- `blank_config()` / `_complete()` introduced so that `parameters` has no default value at parse time and a merge never mutates or aliases its inputs.

- Parent acceptance: implementer sonnet; reviewer opus/high found one P3 (circular includes recursed to RecursionError), fixed by parent with a _chain guard plus test. Parent re-ran full suite (502 passed), ruff 0, mypy 0. Golden unchanged; k-reaction order follows include order (float-order diffs only, evidence/14-before-after.md). USER APPROVED; catalase HDF5 doubling traced to MichaelisMenten.to_dict raw_k_value, to fix in ticket 17 (hdf_expected only).

## Blockers / Challenges

- HDF5 golden did NOT flip (the brief expected it to). The catalase doubling pinned in `tests/golden` (1.324e8 -> 2.648e8, `compute_hdf` docstring) is not produced by the include merge: `MichaelisMenten.to_dict()` writes `raw_k_value()` (k x stoichiometry 2) while `MichaelisMenten.from_generic` takes the file value as is, and the H5 json store goes `ODEResult.to_json_config` -> `from_json_config` -> `SimEnv._chemistry_builder`, which never calls `resolve_includes` / `combine_config`. Measured with the base tree and the working tree: in memory 1.324e8, written 2.648e8, read back 2.648e8, both before and after. No golden array differs, so no number to approve. Fixing the round trip means changing `radiopyo/chemistry/michaelis.py` (not in this ticket's files); the parent should decide whether to open a ticket (it would change the pinned rebuilt catalase rates by a factor 0.5 in that golden).
- The old merge only doubled k when an enzymatic reaction was OVERRIDDEN (`raw_k_value()` written back into the config); that is what this ticket fixes. Several "before" cases in the evidence file raise instead of returning wrong numbers (enzymatic override without `k_micha`, acid/base override without a full set of constants): the old code required complete overrides.
- "TOMLFileParser overwriting default parameters": `parse()` already replaced the default dict by the file content, so there was no observable bug on this path (before and after give the same resolved parameters, shown in the evidence file). I made it explicit (`blank_config()`, no default parameter at parse time) and pinned it with a test. The `toml.py:271-272` reference of the brief does not match the current 68-line file.
- Tooling: heredoc with Python code containing quotes broke once in the shell; the test file was written with the file tool instead. Line endings kept (CRLF in the edited sources, README, CHANGELOG and this log; LF in `tests/parser/test_includes.py`).

## Commits

- 662277a fix: generic include merge, later includes win, relative include paths (ticket 14)

## Time Spent

54m (ticket-start.js to ticket-complete.js)
