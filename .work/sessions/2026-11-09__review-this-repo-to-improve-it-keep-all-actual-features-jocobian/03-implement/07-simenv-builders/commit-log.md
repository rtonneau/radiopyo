# Ticket 07: simenv-builders

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo` -> 384 passed (316 baseline + 68 new in tests/simulation/test_env_queries.py), ruff "All checks passed!", mypy "no issues found in 36 source files".
`tests/golden` + new tests also pass under PYTHONHASHSEED=0/1/7 (73 passed each). Evidence: .scratch/tests/<session>/evidence/07-verify.txt, 07-golden-seeds.txt.
Before/after snapshot (snap_queries.py / snap_index.py, evidence/queries_before.json vs queries_after.json, index_before.json vs index_after.json): the 10 old query names for every species of the Labarbe and pure-water configs, and the (global index, reaction) list, are byte-identical.

## Review Notes

- Parent: base 06f4dbea0b1f4f0bd06d1bfa587da2e6801bac37; pre-existing changes: gps session state only.

- AC1 (table-driven reader): `_REACTION_BUILDERS` table + `reactions_from_dict(section, ...)` replace the six copy-pasted loops; `acid_base_from_dict` & co. stay as thin wrappers (used by parser/utils.py); `_chemistry_builder` reads sections in the same order, so species creation order and `ReactionsList` index order (ab, eq, k, enzymatic, radiolytic) are unchanged (test_global_index_order + snapshot). `radiolytic_from_dict` is not reaction-string based and was left as is. Minor: `michaelis_from_dict({})` now returns [] (was KeyError).
- AC2: `reactions_involving(species, role="any"|"reactant"|"product", kinds="k_all"|"k"|"eq"|"enzymatic"|"radiolytic")` (typing overloads). The default (role any, k_all) still returns (reaction, k) tuples, so the historical `reactions_involving(species)` is unchanged; the other 9 historical names are thin aliases. Lists of species mean "all of them" for every kind. Species given as SimSpecies now work for role="any" (labels are passed to `in`, which used to raise NotImplementedError for SimSpecies). test_env_queries.py re-implements the 10 old queries verbatim and compares identity + labels on both bundled configs.
- AC3: removed `SimEnv.pKw`, `pKH2O`, `update_temperature`, `update_flag_positive` (grep of radiopyo/, tests/, examples/*.ipynb, README: no user; the `temperature`/`force_positive` keys of `parameters` are untouched, ticket 12/13). `update_flag_GValues` kept (read by UnitCell). `pH or ...` in `update_flag_pH` replaced by an `is None` test, so pH=0 is honoured; no golden touches it (goldens pass).
- AC4: `get_radiolytic_reaction_from` on an empty list now raises the intended ValueError (was UnboundLocalError); `from_dict` raises ValueError naming the missing `concentrations`/`parameters` section (was KeyError after a misleading warning); log level restored in a `finally` (from_dict split into from_dict + _from_validated_dict).
- AC5: golden tests pass unchanged.
- Noticed, not fixed (out of scope): `SimEnv.make_hash` raises TypeError (`frozenset(parameters.items())` with the unhashable `species_order` dict); the method has no caller. Candidate for a later ticket.
- Process note: I ran `ruff format radiopyo` once by mistake, which reformatted 27 files; I restored all of them except environment_build.py with `git checkout -- <file>` (they had no other uncommitted work). environment_build.py was formatted with ruff (it had 2-3 pre-existing format deviations); the working copy is CRLF as before.

- Parent acceptance: independent reviewer (sonnet/high): no findings; parent re-ran pytest (384 passed), ruff 0, mypy 0. pH=0 now honoured (no golden affected). Seed for ticket 18: SimEnv.make_hash always raises TypeError (unhashable species_order), unused.

## Blockers / Challenges

None blocking. One test-design finding: `species in reaction` raises NotImplementedError for SimSpecies, so the SimSpecies comparison test uses the label for the legacy reference.

## Commits

- 3f8035f refactor: table-driven SimEnv readers and reactions_involving query (ticket 07)

## Time Spent

11m (ticket-start.js to ticket-complete.js)
