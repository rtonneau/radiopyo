# Ticket 09: chemistry-cleanup

**Status:** ✅ Done

## Local Test Result

```
$ python -m pytest -q && ruff check radiopyo && mypy radiopyo
410 passed in 29.98s
All checks passed!
Success: no issues found in 39 source files
```

- Baseline at 962c6b4: 391 passed, ruff clean, mypy clean. The 19 new tests account for the difference.
- Golden tests: `PYTHONHASHSEED=0` and `PYTHONHASHSEED=12345` both give `tests/golden` 5 passed. The golden files, tolerances and tests are unchanged.
- Evidence (in `.scratch/tests/2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian/evidence/`): `09-verify.txt`, `09-golden-seeds.txt`, `09-snapshot-before.json`, `09-snapshot-after.json` and `09-snapshot-diff.txt`.

## Review Notes

- Parent: base 962c6b443c7e04cb4f13169694a1bdd2181e6c63; pre-existing changes: gps session state only.

- AC1 (one helper): `BasicReaction._add_species_from_generic(generic, list_sp, reactants=, products=)` adds reactants, then products, in the old order. KReaction, Arrhenius, CstReaction (products only) and MichaelisMenten all use it in place of their 4 loops. MichaelisMenten still rejects more than 1 substrate with the same ValueError, now checked before the helper runs.
- AC2 (BiReaction owns its k):
  - `BiReaction._add_reaction_pair` sets the k of the forward KReaction and of the reverted backward KReaction. `_k_forward`/`_k_backward` are gone from AcidBase and EReaction.
  - `BiReaction.k_forward`/`k_backward` are now concrete and return `forward.k_value(T)`/`backward.k_value(T)`. This also merges the identical H2O branches.
  - `AcidBase.get_couple` uses one `_couple_member` helper for both the acid and base branches; the error texts are unchanged. `AcidBase.to_dict` reads `KReaction.k_polynomial`.
  - The snapshot shows AcidBase and EReaction kf/kb/keq/to_dict/get_couple are bit-identical.
- AC3 (collection lookups):
  - `KReactionCollection._reindex()` now runs after `sort_by_label` and `sort_by_reaction_index`.
  - `update_label` now renames **in place**, so the reaction keeps its position. Before, it appended at the end and left the indexer stale. No caller relied on the old position.
  - Tests are in `tests/chemistry/test_reaction_collections.py`.
- AC4 (species):
  - An empty `SimSpecies.cc_init` is now `None` instead of the `Polynomial(-inf)` sentinel; `initial_cc()` still returns 0.
  - `BasicSpecies.__eq__` returns NotImplemented for foreign types, and `__hash__ = hash(label)` was added.
  - `ListSpecies.sort`/`set_order` set `next_index` from lengths, so empty lists no longer hit an unbound variable.
  - `ListSpecies.copy` copies each SimSpecies (index and cc kept) and keeps `next_index`/`sorted`.
  - Outside the listed files, `simulation/environment_build.py` needed 3 None guards for the Optional `cc_init` (the H2O constant and two `len(cc_init) > 1` checks in `convert_to_config_dict`).
- AC5 (shared make_hash/as_label):
  - `BasicReaction.make_hash` is concrete and digests `_hash_source()`, which defaults to `repr(self)`. The 3 copies were removed.
  - `BasicReaction._index_label(prefix)` backs the KR_/CR_/EZ_ labels. The EZ error text is now the shared KeyError text.
  - MichaelisMenten `__repr__` now has its closing `)`. Its `_hash_source` strips it, so `make_hash()` is unchanged (checked in the snapshot and by a test).
  - Decision for `RadiolyticReaction.label`: it is now a property, consistent with `species.label`. It returns a `_CallableLabel` str subclass, so the old `reaction.label()` still works and returns a plain str (used in unit_cell.py:375 and test_radiolytic). I chose this as the least surprising option: both call forms keep working and str/repr stay `RadiolyticReaction('e_aq')`.
- AC6 (docs): the RawSpecies/ChemSpecies/BaseChemicalReaction/base_chemistry mentions are gone from radiopyo/chemistry. The KReaction docstring documents the log10(k) storage and the n·k convention (ADR 0001, scalar k only).
- Snapshot at base vs after (`09-snapshot-diff.txt`), checked on all 3 bundled configs:
  - Identical: every `str(reaction)`, every `as_label()`, all make_hash values except the AcidBase sub-reactions, `convert_to_config_dict(True)`, and all bireaction values.
  - Expected differences:
    - The repr/make_hash/to_dict of the 6 AcidBase sub-KReactions now show their k instead of `k_value=None` (AC2).
    - The MichaelisMenten repr gains its `)`; its hash is unchanged.
    - The SimSpecies repr shows `None` for an empty concentration instead of `inf` (AC4).
  - No stored key depends on these: the only consumer of reaction/species repr is `SimEnv.make_hash`, which has no callers and already raises `TypeError: unhashable type: 'dict'` at base and after.
- Tests adapted (they tested removed internals or sentinels):
  - `test_species.py::test_simspecies_init` now expects `cc_init is None` instead of `Polynomial(-inf)`.
  - `test_raw_species_eq_unsupported_type` now expects `False`/NotImplemented instead of a raised NotImplementedError.
  - Two fixture docstrings now say SimSpecies.
- Tests added:
  - `test_raw_species_hash`.
  - In `test_species_collection.py`: the empty sort/set_order tests and `test_copy_keeps_indices`.
  - `tests/chemistry/test_reaction_collections.py` (15 tests).
- Not touched, as the ticket says: Arrhenius units/temperature, enzymatic raw_k n·k (17), temperature propagation (13), set-iteration species order (16), include merge (14). `radiopyo/chemistry/utils/factories.py` needed no change.

- Parent acceptance: independent reviewer (opus/high): no findings, old-vs-new bundled config builds/matrices/config dict identical; parent re-ran pytest (410 passed), ruff 0, mypy 0. Intended diffs: acid/base sub-reactions now carry k (repr/make_hash changed, no consumer), MM repr closes ')' with make_hash unchanged, RadiolyticReaction.label is now a callable str property. 2 tests adapted (removed internals).

## Blockers / Challenges

- None blocking. A bash heredoc with a large Python edit failed on quoting, so the source edits were applied with small Python scripts in the scratch dir (`edit_09_*.py`).
- The enzymatic grammar rejects 2 substrates, so the MichaelisMenten single-substrate test builds its generic by hand.

## Commits

- a0bee93 refactor: chemistry cleanup, shared reaction population and BiReaction owning KReactions (ticket 09)

## Time Spent

20m (ticket-start.js to ticket-complete.js)
