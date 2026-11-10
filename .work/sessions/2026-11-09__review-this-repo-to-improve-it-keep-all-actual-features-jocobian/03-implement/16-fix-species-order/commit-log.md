# Ticket 16: fix-species-order

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`: 530 passed (529 + 1 new), ruff all checks passed, mypy no issues in 42 files. Goldens untouched and passing.

## Review Notes

- Parent: base 9d296b69f5174d9d00ba7e065b9902e3ced13160; pre-existing changes: gps session state only.

- AC1: cause was `ListSpecies.sort` numbering `all_dyn - set(starter)` (a set). Now iterates `out.dyn_sp` (insertion order = config-file order, chosen over label sorting because it keeps the config's meaning and needs no extra rule); the starter set is used for membership only. Other set/frozenset uses checked: `parser/utils.reaction_key` (order-insensitive key, never iterated for indices) and `SimEnv.make_hash` (repr of a frozenset of items: hash-seed dependent string, not an index; left unchanged, out of ticket scope, worth a follow-up). Reaction global-index order is untouched. New test `test_species_and_reaction_order_independent_of_hash_seed` (subprocess, seeds 0/1/12345, Labarbe config, species in index order + reactions via `iter_all`) fails on the old code and passes on the new.
- AC2: report at `.scratch/tests/2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian/evidence/16-before-after.md`. Before (HEAD, seed 1): 2441 points; after seeds 1/7/99: 2432 each. Max rel. diff before vs after (entries > 1e-6 uM): 1.3e-7 on R_r, others <= 2.8e-9; after seeds identical (0.0). CHANGELOG `### Fixed` entry added.

- Parent acceptance: implementer haiku/medium; reviewer skipped (one function, covered by a subprocess seed test). Parent re-read the diff, re-ran full suite (530 passed), ruff 0, mypy 0. USER APPROVED before/after (evidence/16-before-after.md): Labarbe run 2432 points at seeds 1/7/99, max rel diff vs old code 1.3e-7. Seed for 18: SimEnv.make_hash hashes repr(frozenset).

## Blockers / Challenges

None (first draft of the test snippet had a newline-escape bug and the logger prints to stdout; fixed with a marker line and chr(10)).

## Commits

- 5b3b9e2 fix: deterministic species order, independent of PYTHONHASHSEED (ticket 16)

## Time Spent

9m (ticket-start.js to ticket-complete.js)
