# Ticket 05: poly-k-warning

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q -p no:cacheprovider && ruff check radiopyo && mypy radiopyo`: 626 passed; ruff: All checks passed; mypy: no issues in 42 files.

## Review Notes

Added `KReaction._warn_nk_not_divided` (logger `radiopyo`, one warning when single reactant has stoi > 1 and k is not scalar). Called from `_apply_nk_convention` (polynomial k), `Arrhenius.from_generic`, and `BiReaction._add_reaction_pair` (shared path for EReaction/AcidBase, forward and backward). Scalar k: no warning. Bundled configs (advanced_water_radiolysis, config_Labarbe incl. reactions_biology) build with no such warning (tested). Tests in tests/chemistry/test_nk_warning.py. Note: equilibrium path is not directly unit-tested; it goes through the shared helper. 0 review rounds.

- Parent: base c6c4c8ab803dfb3c0afce03b59c0fa2f3cde9c01; pre-existing changes: only gps session state, owned by the parent.
- Parent review (inline, the diff is small: about 30 lines in k_reactions.py plus the test file): the warning sits at the three places where the n·k convention applies (`_apply_nk_convention` when k is not a scalar, `Arrhenius.from_generic`, `BiReaction._add_reaction_pair`); there is no double warning, because the pair's sub-reactions are built with k None. One gap: the equilibrium path had no direct test. The parent added `test_equilibrium_warns_on_forward_only` ("2 A <-> B" with an EReaction: one warning, from the forward reaction). Verification: 627 passed, ruff clean, mypy clean. Accepted.

## Blockers / Challenges

None

## Commits

- 8c094c7 feat: warn when a non-scalar k on n A -> B is used as given (ticket 05)

## Time Spent

5m (ticket-start.js to ticket-complete.js)
