# Ticket 01: enzymatic-depletion-block

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`: 558 passed (golden included), ruff "All checks passed!", mypy "no issues found in 42 source files".

## Review Notes

AC1: EzMatrices.negative_products (n_species, n_ez) bool, read-only, built from product `ope` in make_enzymatic_matrix; a negative product is now also subtracted in the enzymatic stoichiometry (as for k-reactions; no effect on bundled catalase). AC2: negative_product_block now accepts KMatrices or EzMatrices and is used for both in reaction_rates. AC3: tests/simulation/test_enzymatic_depletion.py (block on/off, vectorised, zero-reaction/bundled catalase); goldens unchanged. Also updated the EzMatrices constructor in tests/simulation/test_reaction_rates.py (outside the listed files, required by the new field). No review rounds.

Parent: base d4a17a796c25ed1f4f7fc80e5c3d2b401d483a6f, no pre-existing changes. Independent review skipped: localized change (one matrix field, one np.where) whose diff the parent read in full; verification re-run by the parent (558 passed, ruff and mypy clean). Stoichiometry sign fix accepted: a '- X' enzymatic product was previously added (+stoi), inconsistent with k-reactions; no bundled config affected. Accepted (implementer: sonnet, medium).

## Blockers / Challenges

None

## Commits

- c9052c6 fix: block depleted enzymatic reactions on negative products (ticket 01)

## Time Spent

5m (ticket-start.js to ticket-complete.js)
