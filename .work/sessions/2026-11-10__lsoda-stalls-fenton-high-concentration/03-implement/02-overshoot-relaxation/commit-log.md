# Ticket 02: overshoot-relaxation

**Status:** ✅ Done

## Local Test Result

Command: `python -m pytest -q && ruff check radiopyo && mypy radiopyo`

Result: `656 passed in 96.11s` (golden tests included, no drift, nothing re-recorded), `All checks passed!` (ruff), `Success: no issues found in 42 source files` (mypy). Output kept in `.scratch/tests/2026-11-10__lsoda-stalls-fenton-high-concentration/evidence/ticket02-verification.txt`.

Extra check: on the Labarbe and pure-water configs, `derive` for a random state with ~30% negative species equals the prototype formula (`derive` of the clipped state + diag(jacobian at the clipped state) * y on the negative rows) to 4e-16 relative.

## Review Notes

Review rounds: 0 (implementer self-check only).

- AC1 (`derive`): `derive` adds `overshoot_slopes(...) * c_j` on the rows of the dynamic species with y_j < 0; J_jj(0+) is computed with c_j set to 0 (0+ side), other negative species keep reading as 0. Other rows untouched (`test_other_rows_are_unchanged`, exact equality with `derive` of the clipped state).
- AC2 (`jacobian`): diagonal = J_jj(0+), row j gets c_j * dJ_jj(0+)/dc_k (exact second derivatives, including negative-product switch, enzymatic E/S terms, repeated reactant). Matches central finite differences on toy k/enzymatic envs with one or several negative species (`test_overshoot_relaxation_k_reactions`, `_enzymatic`, `test_negative_cell_with_clipped_species`) and entry by entry (`test_overshoot_relaxation_cross_derivative_value`).
- AC3 (no full Jacobian): `_overshoot` loops over the negative species only and evaluates only the reactions with S_ji != 0 (sub-matrices), plus the derivatives only when called from `jacobian`. A state with no negative species skips it (`np.flatnonzero(y < 0)` empty): identical to S @ rates (`test_state_without_negative_species_is_s_times_rates`). Cost with 6 negative species on the pure-water config: about 0.7 ms per `derive` against 0.06 ms without negative species.
- AC4: `reaction_rates` / `ODEResult.reconstruct_reactions` untouched; module docstring and `derive` / `jacobian` docstrings explain the term and cite ADR 0007.
- AC5 tests (new `tests/simulation/test_overshoot_relaxation.py`): stiff A -> B decay below 0 restored (derive value and LSODA integration), species with no consumption stays put, 2 A -> E has zero slope, other rows unchanged, negative-product and enzymatic-substrate slopes.

Changed after the change of behaviour (ADR 0007 amends ADR 0004 for the rows of negative species), existing tests updated: `test_clipped_concentrations.py` (`test_reverse_reaction_not_run_backwards` now expects the C row = -k*C and 0 elsewhere; the stoichiometry test adds the relaxation term using a forward-difference slope, rtol 1e-5 on states with negative species; the conservation test skips combinations containing a negative species) and `test_jacobian.py` (`test_negative_species_contributes_zero`, `test_negative_enzymatic_species_contribute_zero`, `test_negative_cell_with_clipped_species` allow the own-row diagonal; new FD tests for the relaxation). Note: the ticket lists `test_jacobian.py` but not `test_clipped_concentrations.py` under Files to Touch; the latter had to change since its assertions encoded the old behaviour.

Remaining gap: the chunk / beam parts of ADR 0007 are not in this ticket.

## Blockers / Challenges

- The ticket named `proto.py`; the prototype is in `proto3.py` (as the parent said). It omits the y_j * dJ_jj/dy_k term; the exact term was derived (effective law T_i = law_i/eps if j is a negative product of reaction i else d law_i/dc_j, times the switch without the factor of j) and verified against finite differences.
- Two existing test files encoded the pre-ADR-0007 behaviour (see Review Notes) and were updated.
- `derive` is about 10x slower on states with several negative species (Python loop over the negative species); not slower otherwise.

## Commits

- d82357a fix: overshoot relaxation of species below 0 in derive and jacobian (ticket 02)

## Time Spent

11m (ticket-start.js to ticket-complete.js)
