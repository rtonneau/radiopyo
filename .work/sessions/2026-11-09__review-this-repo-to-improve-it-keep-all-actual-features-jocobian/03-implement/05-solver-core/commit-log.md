# Ticket 05: solver-core

**Status:** ✅ Done

## Local Test Result

Command (repo root, Anaconda Python on PATH):

```
$ python -m pytest -q && ruff check radiopyo && mypy radiopyo
287 passed in 23.09s
All checks passed!
Success: no issues found in 36 source files
exit=0
```

Evidence is in `.scratch/tests/2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian/evidence/`:

- `05-verify.txt`: the run above.
- `05-golden-seeds.txt`: `tests/golden` with PYTHONHASHSEED=0 and 7 gives 5 passed both times.
- `05-rates_old_vs_new.txt` (plus `rates_old.npz` and `rates_new.npz`): the base f1e3a99 tree (git archive) and the new tree are compared on Labarbe and pure water, 200 random states each. I checked `derive` for every force_positive × beam on/off combination, and `reconstruct_reactions` values, columns and column kinds. Every comparison is `bit_identical=True n_diff=0`.
- `bench_before.txt` / `bench_after.txt` (Labarbe):

  | | before | after |
  |---|---|---|
  | `derive` | 77.5 µs/call | 45.8 µs/call |
  | `reconstruct_reactions` | 0.975 s for 2409 pts | 0.011 s for 2385 pts |
  | `run` | 0.65 s | 0.39 s |

## Review Notes

- Parent: base f1e3a99933cb6b88968aa58b295362e47fb992e7; pre-existing changes: gps session state only.

Acceptance Criteria check:

- **Named matrices.** `SimMatrices` is now a frozen dataclass in `radiopyo/utils/sim_types.py`, so the import path is unchanged and no import cycle is introduced. It holds `k_reactions: KMatrices`, `radiolytic: GVector`, `enzymatic: EzMatrices` and `state: StateLayout`. All four are frozen and their arrays are made read-only in `__post_init__`. The old positional tuples and `mat["..."]` / `mat[i]` indexing are gone: grep finds no `Ki_mat`, `Gi_mat`, `Ez_mat`, `JRki` or `mat[0]` in `radiopyo/` or `tests/`.
- **One pure rate function.** `reaction_rates(t, y_full, mat, beam) -> ReactionRates` is used by `derive()` and by a vectorised `ODEResult.reconstruct_reactions()`. The latter makes one call over all time points instead of a per-time-point loop. The depletion threshold is the single constant `DEPLETION_LIMIT = 1e-12` in `ode_functions_solver.py`; `ode_result.py` no longer hardcodes it. The remaining `1e-12` in `unit_cell.py` are solver tolerances and docstrings, unrelated.
- **One unit/fixed-species helper.** `StateLayout.to_full_molar` / `to_dynamic_micromolar` handle 1-D and 2-D input and zero fixed species. Tests cover a round trip with 0 and 2 fixed species × 1-D/2-D, and a hand-built 2A→B system with no fixed species through `derive`.
- **Stale Jacobian removed.** `jacobian()`, `EPSILON` and `JRki` are removed. `UnitCell._raw_run` raises `NotImplementedError("The analytical Jacobian is being rewritten ...")` when `use_jac=True`. This covers `run` and the chunked runs, and it is tested for `run` and `prepare_chunked_run(...).run()`. The `use_jac` docstrings were updated.
- **Golden tests.** `tests/golden/test_golden.py` changes only in `_derive_adapter`, which now calls `derive(t, y, beam, mat, force_positive)`. All golden tests pass with no tolerance change.
- **Focused tests.** The new `tests/simulation/test_reaction_rates.py` has 25 test cases, some parametrized over Labarbe and pure water. They cover:
  - frozen/read-only matrices and shapes
  - stoichiometry vs reaction objects (k and G)
  - mass action and Michaelis-Menten laws
  - `dose_rates`
  - `derive` == stoichiometry @ rates (beam on/off)
  - `reconstruct_reactions` vs per-time-point `reaction_rates`, plus column order k→G→Ez, zero G when the beam is off, the cache, and the no-time-point case
  - the negative-product block at the `DEPLETION_LIMIT` boundary (1-D and 2-D)
  - `zero_depleted_consumption`
  - `use_jac=True`
- **Behaviour preserved.** Results are bit-identical to base f1e3a99 (see Local Test Result).
- **Public API.** The `UnitCell` methods and kwargs are unchanged. `ODEResult` outputs are unchanged, including the `reconstruct_reactions` column order k→G→Ez. That order differs from the global index order k→Ez→G; this is documented in its docstring and in `ReactionRates.stacked()`. The HDF5 format is untouched.

Notes for review:

- **Ez stoichiometry bug, recorded and not fixed.** The old `Ez_iy` is built with `=` instead of `+=`, and this is kept in `make_enzymatic_matrix` with a NOTE comment. It is wrong only when one species sits on both sides of an enzymatic reaction; the bundled catalase reaction is unaffected. Enzymatic reactions also ignore the negative-product `ope`, as before. Both are candidates for a Phase C ticket.
- **NaN/inf corner case.** The old one-hot einsum turned NaN/inf anywhere in the state into a factor of 1 for every reaction. The new index gather propagates NaN/inf only to the reactions that involve the species. This only matters for an already divergent state.
- **Constant species in the G stoichiometry.** Their rows are now 0, as for k-reactions. Constant species are stripped from the output, so nothing changes.
- **G term formula.** The G term is computed as `S_g @ r_g` rather than `(sum of yields) * dose_rate`. The two are identical when a species has at most one radiolytic reaction, which holds for every bundled config.
- **Deepcopy.** `ODEResult` deep-copies `mat`, and the copied arrays lose the read-only flag. The dataclasses stay frozen, so this is harmless.
- **Memory.** The 2-D gather in `reconstruct_reactions` uses about 4× the output size as temporary memory, against roughly 90× faster runtime.
- **Internal API break.** The signature of `derive` changed to `(t, y, beam, mat, force_positive)`. `jacobian` is removed, `SimMatrices` is no longer a TypedDict, and the `make_*` builders return the new dataclasses.
- **Deliberately left alone (other tickets):**
  - `make_cst_species_vector` is still called without `T` (ticket 13).
  - The missing `raise` on `NoBeamDefinedError` in `_raw_run` and the `run()` argument-order bug are untouched (ticket 12).
- **Docs now stale (outside this ticket's file list, for the parent).** `CLAUDE.md` architecture items 3 and 4 still name `Ki_mat`/`Gi_mat`/`Ez_mat`/`cst_sp` and the stale `jacobian()`. `CHANGELOG` has no entry for this ticket.
- **Formatting.** I ran `ruff format` only on the new test file. Package formatting was already inconsistent and I did not reformat it, to keep the diff readable.

Nothing changed after review.

- Parent acceptance: independent reviewer (opus/xhigh): no findings; parent re-ran pytest (287 passed), ruff 0, mypy 0; old-vs-new bit identity evidence in evidence/05-rates_old_vs_new.txt. Carried forward: CLAUDE.md/CHANGELOG stale (ticket 18), Ez_iy = vs += (Phase C), NaN now propagates per-reaction only (diverged state only).

## Blockers / Challenges

None blocking. Minor issues along the way:

- **Bit identity.** I proved it against a `git archive` copy of f1e3a99 under a fixed PYTHONHASHSEED, because set/dict-ordered species can change reaction order between processes.
- **Gather vs einsum.** I chose the index gather with a padding row of 1.0 over the old NaN one-hot einsum because it gives the exact same product order, and so bit-identical floats.
- **Scratch-script slips (fixed).** One script passed the wrong keyword to `PulsedBeam`: `dose_rate_peak` is correct, not `dose_rate`. Another had a docstring escape `/!\\` mangled by a shell heredoc.

## Commits

- f567a6f refactor: named solver matrices and shared reaction_rates (ticket 05)

## Time Spent

35m (ticket-start.js to ticket-complete.js)
