# Ticket 17: fix-arrhenius-enzymes-add

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`

- pytest: 547 passed (530 baseline + 17 new tests), golden tests unchanged and passing.
- ruff: All checks passed. mypy: Success, no issues in 42 source files.
- `git status`: `golden_results.h5`, `.golden_results/*.json` and `hdf_expected.npz` are unmodified.
- Evidence: `.scratch/tests/2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian/evidence/17-before-after.md`

## Review Notes

- Parent: base cb56f469edf693087161c3eaf4a682d92aba25f4; pre-existing changes: gps session state only.

AC1 Arrhenius: `Arrhenius.k_value` = A*exp(-Ea*1000/(R*T)), T in K from the simulation temperature (same `Temperature` path as before); tests: hand-computed 3.1344e6 at 25 degC and 4.2827e6 at 37 degC (class, from_generic, TOML file through `UnitCell.from_toml`). Extra bug found and fixed: the "arrhenius" section builder required `k_value`, so a config entry never loaded (now `("A", "Ea")`). The old meaning of `Ea` (Ea/R in K, used as `exp(-Ea/T)`) cannot coexist with kJ/mol under the same key; the README and CHANGELOG give the conversion (x 8.314e-3). No bundled config uses arrhenius.
AC2 Enzyme names: grammar `enzyme: element` (a name only). `SOD`, `GPx`, `Cat2`, `E` and lowercase `catalase` now give a `BasicSpecies` with the right label; tests in `tests/parser/test_reaction_parser.py`.
AC3 `ODEResult.__add__` returns a new object (shallow copy, every list/array rebuilt); test with deepcopy snapshots of both operands.
AC4 Goldens: all pass unchanged. Extra item `MichaelisMenten.to_dict` now writes `k_value` (as `from_generic` reads it): round trip 1.324e8 -> 1.324e8 (was 2.648e8, rebuilt rate ratio 2.0 -> 1.0), test added.
DEVIATION: `hdf_expected.npz` was NOT modified. The HDF golden test reads the committed old-code json (k_value 2.648e8), so its expected rebuilt rate is still the doubled one and halving it would make the test fail. Only the stale `compute_hdf` docstring in `tests/golden/test_golden.py` was updated. Regenerating the npz via compute_hdf would also change unrelated arrays (n_points 40 vs 38).
FOLLOW-UP (user decisions): (1) `species_balance_evolution` FIXED: it now applies the k, radiolytic and enzymatic stoichiometry matrices to the reaction rates and honours `start`/`stop`; H2O2 at 2.59e-3 s went from +3.53e-9 to -1.1518e-5 mol/L/s (finite difference -1.1520e-5); tests compare it with the slope of the solution, check the catalase term and `start`/`stop`. (2) Enzymatic stoichiometry `=` changed to `+=`/`-=`: `Aa -- Ee >> Aa + Bb` gives net 0 on Aa (test added), NOTE comment updated; goldens unchanged. (3) NOT fixed, deferred to the Jacobian work: enzymatic reactions ignore the negative-product flag (no block in `michaelis_menten_rates`, `ope` not read in `make_enzymatic_matrix`; `Aa -- Ee >> Bb - Aa` gives Aa = +1, Bb = +1). Numbers in the evidence file.
Tooling note: edits keep the files' CRLF working-copy line endings.

- Parent acceptance: implementer sonnet/medium; USER chose to fix findings 1 (species_balance_evolution missing radiolytic/enzymatic terms, start/stop ignored) and 2 (Ez stoichiometry = -> +=/-=), defer 3 (enzymatic negative-product flag) to the Jacobian work. Catalase to_dict fix included as approved; hdf_expected unchanged since the golden reads old-code json. Reviewer sonnet/medium: no findings. Parent re-ran full suite (551 passed), ruff 0, mypy 0, golden data untouched.

## Blockers / Challenges

None blocking. A first Python edit script wrote a non-ASCII character with the default cp1252 encoding in `michaelis.py` (SyntaxError at collection); replaced by ASCII, all files re-checked. Review rounds: 0.

## Commits

- 2ad881e fix: Arrhenius A/Ea in kJ/mol, enzyme names, ODEResult.__add__, enzymatic to_dict and balance (ticket 17)

## Time Spent

1h 00m (ticket-start.js to ticket-complete.js)
