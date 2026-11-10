# Ticket 18: docs

**Status:** ✅ Done

## Local Test Result

- `python -m pytest -q && ruff check radiopyo && mypy radiopyo` -> `551 passed in 98.49s`, `All checks passed!`, `Success: no issues found in 42 source files` (exit 0; evidence/18-verify.txt).
- README check: `.scratch/tests/<session>/t18/check_readme.py` runs each of the 3 python blocks of README.md in a fresh subprocess (cwd = temp dir): 3 PASS (evidence/18-readme-check.txt, exit 0). The 16 toml blocks are concatenated into one user file, loaded with `UnitCell.from_toml` (all beam types incl. `macro_pulsed`, acid_base, equilibrium, arrhenius, enzymatic, parameters): PASS. Skipped by language tag: 2 bash blocks (pip commands). Harness-only change: the chunked block (max_dose=10 Gy = 2000 pulses, about 1e6 solver steps) is run with `t_span=[1e-7, 0.05]` instead of `[1e-7, 10]`; the README text keeps 10 s.
- Seeds: evidence/18-followups.md (F1-F14) and evidence/18-seeds.json (8 candidates, not merged).

## Review Notes

- Parent: base a09d989de36e58ebaf415688850896ff4cbfae8b; pre-existing changes: gps session state only.

- AC1 (README): version badge 0.6.0; `from_toml` (was `from_toml_file`); `basic_water_radiolysis.toml` references replaced (reactions_biology.toml has no includes; Labarbe includes advanced_water_radiolysis + reactions_biology); `max_dose=` keyword; `[[beam.macro_pulsed]]` documented (and `frequency` / `pulse_dose` variants); acid_base keys checked against `_build_acid_base` / `AcidBase` (reaction, pKa, k_forward and k_backward both required, k_value rejected); `[[reactions.equilibrium]]` and the `[parameters]` keys added; G-value list with H3O_plus = e_aq (ADR 0002) and `variable_GValues` (factor and unchecked range); dependencies (h5py added, pandas-stubs to dev extra, tqdm optional extra, Python >= 3.11); GitHub URL; n*k convention (scalar only, ADR 0001); stale "no tests" remark replaced.
- AC2: see Local Test Result.
- AC3: CHANGELOG got Added / Changed / Known issues sections and the Removed list was completed (refactor tickets 02-10: matrices, RunOptions/_raw_run, table-driven builders, utils split, chemistry cleanup, beam split, exceptions, hygiene/logging, removed jacobian and unused members); checked against `git log main..HEAD`. A stale sentence (json round trip of enzymatic reaction still doubles k_value) was removed, ticket 17 fixed it. CLAUDE.md: Current state (551 tests, ruff/mypy zero) and architecture items 1-7 + Logging rewritten and verified against the code (SimMatrices/KMatrices/GVector/EzMatrices/StateLayout, `derive(t, y, beam, mat, force_positive)`, `reaction_rates`, no jacobian, module splits, include semantics, single temperature source, deterministic species order, `make_hash` broken).
- AC4: seeds recorded, each finding re-checked against the code first (NaN timer, make_hash TypeError, strip, H2O polynomial values, `_scale_factor(0)` = inf all reproduced; none already fixed). Ticket 17 fixed the enzymatic `k_value` doubling, so it is only listed as an old-file caveat (F12).
- Notebook (read-only): no removed or renamed API is called. Cell 19 passes `max_dose` positionally to `PulsedBeam.from_dose_rate_peak` (lands in `n_pulse`); cells 25/27 contain a corrupted micro sign; outputs are stale. Not edited.
- README.md and CHANGELOG.md keep their CRLF line endings. Nothing staged or committed. Review rounds: 0 (no independent reviewer run by the implementer).

- Parent acceptance: implementer sonnet/medium (plan said fast; chosen for accuracy of CLAUDE.md and README checks). Reviewer skipped (docs only); parent read the CLAUDE.md rewrite, re-ran the README snippet check (exit 0) and the full suite (551 passed). Seeds merged by parent with scout-merge.js --from evidence/18-followups.md: 9 seeds (incl. refresh-example-notebook for cell 19 positional max_dose).

## Blockers / Challenges

None. The first run of the README check exceeded the tool timeout (the README chunked example now really delivers 10 Gy, i.e. 2000 pulses, because `max_dose=` replaced the positional `n_pulse`); solved by shortening its `t_span` in the harness only. A bash heredoc holding quoted Python failed, the edit scripts were written with the file tool under the scratch dir (`t18/edit_*.py`).

## Commits

- 4ea057d docs: README, CHANGELOG and CLAUDE.md for the refactor session (ticket 18)

## Time Spent

20m (ticket-start.js to ticket-complete.js)
