# Ticket 02: adapter-and-single-run

**Status:** ✅ Done

## Local Test Result

`pytest tests/rust -q`: 24 passed, none skipped (radiopyo-rs built, rev 5e0d181).
`pytest -q`: 834 passed (output in `.scratch/tests/2026-11-13__ship-rust-radkin-backend/02-pytest.txt`).
`ruff check radiopyo`: all checks passed. `mypy radiopyo`: no issues in 43 files.
Ad-hoc check (`02-paper.py`, pure water with GF(DR) on, O2 changed on the env, single pulse at 1e6 Gy/s): max relative difference of the final concentrations between the Python and Rust backends is 2.9e-10.

## Review Notes

Base: ac7fb3b (session record of ticket 01, on top of d45a37b). Pre-existing changes: none besides the session files; `.gitignore` also gains the rule ignoring the built `_core.pyd`/`.pdb` that ticket 01's `maturin develop` leaves in `rust/python/radiopyo_rs/`.
AC 1: `UnitCell.config` holds a deepcopy of the merged ConfigDict when the cell is built by `from_config_dict` (and so by `from_parser`/`from_toml`), and is None for `UnitCell(env=...)`. `copy.deepcopy` copies it (tested).
AC 2: `radiopyo/rust.py` provides:
- a lazy `_core()` import, whose `ImportError` gives both install commands (tested with the extension hidden);
- `RustBackendError` and `ParityError` in `simulation/exceptions.py`;
- LSODA only. A cell that overrides `FIRST_STEP`, `T_ZERO` or `MIN_BEAM_ON_TIME` (constants radkin hard-codes) raises, and radkin rejects unknown options;
- a parity guard comparing species labels, reaction labels (radkin's groups reordered to Python's global order: k-reactions, enzymatic, radiolytic), k, G, enzymatic k and Km, fixed concentrations and initial values at rtol 1e-12;
- the result wrapped as an `ODEResult`, one block per radkin step, with neutral scipy fields;
- the error mapping solver_failed→SolverFailedError, no_beam→NoBeamDefinedError, short_beam_on→ShortBeamOnError, everything else→RustBackendError.
Added beyond the ticket: `radkin_config(uc)` folds `env.parameters` and the scalar initial concentrations into the config sent. The paper workflow edits the env after the build (`update_flag_GValues(True)`, `set_initial_cc`); without this the guard would reject it, and the guard still checks the outcome.
AC 3: `run(..., backend=)` is appended last with default "python" and checked by `check_backend`. `run_till_max_dose` forwards it through `KwargsRun` (and so runs radkin's run mode over the same t_span as Python). `test_public_signatures_unchanged` was extended for the appended keyword.
AC 4: the beams cover constant, single pulse, pulsed, macro-pulsed and `from_dose_rate_peak`, each on pure water and Labarbe. They also cover the paper setup (pure water with GF on and O2 edited, the paper reproduction's config), `t_eval`, a temperature change and an LET change, `y0` (not modified), `use_beam=False`, `run_till_max_dose` and the short ON spans. The parity guard is tested by doubling a k. The no-config cell, a BDF method and the no-beam case each raise; `to_pandas`, `reconstruct_reactions` and an HDF5 push/read also run on a Rust result.
The pulsed cases pass `max_step` below the ON time: a plain run reads the beam at every evaluation and LSODA otherwise steps over pulses in a solver-dependent way, which differed between the backends by up to 50 %. That is a property of plain runs on pulsed beams, not a backend bug.
AC 5: green (see above).

## Blockers / Challenges

- Python's global reaction order is k-reactions, enzymatic, radiolytic, while radkin lists k-reactions, radiolytic, enzymatic. The guard compares the groups after reordering.
- `beam_type()` returns display names ("Constant"), so the BEAM_TYPES key of the class is used instead (radkin reads `macro_pulsed`).

## Commits

- 98f5e6a feat: backend="rust" on UnitCell.run via radiopyo.rust, with parity guard (ticket 02)

## Time Spent

12m (ticket-start.js to ticket-complete.js)
