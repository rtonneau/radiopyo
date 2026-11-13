# Implementation Plan

**Session:** ship rust radkin backend
**Date:** 2026-11-13T10:47:17.988Z
**Estimated effort:** 2 days

## Strategy

Build from the bottom up. First the `rust/` PyO3 crate (`radiopyo-rs`): it pins radkin by git rev and exposes one `run` covering all of radkin's modes. Next, the Python adapter for the single-run entry points with the parity guard and wrapping into `ODEResult`. Then the chunked and single-pulse wiring, then the parallel batch, then CI and release wheels, then the notebook and docs. Every Rust-dependent test is skipped when the extension is not installed, so each ticket leaves the pure-Python suite green. A radkin addition, if one turns out to be needed, is made in a branch + PR in `rtonneau/radkin` during ticket 01 and the run waits for the user's merge before pinning that rev.

## Tickets Overview

- **Ticket 1:** `rust/` maturin crate `radiopyo-rs` (abi3, py>=3.11) pinned to a radkin rev: `_core.run(config_toml, state, mode, options)` for run / run-till-max-dose / chunked / single-pulse, returning numpy arrays + metadata + labels and k/G vectors; radkin logs forwarded to Python logging; typed error kinds.
- **Ticket 2:** `UnitCell` keeps its merged `ConfigDict`; `radiopyo/rust.py` adapter (lazy import with install hint, parity guard, option checks, wrapping into `ODEResult`, error mapping); `backend=` on `UnitCell.run` and `run_till_max_dose`; equivalence tests vs Python.
- **Ticket 3:** `backend=` on `prepare_chunked_run(...).run()` and `prepare_single_pulse(...).run()` (point thinning, short ON spans, step dose rates), equivalence tests across configs × beams.
- **Ticket 4:** `radiopyo.rust.run_batch(cells, t_span, mode=, n_threads=, return_exceptions=, **options)`: rayon pool with the GIL released, input order, raise-after-all semantics; tests and a speed check.
- **Ticket 5:** CI job building the extension and running the Rust tests; release workflow building Windows/Linux abi3 wheels attached to the GitHub release; `radiopyo[rust]` extra and lockstep version files.
- **Ticket 6:** `examples/notebook_rust_backend.ipynb` (single run, batch LET sweep, timing), README and `docs/configuration.md` sections, CLAUDE.md architecture note, CHANGELOG fragment.

## Sequencing Rationale

Ticket 1 is the only Rust-facing interface and everything else calls it. Ticket 2 sets up the adapter, the guard and the `ODEResult` wrapping that tickets 3 and 4 reuse. Batch (4) depends on the per-cell preparation of 2–3. CI (5) needs the crate and the tests in place. Docs and the notebook (6) describe the finished API.

## Risks & Mitigation

- **Risk:** Python's merged `ConfigDict` does not round-trip into radkin's `Config` (shapes, numbers, missing keys) → **Mitigation:** the binding serializes the dict to TOML text (includes removed, already resolved) for `Config::from_toml_str`; the parity guard and the equivalence tests catch any drift. A radkin change goes through a radkin PR.
- **Risk:** git dependency chain (radkin → rsoda) and a long first compile slow CI and source installs → **Mitigation:** `Swatinem/rust-cache`, release builds only, abi3 wheels (one per OS).
- **Risk:** results differ between scipy LSODA and rsoda beyond rtol 1e-6 on some configs/modes (radkin's acceptance matrix passes 27 cases at 1e-6) → **Mitigation:** compare the same case matrix; a failing case is investigated against radkin's own golden data, never by loosening the tolerance.
- **Risk:** `ODEResult` features expect scipy fields (`sol`, `nlu`, `status`, `message`) or deepcopies, which breaks HDF5 persistence or `__add__` → **Mitigation:** fill them with neutral values in the wrapper; a test saves and reloads a Rust result through `H5Manager`.
- **Risk:** radkin and Python disagree on the beam built from `as_dict()` → **Mitigation:** tests cover each beam type built in code (e.g. `PulsedBeam.from_dose_rate_peak`) as well as from config.

## Assumptions

- Rust stable toolchain available locally (cargo 1.99) and on GitHub runners; maturin installed into the `radiobio_env` env (`C:\Users\rtonneau\AppData\Local\miniconda3\envs\radiobio_env\python.exe` on this machine).
- Network access to GitHub for the radkin and rsoda git dependencies.
- The user merges any radkin PR; the run waits for that merge before pinning.
- Python env edits after construction other than temperature / beam / LET are not forwarded (the parity guard catches those that change labels, k or G).
