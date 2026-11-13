---
bump: minor
floor: minor
---
### Added

- Optional Rust backend: the `radiopyo-rs` extension (`rust/`, PyO3 binding of radkin, the Rust port of
  radiopyo, pinned by git revision) solves a simulation from the cell's config and current state
  (temperature, beam, LET, initial concentrations) and returns a normal `ODEResult`. Select it with
  `backend="rust"` on `UnitCell.run`, `run_till_max_dose` and `SimChunksIter.run` (chunked and single-pulse
  runs, same chunk options, thinning and short ON spans). LSODA only; cells built from a config only.
- Parity guard (`radiopyo.rust.check_parity`): before solving, the species, reactions, k values, G-values,
  enzymatic constants, fixed and initial concentrations built by radkin are compared with Python's
  (rtol 1e-12); a mismatch raises `ParityError` (`RustBackendError` for the other radkin errors).
- `radiopyo.rust.run_batch(cells, t_span, mode=..., n_threads=None, return_exceptions=False, **options)`:
  parallel runs on a rayon pool with the GIL released. By default the first error is raised once every cell
  has finished; with `return_exceptions=True` each error takes its cell's slot.
- `rust` extra (`radiopyo-rs`; install it from a release wheel or `pip install ./rust`). Release wheels
  (abi3, Windows x64 and Linux x86_64) and an sdist are built on `v*` tags (`.github/workflows/wheels.yml`);
  the tests workflow gets a Rust job on Ubuntu and Windows.
- `examples/notebook_rust_backend.ipynb`: both backends side by side, a chunked pulsed run, a `run_batch`
  LET sweep and timings.
