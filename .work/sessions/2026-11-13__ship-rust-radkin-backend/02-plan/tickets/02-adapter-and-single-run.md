# Ticket 02: adapter-and-single-run

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] `UnitCell` stores the merged `ConfigDict` built by `from_toml` / `from_parser` / `from_config_dict` (`None` for `UnitCell(env=...)`); deepcopies and existing behaviour unchanged.
- [ ] `radiopyo/rust.py` (public `radiopyo.rust`): lazy import of `radiopyo_rs` raising `ImportError` with both install commands; `RustBackendError` (deriving from `RadioPyoError`) and `ParityError`; option checks (LSODA only, unsupported options raise); parity guard (labels exact, k / G at rtol 1e-12 vs `uc.mat`); wrap into an `ODEResult` with the Python env/mat/beam/LET; radkin error kinds mapped to `SolverFailedError`, `NoBeamDefinedError`, `ShortBeamOnError`, ...
- [ ] `UnitCell.run(..., backend="python"|"rust")` and `run_till_max_dose(..., backend=)`; default stays `"python"`.
- [ ] Tests (`tests/rust/`, skipped without `radiopyo_rs`): run mode vs Python at rtol 1e-6 on pure water / Labarbe / paper config × constant / single-pulse / pulsed / macro-pulsed beams, `t_eval`, temperature and LET changes, a beam built with `from_dose_rate_peak`; parity guard raises after editing a k; `UnitCell(env=...)` raises; `to_pandas`, `reconstruct_reactions` and an HDF5 save/reload work on a Rust result. A test without the extension checks the `ImportError` hint.
- [ ] `ruff check radiopyo`, `mypy radiopyo`, full pytest green.

**Files to Touch:**
- `radiopyo/rust.py`
- `radiopyo/simulation/unit_cell.py`
- `radiopyo/simulation/exceptions.py`
- `radiopyo/__init__.py`
- `tests/rust/test_rust_run.py`
- `tests/rust/conftest.py`

**Verification Step:**

Run:
```bash
<radiobio_env python> -m pytest tests/rust -q && <radiobio_env python> -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
all Rust tests pass (none skipped with the extension built), the full suite passes, ruff and mypy report zero errors.

**Notes:**

Keep the module banner layout and full typing; `radiopyo_rs` has no stubs, so add a small `.pyi` in `rust/` or a typed wrapper. Use radiopyo's `as_dict()` and the beam type name from `BEAM_TYPES`. Fill the scipy-only `ODEResult` fields (`nlu`, `status`, `message`, `success`) with neutral values.
