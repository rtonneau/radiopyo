---
bump: minor
floor: patch
reason: breaking removals (RON format, run_multi_doses, e_LET, use_jac) in a 0.x version
---
### Changed

- Internal refactoring for readability (solver matrices, run flow, config readers, beams, results); public names and call signatures are kept and results are unchanged unless listed under Fixed. See the Unreleased section for the details.

### Fixed

- `UnitCell.run` now applies the initial radiolytic burst and passes `force_positive` to the solver; a `y0` given by the caller is never modified.
- One temperature source of truth: the `temperature` of a config and `set_temperature` now reach fixed species, initial values and results.
- `[includes]`: every section is merged, a later include wins, and relative paths resolve next to the including file.
- Beams: `MacroPulsedBeam` no longer hangs without `max_dose`, respects `stop` and counts partial micro pulses; average dose rates account for `start_time`.
- Species order no longer depends on `PYTHONHASHSEED`, so runs are reproducible.
- Arrhenius reactions load from config (`A`, `Ea` in kJ/mol); enzyme names with capitals parse; enzymatic `k_value` survives the HDF5 round trip.

### Removed

- RON config format, `UnitCell.run_multi_doses`, `radiopyo.utils.e_LET` and the stale Jacobian (`use_jac=True` raises `NotImplementedError`).
