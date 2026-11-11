---
bump: major
floor: minor
reason: beams with a configured ON time below 1e-6 s, accepted before, now raise ShortBeamOnError
---
### Changed

- A beam whose configured ON time (pulse `on_time`, `micro_on_time`, macro train `on_time`, single pulse, constant beam limited by `max_dose`) is shorter than `MIN_BEAM_ON_TIME` = 1e-6 s (`radiopyo.physics.beam`, also `UnitCell.MIN_BEAM_ON_TIME`) raises `ShortBeamOnError` (`radiopyo.physics.exceptions`, a `BeamError` and a `ValueError`) when it is built or its parameters are updated. Such beams were accepted before.
- `prepare_single_pulse`: the OFF chunks before the time of max dose use `max_step_size_trans` (they used `max_step_size_final`).

### Fixed

- A beam ON span shorter than `MIN_BEAM_ON_TIME` met at run time (last pulse cut by `max_dose`, pulse clipped by `t_span`) no longer raises `ValueError("Too Short beam pulse")` in chunked, single-pulse and plain runs: it is not integrated, its radiolytic yield is added analytically, a warning is logged, and it is recorded in the new `ODEResult.short_on_spans` (`ShortOnSpan`: start, stop, dose rate; not stored in HDF5).
