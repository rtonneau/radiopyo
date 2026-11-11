---
bump: minor
floor: minor
---
### Added

- `SolverFailedError` (`radiopyo.simulation.exceptions`): `run`, chunked and single-pulse runs raise it when `solve_ivp` returns `success=False` (with the span, the last reached t and scipy's message), instead of continuing silently from a truncated solution.
- `ODEResult.step_dose_rates`: the dose rate each step of a chunked run was integrated with, used by `reconstruct_reactions` (not stored in HDF5; results read from a file keep reading the beam). New `ConstantDoseRate` / `DoseRateSource` in `radiopyo.physics.beam` and `SimChunk.dose_rate`.
- Benchmark `examples/bench_fenton_stall.py` (Fenton reproduction grid) and the regression test `tests/simulation/test_fenton_stall.py`.

### Changed

- Overshoot relaxation (ADR 0007, amends ADR 0004): a dynamic species driven below 0 is pulled back to 0 at the rate it would be consumed at c = 0+ (`dc/dt += J_jj(0+) * c`, included in `jacobian`); every other species still reads it as 0. Results are identical when all concentrations are >= 0. `derive`/`jacobian` are about 15x slower on states with species below 0.
- Chunks follow the beam exactly: `make_time_intervals` no longer extends a pulse by `GROWTH_TIME` (removed), and chunked runs (`prepare_chunked_run`, `prepare_single_pulse`) hold the dose rate of each chunk constant. `prepare_single_pulse` now builds one chunk per pulse, adds a leading OFF chunk when the beam starts after `t_span[0]` and clips to `t_span[1]`. Golden beam/run data re-recorded (time grid only; species agree to 7.6e-9).

### Fixed

- LSODA stalled for minutes (and runs froze, then resumed) with a stiff radical near 0, e.g. Fenton with RH = 1 M and O2 = 50 µM; the dose-per-pulse × dose-rate maps of ROOH showed "islands". The Fenton benchmark now completes every point in under 30 s (a point that took 81 s takes 0.9 s).
- `prepare_single_pulse` with several pulses before the max dose, or with a macro-pulsed beam, now delivers the right dose.
