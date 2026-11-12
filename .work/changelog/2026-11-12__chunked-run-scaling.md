---
bump: minor
floor: minor
---
### Added

- Point thinning of the chunked runs (opt-in): `prepare_chunked_run` and `prepare_single_pulse` take
  `points_per_chunk_on` / `points_per_chunk_off` (None by default: every solver step is stored). With an
  integer N >= 2, N log-spaced points are stored per chunk (start and end kept, first point 1e-6 of the span
  after the start); the solver steps, hence the state passed from chunk to chunk, are unchanged. The
  production / consumption / species integrals of such a result are computed over the stored points only
  (coarser): `ODEResult.thinned` is set and a warning is logged once per result.
- `examples/bench_chunked_scaling.py`: wall time and cost per chunk of pulsed chunked runs versus the number
  of pulses (config_Labarbe), with optional thinning.

### Changed

- Chunked runs scale linearly with the number of chunks: `ODEResult` stores every chunk as a block and
  joins `t` / `y` only when they are read (it used to copy the whole history at every chunk), and the next
  chunk starts from `final_cc`, read from the last block. `derive` allocates less per call (1-D fast paths,
  no depletion switch without negative products, no `flatnonzero` when no species is below 0): about 35%
  faster per call, results bit-identical. 1600 pulses on config_Labarbe: 455 s -> 160 s, about 50 ms per
  chunk whatever the number of chunks.
