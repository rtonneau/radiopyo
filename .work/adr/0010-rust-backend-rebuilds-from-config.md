# The Rust backend rebuilds the system from the config, guarded by a parity check

_ADR 0010, recorded 2026-11-13 (session 2026-11-13__ship-rust-radkin-backend)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

With `backend="rust"`, radkin receives the merged config dict plus the current temperature, beam and LET and runs the whole simulation itself (bursts, chunking, short ON spans), rather than receiving Python's `SimMatrices` and only replacing `solve_ivp`: this reuses radkin's golden-tested drivers instead of adding a matrices-from-arrays path. The cost is two builds of the same system, so before every Rust run a parity guard compares radkin's species/reaction labels (exact) and k/G vectors (rtol 1e-12) with Python's matrices and raises on any mismatch; a `UnitCell` built directly from a `SimEnv` (no config) cannot use the Rust backend.
