---
bump: minor
floor: minor
---
### Removed (breaking)

- `force_positive` (`RunOptions`, `run`, `prepare_chunked_run`, `prepare_single_pulse`, `derive`, `jacobian`, `SimParameters`) and `zero_depleted_consumption`: passing `force_positive` now raises `TypeError`. Old configs and stores that contain the key still load.

### Changed

- Rate laws (mass action and Michaelis-Menten) read the clipped concentration max(c, 0), always; the analytic Jacobian column of a species with c < 0 is 0.
- The hard depletion block of "negative" products is replaced by a smooth switch c+/(c+ + eps), eps = `DEPLETION_LIMIT` = 1e-12 mol/L, applied once per negative product to k-reactions and enzymatic reactions; the Jacobian is smooth except at c = 0.
- Bundled catalase k in `reactions_biology.toml`: 1.324e8 -> 6.62e7, the per-reaction value used as given (assumption, not verified against Labarbe's model). Golden derive and runs data re-recorded.

### Added

- Building `n A -> B` with n > 1 and a non-scalar k (polynomial, Arrhenius, equilibrium) logs a warning on the `radiopyo` logger: k is used as given.
