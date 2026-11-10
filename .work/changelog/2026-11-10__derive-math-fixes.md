---
bump: major
floor: major
---
### Removed

- `force_positive` (`RunOptions`, `run`, `prepare_chunked_run`, `prepare_single_pulse`, `derive`, `jacobian`, the `[parameters]` key) and `zero_depleted_consumption`: passing `force_positive` now raises `TypeError`. Old configs and stores that contain the key still load. Negative concentrations are handled by clipping instead (see Changed).

### Changed

- Rate laws (mass action and Michaelis-Menten) read the clipped concentration max(c, 0), always: a species overshooting below 0 is neither consumed further nor runs a reaction backwards. The analytic Jacobian column of a species with c < 0 is 0.
- The hard depletion block of "negative" products (rate set to 0 below `DEPLETION_LIMIT`) is replaced by a smooth switch c+/(c+ + eps), eps = `DEPLETION_LIMIT` = 1e-12 mol/L, applied once per negative product to k-reactions and enzymatic reactions; the Jacobian is smooth except at c = 0. The bundled configs have no negative product, so their results do not change.
- Bundled catalase k in `reactions_biology.toml`: 1.324e8 -> 6.62e7, the per-reaction value, used as given (Michaelis-Menten k is never divided by the stoichiometry). This reading of Labarbe's value is an assumption, not verified against his model. Results of configs including `reactions_biology.toml` change: in the Labarbe run the H2O2 plateau roughly doubles (e.g. 0.177 -> 0.334 umol/L).

### Added

- Building `n A -> B` with n > 1 and a non-scalar k (polynomial, Arrhenius, equilibrium) logs a warning on the `radiopyo` logger: such a k is used as given, not divided by n (n*k convention).
