---
bump: patch
floor: patch
---
### Fixed

- `radiopyo/data/reactions_biology.toml`: the superoxide dismutase reaction used `O2_r_minus`, a species
  that no reaction produces (the water radiolysis config names superoxide `O2_minus`), so SOD never
  acted. It now reads `2 O2_minus -> H2O2 + O2`. Results of configs that include this file change.
