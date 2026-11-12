---
bump: patch
floor: patch
---
### Added

- Logo (`docs/logo.svg`), a citation section (Tonneau et al. 2025, Phys. Med. Biol. 70 235021, DOI 10.1088/1361-6560/ae1fcc, with BibTeX) and a GitHub Actions `Tests` workflow (ruff, mypy, pytest on Python 3.11).

### Changed

- README shortened to install, update, usage, a condensed configuration table and citation; the full TOML, beam and run reference moved to `docs/configuration.md`.
- Ruff rule set is now an explicit `select` (stable across ruff versions); mypy ignores missing h5py stubs and no longer reports unused ignores.

### Fixed

- `ODEResult.integrate_production_of` / `detail_production_of` raised `IndexError` under pandas 3 (columns dropped while iterating).
