---
bump: patch
floor: patch
---
### Changed

- Faster runs when a species dips below 0 (overshoot relaxation, ADR 0007): the per-species tables of that path are built once with the matrices (`SimMatrices.overshoot`) instead of on every solver call, and the depletion switch and the mass-action / switch derivatives are vectorised. A pulsed 200-chunk run with `config_Labarbe` (`examples/bench_overshoot_pulsed.py`) goes from about 22.3 s to 14.9 s. Results are unchanged (checked against the former code at rtol 1e-12).
