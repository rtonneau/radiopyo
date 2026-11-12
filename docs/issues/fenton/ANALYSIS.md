# Fenton "runaway" with Fe2+ = 2.2e-8 M: analysis (2026-11-12)

Case: `Configuration_Bordet_Emballement.toml` + this folder's `reactions_biology.toml`, 1200 s,
`max_step_size_on=1e-4`, `max_step_size_off=0.5`, run on radiopyo main after the PR
"Fix LSODA stalls near zero concentration". Reproduction scripts:
`.scratch/tests/2026-11-12__fenton-runaway-after-overshoot-fix/`.

## 1. The freeze and the islands are fixed

No LSODA convergence failure in any run. Scan around `UHDR_Freeze` (0.4 Gy/s peak):

| Frequency (Hz) | ROOH at 1200 s, Fe = 0 (µM) | ROOH at 1200 s, Fe = 2.2e-8 (µM) |
|---|---|---|
| 0.10 | 1.56 | 2.82 |
| 0.12 | 1.53 | 2.79 |
| 0.13 (`UHDR_Freeze`) | 1.76 | 3.01 |
| 0.14 | 1.77 | 3.02 |
| 0.16 | 1.74 | 2.99 |

The small step between 0.12 and 0.13 Hz is there without Fe as well, so it does not come from Fenton. The islands and the freeze were solver artifacts.

## 2. The ROOH growth is chemistry, and it is linear

- Catalase = 0 in this file, so H2O2 stays at about 1.06 µM for the whole run.
- Fe2+ is a fixed species and is never consumed: the Fenton OH• flux is 1e3 × 2.2e-8 × 1.06e-6
  ≈ 2.3e-5 µM/s, and it stays constant (H2O2 lifetime ≈ 4.5e4 s).
- `ROO_r + RH -> ROOH + R_r` (1.25 /s per ROO•) regenerates R_r, while the exits are slow
  (`ROO_r -> x`, 0.01 /s, and self-reaction). Chain length ≈ 50, ROO• ≈ 1e-3 µM,
  so ROOH grows at ≈ 1.2e-3 µM/s (3e-6 µM/s without Fe).
- Fe3+ accumulates (0.028 µM at 1200 s): nothing recycles it.

## 3. Levers (UHDR, 50 Gy/s peak, 25 Hz)

| Variant | ROOH at 1200 s (µM) | End slope (µM/s) |
|---|---|---|
| As is | 2.40 | 1.2e-3 |
| Fe2+ consumed (moved to `[concentrations.initial]`) | 1.96 | 6.8e-4 |
| + Fe3+ + O2_minus -> Fe2+ + O2 (k = 1.5e8, assumed) | 2.01 | 7.1e-4 |
| GSH = 6.5e-3 M | 1.49 | 5.2e-4 |
| Catalase = 8e-8 M | 1.11 | 2.7e-6 |
| All of the above | 0.98 | 2.7e-6 |

CONV (0.16 Gy/s, 30 Hz), as is: 6.61 µM, slope 1.1e-3 µM/s. The CONV lever variants were not finished when this was written.

Catalase is the lever that stabilises ROOH, because it removes the H2O2 that feeds Fenton.

## 4. Can the runaway results be used?

Not as they are. With no H2O2 sink and Fe2+ that never runs out, the final ROOH depends on when the run stops,
not on the beam. Comparing UHDR and CONV at the same period would only compare two linear drifts.
Use a config whose ROOH settles (catalase on, or more antioxidants), and scale a tolerable Fe2+
concentration by the healthy/tumour ratio.

## 5. Bug found in the bundled config

`radiopyo/data/reactions_biology.toml` wrote SOD as `2 O2_r_minus -> H2O2 + O2`, but the water config
names superoxide `O2_minus`. SOD never acted. This is fixed in radiopyo. This folder's copy still has the old
label: replace `O2_r_minus` with `O2_minus` there. In this case the effect is small (2.38 vs 2.40 µM).
