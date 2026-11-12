# Ticket 02: showcase-notebook

**Status:** ✅ Done

## Local Test Result

Full map first (outside the notebook, same arguments): `python .scratch/tests/2026-11-12__fenton-fix-showcase-notebook/run_map.py` -> `done 168.4 (12, 12)`, cache `docs/issues/fenton/.map_cache/fenton_map_4ace6041a0a8.npz` (144 points, 0 NaN).

Verification Step, run in `docs/issues/fenton`:
`python -m jupyter nbconvert --to notebook --execute Fenton_Fixed_Showcase.ipynb --output .scratch_check.ipynb` -> exit 0 ("Writing ... bytes to .scratch_check.ipynb"; log in `.scratch/tests/2026-11-12__fenton-fix-showcase-notebook/evidence/nbconvert.log`). Map loaded from cache ("Map ready in 0.0 s"). All code cells have outputs, no error output, all chunks `success=True`, 0 failed scan runs. 1200 s ROOH: UHDR 2.402, conv 6.612, UHDR_Freeze 3.01 umol/L; scan values match ANALYSIS.md section 1 (e.g. 0.13 Hz: 1.759 / 3.010). The executed `.scratch_check.ipynb` was then moved over `Fenton_Fixed_Showcase.ipynb` (outputs saved); no stray output left in the fenton folder.

## Review Notes

- Sections present: recap, 3 runs (UHDR/conv/UHDR_Freeze) with ROOH/H2O2/Fe3+ plots and solver stats table, frequency scan 0.10-0.16 Hz Fe=0 vs 2.2e-8, normalized ROOH map at max dose (jet, log axes, white NaN, text box, colorbar 0-1), ROOH-at-1200 s map with the drift note.
- Executed top to bottom with outputs saved, no solver failure.
- Folder TOMLs unchanged in content (git shows only the pre-existing line-ending diff on the Configuration file). SOD `O2_r_minus` label mentioned in markdown.
- Changed `docs/issues/fenton/fenton_map.py` (ticket 01 file): added `load_config(fe)` (exported) and made `_get_cell` use it. See Blockers.
- Cache file is `docs/issues/fenton/.map_cache/fenton_map_4ace6041a0a8.npz` (not `fenton_map_cache.npz` as named in the ticket); it is untracked and must be committed. Map grid: n=12, full size.
- Scan uses the beam's mean dose rate 0.4 Gy/s (TOML `dose_rate` is the mean; peak = mean / (on_time x f)); ANALYSIS.md calls it "peak", the notebook says "mean".

## Blockers / Challenges

- Bug in ticket 01 helper: `Configuration_Bordet_Emballement.toml` includes `radiopyo/config_Labarbe.toml`, hence the PACKAGED biology (catalase 8e-8, GSH, Fe 8.9e-7), not this folder's `reactions_biology.toml`; ROOH(UHDR_Freeze, 1200 s) was ~1e-3 umol/L instead of 3.01. A first full map run on that config was killed and discarded. Fix: `load_config` merges water radiolysis + folder `reactions_biology.toml` (with its `H2O = 55` line dropped in memory, since it clashes with the water config) + the Configuration's beams, via `combine_config`. Reproduces the ANALYSIS.md values. The ticket 01 tests/behaviour of `fenton_map` for the wrong config are superseded; the parent may want to re-check any ticket 01 evidence that used the old cell.
- The notebook must import `fenton_map` before `radiopyo` (it puts the repo checkout on `sys.path`; an installed radiopyo lacked `radiopyo.physics`).
- tqdm notebook widget crashed when disabled; the notebook swaps `tqdm.notebook.tqdm` for the std tqdm.
- conv run takes ~5 min, so the notebook runs in ~6 min (map from cache).

## Commits

- 5cbf6c4 docs: add Fenton fix showcase notebook with ROOH maps (ticket 02)

## Time Spent

21m (ticket-start.js to ticket-complete.js)
