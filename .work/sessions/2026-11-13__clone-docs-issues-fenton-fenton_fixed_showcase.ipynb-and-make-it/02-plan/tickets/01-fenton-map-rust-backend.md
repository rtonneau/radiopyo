# Ticket 01: fenton-map-rust-backend

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `compute_map(..., backend="rust")` returns the same keys and shapes as the Python path, computed with `rp.rust.run_batch(mode="chunked")` and the same max step sizes.
- [ ] The feasibility rules (pulse period <= pulse, max dose not reached within `period`) and the ROOH extraction (at max dose, at the end) are shared helpers used by both paths. A failed Rust point is NaN.
- [ ] For `backend="python"` the cache key is unchanged (the existing `.map_cache` file is still found). A Rust map has its own key.
- [ ] An unknown backend raises `ValueError`.
- [ ] The Rust map matches the cached Python map: same NaN mask, values within rtol 1e-3.

**Files to Touch:**
- `docs/issues/fenton/fenton_map.py`

**Verification Step:**

Run:
```bash
cd docs/issues/fenton && $PY -c "import numpy as np, fenton_map as f; p=f.compute_map(); r=f.compute_map(backend='rust', cache=None); m=np.isnan(p['rooh_max_dose']); print((m==np.isnan(r['rooh_max_dose'])).all(), np.nanmax(abs(r['rooh_max_dose']/p['rooh_max_dose']-1)), np.nanmax(abs(r['rooh_end']/p['rooh_end']-1)))"
```

Expected:
`True` followed by two relative differences well below 1e-3. The Python map loads from the existing cache (no recompute).

**Notes:**

`$PY` = `C:/Users/rtonneau/AppData/Local/miniconda3/envs/radiobio_env/python.exe`. Keep the module banner layout. Build the cells with `copy.deepcopy` of the per-Fe cell and then `set_beam(beam, LET=LET)`. Pass `return_exceptions=True` so that one bad point does not kill the map.
