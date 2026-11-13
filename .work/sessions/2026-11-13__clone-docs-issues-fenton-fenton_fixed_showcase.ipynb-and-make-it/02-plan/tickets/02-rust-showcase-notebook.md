# Ticket 02: rust-showcase-notebook

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `docs/issues/fenton/Fenton_Fixed_Showcase_Rust.ipynb` has the same sections, cells and figures as the original, with every simulation run with `backend="rust"` and no Python-backend run.
- [ ] Section 4 calls `compute_map(n=12, period=PERIOD, backend="rust", cache=None)` and prints its time. The markdown says that the Rust backend is used, and the cache and process-pool wording is updated.
- [ ] The notebook executes top to bottom without errors in `radiobio_env` and is committed with its outputs.
- [ ] The original `Fenton_Fixed_Showcase.ipynb` is not modified by this session.

**Files to Touch:**
- `docs/issues/fenton/Fenton_Fixed_Showcase_Rust.ipynb`

**Verification Step:**

Run:
```bash
cd docs/issues/fenton && $PY -m jupyter nbconvert --to notebook --execute --inplace Fenton_Fixed_Showcase_Rust.ipynb && $PY -m pytest -q ../../../tests -x -p no:cacheprovider
```

Expected:
nbconvert writes the notebook without an error. The section 2 table shows every chunk OK and ROOH(1200 s) about 2.40 µmol/L for UHDR. Section 3 prints "Runs with a failed chunk: 0". pytest passes.

**Notes:**

Build the clone with a small Python script that edits the JSON (do not hand-edit the source). Keep the original's styling and its non-ASCII characters (µ, ·) correct in UTF-8.
