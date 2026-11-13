# Ticket 06: notebook-and-docs

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `examples/notebook_rust_backend.ipynb` (kernel "Python (radiobio)") runs top to bottom: install check, a Labarbe run with both backends and an overlay plot, a chunked pulsed run, a `run_batch` LET sweep, a timing comparison table; executed with outputs saved.
- [ ] README: "Rust backend (optional)" section (install from release wheel or git subdirectory, `backend="rust"`, `run_batch`, limits: LSODA only, parity guard, config-built UnitCells only); `docs/configuration.md` running section updated.
- [ ] CLAUDE.md architecture gets a short Rust backend entry (crate location, pinned radkin rev, parity guard, tests skipped without the extension); CHANGELOG fragment in `.work/changelog/`.

**Files to Touch:**
- `examples/notebook_rust_backend.ipynb`
- `README.md`
- `docs/configuration.md`
- `CLAUDE.md`
- `.work/changelog/2026-11-13__ship-rust-radkin-backend.md`

**Verification Step:**

Run:
```bash
<radiobio_env python> -m jupyter nbconvert --to notebook --execute --inplace examples/notebook_rust_backend.ipynb
```

Expected:
the notebook executes without error and saves its outputs.

**Notes:**

Follow the existing notebook style (`examples/notebook_basic_usage.ipynb`). Keep prose plain; the CHANGELOG fragment format follows the existing ones in `.work/changelog/`.
