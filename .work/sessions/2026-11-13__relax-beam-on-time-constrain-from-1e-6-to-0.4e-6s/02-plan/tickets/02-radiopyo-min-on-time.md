# Ticket 02: radiopyo-min-on-time

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `MIN_BEAM_ON_TIME == 0.4e-6` in `radiopyo/physics/beam.py`, UnitCell and `_RADKIN_CONSTANTS`; docstrings (unit_cell.py, exceptions.py, ode_result.py), README and CHANGELOG updated.
- [ ] `tests/physics/test_min_beam_on_time.py` boundary at 0.4 µs, plus a test that a 0.5 µs pulse runs and a 0.3 µs pulse raises ShortBeamOnError.
- [ ] `.work/adr/0008-*.md` has a dated amendment (0.4e-6 s, still 4x T_ZERO, because real UHDR machines deliver 0.4–1 µs pulses).
- [ ] `rust/Cargo.toml` radkin rev = ticket 1 SHA; extension rebuilt; pytest (incl. tests/rust), ruff, mypy green.

**Files to Touch:**
- `radiopyo/physics/beam.py`
- `radiopyo/simulation/unit_cell.py`
- `radiopyo/physics/exceptions.py`
- `radiopyo/simulation/ode_result.py`
- `radiopyo/rust.py`
- `rust/Cargo.toml`
- `tests/physics/test_min_beam_on_time.py`
- `.work/adr/0008-short-on-spans-added-analytically.md`
- `README.md`
- `CHANGELOG.md`

**Verification Step:**

Run:
```bash
conda run -n radiobio_env python -m pytest -q && conda run -n radiobio_env ruff check radiopyo && conda run -n radiobio_env mypy radiopyo
```

Expected:
all pass, zero lint/type errors

**Notes:**

Env name on this machine: see memory (radiobio_env under AppData miniconda). Rebuild the extension in `rust/` with maturin develop --release. Do not re-record golden data; report any golden change.
