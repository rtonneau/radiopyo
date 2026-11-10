# Ticket 03: packaging-hygiene

**Model:** fast
**Effort:** medium

**Acceptance Criteria:**
- [ ] `mypy`, `pytest`, `pandas-stubs` move to `[project.optional-dependencies] dev`; `tqdm` to an optional extra; runtime deps unchanged otherwise.
- [ ] `dist/`, `.coverage`, `tests/output/` are untracked (`git rm --cached`) and listed in `.gitignore`; history untouched.
- [ ] `radiopyo/data/IBA/` is excluded from the wheel and sdist (files stay in the repo); `python -m build` produces a wheel without it.
- [ ] ruff settings move under `[tool.ruff.lint]` (no deprecation warning).
- [ ] Top-level `__all__` has no duplicates and exports `MacroPulsedBeam`; `rp.MacroPulsedBeam` works; `utils/logging.py` and `physics/beam.py` `__all__` list their public names.
- [ ] Logging: `import radiopyo` no longer disables pre-existing loggers (`disable_existing_loggers = false`); INFO console output from radiopyo unchanged; one config source (pyproject and `default_log_config.toml` no longer diverge, the dangling `InfoOnly` filter is removed or implemented); importing next to an unrelated `pyproject.toml` without `[tool.logging]` falls back instead of raising.
- [ ] `project.urls` points to GitHub.

**Files to Touch:**
- `pyproject.toml`
- `.gitignore`
- `radiopyo/__init__.py`
- `radiopyo/utils/logging.py`
- `radiopyo/utils/default_log_config.toml`
- `radiopyo/physics/beam.py` (`__all__` only)
- `tests/test_logging.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && python -m build && python -c "import zipfile,glob;n=zipfile.ZipFile(sorted(glob.glob('dist/*.whl'))[-1]).namelist();print(any('IBA' in x for x in n))"
```

Expected:
Tests pass, build succeeds, prints `False`.

**Notes:**

Add a test creating a logger before `import radiopyo` (subprocess) and asserting it is still enabled afterwards.
