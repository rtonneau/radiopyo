# Ticket 02: ci-workflow

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] .github/workflows/tests.yml named "Tests", on push and pull_request to main, ubuntu-latest, Python 3.11.
- [ ] Installs `.[dev]` and ruff, runs `ruff check radiopyo`, `mypy radiopyo`, `python -m pytest` with PYTHONHASHSEED=0.
- [ ] The same three commands pass locally.

**Files to Touch:**
- `.github/workflows/tests.yml`

**Verification Step:**

Run:
```bash
python -c "import yaml; yaml.safe_load(open('.github/workflows/tests.yml')); print('yaml ok')" && ruff check radiopyo && mypy radiopyo && PYTHONHASHSEED=0 python -m pytest -q
```

Expected:
"yaml ok", ruff and mypy with no errors, all tests pass.

**Notes:**

The full pytest run takes a few minutes.
