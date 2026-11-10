# Ticket 04: dead-code-lint-zero

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] Removed: `parser/ron.py`, `parser/grammar/ron_grammar.lark`, `UnitCell.from_ron`, `run_multi_doses`, `utils/e_LET.py`, `data/energy_LET/`, `radiopyo/test.py`, empty `physics/utils.py`, `make_beam_func`, `RawReader`, `_average_val_before_100ns`, `reconstruct_jacobian`, `push_result`/`read_result`, `KReaction.from_dict`/`Arrhenius.from_dict` stubs, `umol_per_joule_to_ge`, `dict_compare`, `MichaelisMentenDict`, `_set_mass` dead body, commented-out code, the stray `print` in `derive`.
- [ ] A CHANGELOG entry lists the breaking removals (`from_ron`, RON format, `run_multi_doses`, `e_LET`).
- [ ] `ruff check radiopyo` and `mypy radiopyo` report zero errors (fixing remaining errors without behaviour change).
- [ ] Golden tests still pass.

**Files to Touch:**
- `radiopyo/parser/*`
- `radiopyo/simulation/*.py`
- `radiopyo/chemistry/*.py`
- `radiopyo/physics/*.py`
- `radiopyo/utils/*.py`
- `CHANGELOG.md`
- `CLAUDE.md` (remove RON/run_multi_doses mentions)

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
Tests pass; `All checks passed!`; `Success: no issues found`.

**Notes:**

Keep `jacobian()` for now (ticket 05 handles it). If a mypy fix would change behaviour, use a narrow cast and leave a note rather than changing logic. `physics/exceptions.py` importing from `simulation` is a layering inversion: fix it here if trivial.
