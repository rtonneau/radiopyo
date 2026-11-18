# Ticket 02: disabled-reactions

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] A config-level function removes every reaction with a zero scalar rate constant: `k_value` (k_reaction, constant, enzymatic), Arrhenius `A`, `k_forward`/`k_backward`/`k_equilibrium` of acid_base and equilibrium (one zero removes the whole pair); one warning per removed reaction naming its section and reaction string; polynomial (list) constants are never treated as zero.
- [ ] Applied in `SimEnv.from_dict` and in `UnitCell.from_config_dict`, whose stored `config` is the cleaned one; idempotent (no second warning); the input dict is not mutated.
- [ ] An `[includes]` override setting `k_value = 0` removes the included reaction (TOML test); no numpy divide-by-zero warning is emitted.
- [ ] Species only used by a removed reaction are absent from the env; global reaction index order unchanged for the rest.
- [ ] Rust parity unaffected (tests/rust pass or are skipped); full pytest, ruff, mypy pass.

**Files to Touch:**
- `radiopyo/parser/utils.py`
- `radiopyo/simulation/environment_build.py`
- `radiopyo/simulation/unit_cell.py`
- `tests/parser/test_disabled_reactions.py`
- `.work/adr/0011-zero-k-disables-reaction.md` (only if the implementation differs)

**Verification Step:**

Run:
```bash
$PY -W error::RuntimeWarning -m pytest tests/parser/test_disabled_reactions.py -q && $PY -m pytest -q && $PY -m ruff check radiopyo && $PY -m mypy radiopyo
```

Expected:
All tests pass; ruff "All checks passed!"; mypy "Success: no issues found".

**Notes:**

Section names: `_REACTION_LISTS` in parser/utils.py. `radiolytic` (a dict of G-values) is out of scope. A zero is `== 0` on an int/float scalar. Keep `ConfigDict` typing (TypedDicts in utils/sim_types.py).
