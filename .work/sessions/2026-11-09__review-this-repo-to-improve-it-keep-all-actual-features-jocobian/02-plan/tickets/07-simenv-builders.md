# Ticket 07: simenv-builders

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] The `*_from_dict` reaction readers become one table-driven reader; `ReactionsList` global index order (acid/base, equilibrium, k, enzymatic, radiolytic) is unchanged.
- [ ] `reactions_involving(species, role=..., kinds=...)` covers all `*_involving*` queries; old method names remain as thin aliases with identical return values.
- [ ] Unused `SimEnv` fields/flags (`pKw`, `pKH2O`, `update_temperature`/`update_flag_positive` writes nobody reads) removed or wired; `pH=0` no longer treated as unset only if it changes no golden (else record for Phase C).
- [ ] Crash fixes: `get_radiolytic_reaction_from` with no radiolytic reactions; `from_dict` missing `concentrations`/`parameters` raises a clear error; logger level restored on exception.
- [ ] Golden tests pass.

**Files to Touch:**
- `radiopyo/simulation/environment_build.py`
- `tests/simulation/test_env_queries.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass, zero lint/type errors.

**Notes:**

Add a test that every old query name returns the same reactions as before on the Labarbe config (compare by label).
