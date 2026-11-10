# Ticket 14: fix-include-merge

**Model:** heavy
**Effort:** xhigh

**Acceptance Criteria:**
- [ ] `combine_config` becomes one generic merge keyed on an order-insensitive reaction key, without building reaction objects and without importing `radiopyo.simulation` from `radiopyo.parser`.
- [ ] Overriding an enzymatic reaction keeps k (no doubling); overriding an acid/base with only `k_forward` keeps `k_backward`; constant-reaction overrides replace instead of duplicating; `B + A -> C` overrides `A + B -> C`.
- [ ] Precedence: later includes override earlier ones, the including file overrides all includes.
- [ ] Relative include paths resolve against the including file, falling back to the CWD; only the exact `radiopyo/` prefix maps to package data.
- [ ] Bundled configs produce identical environments (goldens unchanged); before/after report for a config hitting each fixed case shown to the user.
- [ ] README include section updated.

**Files to Touch:**
- `radiopyo/parser/utils.py`
- `radiopyo/parser/file_parser.py`
- `radiopyo/parser/toml.py`
- `tests/parser/test_includes.py`
- `README.md`
- `CHANGELOG.md`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass; golden tests unchanged.

**Notes:**

Also fix `TOMLFileParser` overwriting the default parameters it just built (toml.py:271-272) if it is part of the same merge path. Update the ticket-02 include tests that pinned the old behaviour.
