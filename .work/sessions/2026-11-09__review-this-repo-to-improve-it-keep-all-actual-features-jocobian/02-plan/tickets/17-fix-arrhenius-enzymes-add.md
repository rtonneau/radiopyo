# Ticket 17: fix-arrhenius-enzymes-add

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `[reactions.arrhenius]` accepts `A` and `Ea` (kJ/mol); k = A·exp(−Ea/(R·T)) with T in K; a test checks a hand-computed value; README documents it.
- [ ] Enzyme names with capitals (`-- SOD >>`, `-- GPx >>`) parse to the enzyme species, not a `BasicSpecies` holding a dict.
- [ ] `ODEResult.__add__` returns a new object and leaves both operands unchanged.
- [ ] Goldens unchanged.

**Files to Touch:**
- `radiopyo/chemistry/k_reactions.py`
- `radiopyo/simulation/environment_build.py`
- `radiopyo/parser/grammar/reaction_grammar.lark`
- `radiopyo/parser/reaction_parser.py`
- `radiopyo/simulation/ode_result.py`
- `tests/chemistry/test_k_reactions.py`
- `tests/parser/test_reaction_parser.py`
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

If `species_balance_evolution` really omits enzymatic/radiolytic terms (plausible finding), fix it here only with a before/after shown to the user; otherwise record it.
