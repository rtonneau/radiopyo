# Ticket 09: chemistry-cleanup

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] One helper populates reactants/products from a `GenericReaction` for all reaction classes (replaces the 4 copies).
- [ ] `BiReaction` owns two `KReaction`s with their k values set; no parallel `_k_forward`/`_k_backward` storage; `AcidBase.get_couple` deduplicated; identical H2O branches in `k_forward`/`k_backward` merged.
- [ ] `KReactionCollection` lookups stay correct after `update_label`, `sort_by_label`, `sort_by_reaction_index`.
- [ ] Empty initial concentration uses `None` instead of the `-1`/`-inf` sentinel; `BasicSpecies.__eq__` returns `NotImplemented` for foreign types and has a `__hash__`; empty-list `sort`/`set_order` no longer hit an unbound variable; `ListSpecies.copy` keeps indices.
- [ ] `make_hash`/`as_label` checks shared instead of copied; `RadiolyticReaction.label` consistent with other classes while keeping the old call form working.
- [ ] Docstrings no longer mention `RawSpecies`/`ChemSpecies`; the n·k convention (ADR 0001) is documented in `KReaction`.
- [ ] Golden tests pass.

**Files to Touch:**
- `radiopyo/chemistry/*.py`
- `radiopyo/chemistry/utils/factories.py`
- `tests/chemistry/test_*.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass, zero lint/type errors.

**Notes:**

`MichaelisMenten.__repr__` is used for hashing: fixing its missing `)` must not change any stored key (check HDF5/identifier use). Keep `RadiolyticReaction.label()` callable if users call it (e.g. via a callable-str or keep method and add a property alias); choose the least surprising option and state it.
