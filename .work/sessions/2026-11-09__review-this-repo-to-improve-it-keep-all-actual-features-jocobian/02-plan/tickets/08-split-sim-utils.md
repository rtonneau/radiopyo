# Ticket 08: split-sim-utils

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `simulation/utils.py` split into `units.py` (`*Convert`), `collections.py` (`ResCollection`), `persistence.py` (`H5Manager`, `SimIdentifier`, `SimCollection`); `simulation/utils.py` re-exports all of them and its `__all__` lists them.
- [ ] `ODEResult.to_hdf` no longer duplicates `push_ODEResult` (one delegates to the other).
- [ ] Crash fixes: `ResCollection.sort` sorts by `_order`; `SimCollection.__len__` returns the real count; `read_ODEResult` checks `json_config` before use.
- [ ] `SimIdentifier` uses `typing.get_type_hints` so subclasses with `from __future__ import annotations` work; serialized names are byte-identical to before.
- [ ] The HDF5 file written by the old code (ticket 01) loads and matches; golden tests pass.

**Files to Touch:**
- `radiopyo/simulation/utils.py`
- `radiopyo/simulation/units.py`
- `radiopyo/simulation/collections.py`
- `radiopyo/simulation/persistence.py`
- `radiopyo/simulation/ode_result.py`
- `tests/simulation/test_persistence.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass, zero lint/type errors.

**Notes:**

`parse_dataset_name` uses `.strip(".sim")`: only change it if every name written by the old code still parses identically (test it); otherwise record it. Avoid naming the module `collections.py` if it shadows the stdlib in any import path; pick `result_collection.py` in that case.
