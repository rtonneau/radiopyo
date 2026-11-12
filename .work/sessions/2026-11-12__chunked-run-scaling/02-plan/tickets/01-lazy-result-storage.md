# Ticket 01: lazy-result-storage

**Model:** medium
**Effort:** high

**Acceptance Criteria:**
- [ ] `ODEResult` keeps the per-step `t` / `y` arrays in lists; `t` and `y` are properties returning one concatenated ndarray, cached until the next `append_sol`; assignment (`res.t = ...`, `res.y = ...`) still works (used by `_default_init`, `__add__`, loaders).
- [ ] `append_sol` does not copy earlier steps (no np.append / np.hstack of the history).
- [ ] `SimChunksIter.__next__` takes the next y0 from the last appended step's own array, never from the concatenated `self.sol.y`.
- [ ] `__add__`, `final_cc`, `__getitem__`, `to_pandas`, `reconstruct_reactions`, `_stacked_rates`, HDF5 push/read keep their results (bit-identical).
- [ ] New test `tests/simulation/test_result_storage.py`: appending a step leaves the earlier step arrays untouched (identity / no reallocation), cache invalidated on append, assignment round trip, `__add__` operands unchanged.
- [ ] Full pytest, ruff and mypy pass.

**Files to Touch:**
- `radiopyo/simulation/ode_result.py`
- `radiopyo/simulation/unit_cell.py`
- `tests/simulation/test_result_storage.py`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All tests pass (golden included), ruff "All checks passed!", mypy "Success: no issues found".

**Notes:**

`SimChunksIter` yields `self.sol` on every chunk: a user reading `.y` inside the loop triggers one concatenation per read (documented, not avoided). Keep the empty-state shape behaviour of `_default_init` (`y` empty 1-D array before the first step) where callers rely on `len(self.y) == 0`.
