# Ticket 08: split-sim-utils

**Status:** ✅ Done

## Local Test Result

`python -m pytest -q && ruff check radiopyo && mypy radiopyo`: 391 passed (384 baseline + 7 new), ruff clean, mypy clean (39 files). Golden tests unchanged and passing (old HDF5 loads).

## Review Notes

- Parent: base 5c1bec6cdd51b46164be8467b598392d3ec98dae; pre-existing changes: gps session state only.

- Split into `units.py`, `result_collection.py` (chosen over `collections.py` to avoid stdlib shadowing), `persistence.py`; `utils.py` is a facade re-exporting everything (incl. SpeciesConvert, CO2/N2Convert, exception classes) with a full `__all__`. `ode_result` now lazily imports from `persistence`.
- `ODEResult.to_hdf` delegates to `H5Manager.push_ODEResult` (which now also accepts a plain str name). Behaviour note: to_hdf no longer writes the extra solver attrs (sol, nfev, ...) and now also writes the json config like push_ODEResult; it had no callers.
- Crash fixes: `ResCollection.sort` uses `_order`; `SimCollection.__len__` uses `len(_content)`; `read_ODEResult` checks `json_config` before use.
- `SimIdentifier` uses `typing.get_type_hints` (new `field_types()`); byte-identical names verified against strings computed from base code (test_persistence.py), plus a future-annotations subclass test in a separate module.
- `parse_dataset_name` `.strip(".sim")` left unchanged (strips characters, not a suffix; could mangle names ending in s/i/m/. chars) - not fixed, recorded only.
- Added tests: tests/simulation/test_persistence.py, test_persistence_future_annotations.py.

- Parent acceptance: independent reviewer (sonnet/medium): no findings; parent re-ran pytest (391 passed incl. golden HDF5 written by old code), ruff 0, mypy 0. to_hdf now delegates to push_ODEResult (old version wrote OdeSolution attrs, almost certainly failed; extra solver attrs no longer written, never read). .strip('.sim') recorded for ticket 18 seeds. CLAUDE.md item 7 stale (ticket 18).

## Blockers / Challenges

None

## Commits

- aefda6e refactor: split simulation utils into units, result_collection, persistence (ticket 08)

## Time Spent

6m (ticket-start.js to ticket-complete.js)
