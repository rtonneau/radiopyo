# Implementation Plan

**Session:** clone docs/issues/fenton/Fenton_Fixed_Showcase.ipynb and make it run simulation with rust backend
**Date:** 2026-11-13T14:27:17.767Z
**Estimated effort:** 1-2 hours

## Strategy

First give `fenton_map.compute_map` a Rust path built on `rp.rust.run_batch`. Then clone the notebook, switch every run to `backend="rust"`, execute it and commit it with its outputs.

## Tickets Overview

- **Ticket 1:** `fenton_map.compute_map(backend="rust")` via `run_batch`, with helpers shared with the Python path and an unchanged Python cache key.
- **Ticket 2:** `Fenton_Fixed_Showcase_Rust.ipynb`, a pure Rust clone of the showcase, executed and committed with its outputs.

## Sequencing Rationale

The notebook's map section calls `compute_map(backend="rust")`, so the helper must come first.

## Risks & Mitigation

- **Risk:** the Rust map differs from the Python map because of a wrong index extraction or feasibility mask → **Mitigation:** compare it with the cached Python map (`.map_cache/fenton_map_4ace6041a0a8.npz`): same NaN mask, values within solver tolerance.
- **Risk:** refactoring `compute_point` changes the Python cache key or its values → **Mitigation:** the key params dict stays identical for `backend="python"`. Check that the existing cache file still loads (same hash).

## Assumptions

- The `radiopyo-rs` extension is installed in `radiobio_env` (checked: `rp.rust.is_available()` is True).
- nbconvert/ipykernel are available in `radiobio_env`. If not, the notebook is executed with the radiobio kernel by another available means.
