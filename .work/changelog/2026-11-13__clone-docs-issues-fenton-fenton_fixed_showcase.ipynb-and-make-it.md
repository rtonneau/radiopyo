---
bump: patch
floor: patch
---
### Added
- `docs/issues/fenton/Fenton_Fixed_Showcase_Rust.ipynb`: the Fenton runaway showcase run entirely with the Rust backend (`backend="rust"`), with the 12x12 ROOH map computed live in about 15 s.
- `docs/issues/fenton/fenton_map.compute_map(backend="rust")` solves every map point in one `rp.rust.run_batch` call. The Python path and its cache are unchanged.
