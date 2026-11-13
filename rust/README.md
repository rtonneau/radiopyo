# radiopyo-rs

Optional Rust backend of [radiopyo](https://github.com/rtonneau/radiopyo):
an in-process PyO3 binding of [radkin](https://github.com/rtonneau/radkin),
the Rust port of radiopyo. Select it with `backend="rust"` on
`UnitCell.run` / `run_till_max_dose` / `SimChunksIter.run`, or run several
cells in parallel with `radiopyo.rust.run_batch`.

Install from a release wheel (`pip install radiopyo[rust]`) or from source
(needs a Rust toolchain):

```bash
pip install ./rust          # or: cd rust && maturin develop --release
```

radkin is pinned by git revision in `Cargo.toml`. To build against a local
radkin checkout, add a git-ignored `rust/.cargo/config.toml`:

```toml
[patch."https://github.com/rtonneau/radkin"]
radkin = { path = "C:/DEV/RUST/radkin" }
```
