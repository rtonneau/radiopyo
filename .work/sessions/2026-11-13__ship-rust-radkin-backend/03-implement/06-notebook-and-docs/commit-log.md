# Ticket 06: notebook-and-docs

**Status:** ✅ Done

## Local Test Result

`jupyter nbconvert --to notebook --execute --inplace examples/notebook_rust_backend.ipynb`: rc 0, outputs saved. Outputs:
- Labarbe constant beam: Python 0.13 s, Rust 0.05 s; max relative difference 1.8e-8 over the 14 species above 1e-6 µmol/L.
- Pulsed chunked run: Python 1.61 s, Rust 0.05 s.
- run_batch of 5 LETs: 0.18 s.

## Review Notes

Base 739b486.
- Notebook: covers the install check, a Labarbe run with both backends plus an overlay plot, a chunked pulsed run, a run_batch LET sweep and a timing table.
  - Kernelspec name is `python3`, with display name "Python (radiobio)": that kernel name does not exist on this machine, and `python3` resolves to the env.
- README: new "Rust backend (optional)" subsection (install from a wheel, `./rust` or a git subdirectory; backend=; run_batch; limits) and a documentation link.
- docs/configuration.md: new Rust item under "Ways of Running Simulation".
- CLAUDE.md: architecture item 8, naming the GitHub repo rather than a machine path, per the global no-hardcoded-paths rule.
- Changelog fragment, minor bump.

Change after review: the first relative-difference printout (3.2e-2) came from trace radicals around 1e-11 µmol/L, at the solver atol. The notebook now compares species above 1e-6 µmol/L and says why.

## Blockers / Challenges

No `Python (radiobio)` kernelspec is registered on this machine; I used `python3` from the env.

## Commits

- 9b19eb2 docs: Rust backend notebook, README, configuration and CLAUDE.md (ticket 06)

## Time Spent

6m (ticket-start.js to ticket-complete.js)
