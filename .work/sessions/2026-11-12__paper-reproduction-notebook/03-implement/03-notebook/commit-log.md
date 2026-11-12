# Ticket 03: notebook

**Status:** ✅ Done

## Local Test Result

- Verification step: `PYTHONPATH=. python -m jupyter nbconvert --to notebook --execute examples/tonneau_2025_pmb/tonneau_2025_pmb.ipynb --output-dir .scratch/tests/2026-11-12__paper-reproduction-notebook` exits 0.
  - The executed copy holds 7 figures: figs 4, 7, 10, 5, 8, 6, 9.
  - It also holds the GF table and the 17-row comparison table. The figures were extracted with `03-extract.py` as `03-fig1..7.png`.
- Each figure was inspected against the paper pages, and they match qualitatively:
  - ΔO2 and net-H2O2 hot spots at 1e3–1e4 Gy/s and low pO2;
  - the two-step O2 recovery;
  - the bell-shaped gross H2O2 production at UHDR;
  - the pH dip.
- Comparison table: radiopyo is mostly 3–8 % below the values quoted in the paper text, and the pH minimum is about 1 % above. The largest gap is the CDR ΔO2 range: 0.57–1.06 µM against 0.8–1.2 µM.
- `ruff check examples/tonneau_2025_pmb/ --ignore RUF001,T201`: clean.
- Outputs stripped with `nbconvert --clear-output --inplace`; 0 outputs remain (21 kB).

## Review Notes

Checked against the acceptance criteria:
- Scenario: the notebook explains tables 3/5 and the GF gap (ADR 0009). The GF table shows the code's fit against the paper's: +4 % at 0.1 Gy/s, equal near 1 Gy/s, −28 % at 1e7 Gy/s.
- `RECOMPUTE = False` loads `tests/golden/data/paper_tonneau_2025.npz`.
  - `True` re-runs `pr.run_grid(workers=WORKERS)` and displays the max relative difference to the golden file for each key.
  - The True branch was not executed end to end: it takes about 10 min, and the same `run_grid` produced the golden file in ticket 02.
- Figures:
  - figs 4, 7 and 10 as 2×2 panels (CDR/UHDR × 1 %/5 %);
  - fig 5 in %pO2 and µM (top axis and colorbar secondary axis);
  - fig 8 in #/100 eV;
  - figs 6 and 9 as bars with a net marker.
- The informational comparison table makes no assertions.
- The setup cell finds the repository root from the CWD (pyproject.toml plus `radiopyo/`). It adds the repository to `sys.path` only when radiopyo is not installed (`importlib.util.find_spec`).
- After review:
  - the "Reading the gaps" markdown was rewritten to match the data;
  - ruff findings were fixed (import order, unused noqa, an unused import).
  - The remaining RUF001 (× and – in display strings) and T201 (notebook `print`) are deliberate. The project's lint scope is `radiopyo/` only.

## Blockers / Challenges

- My first draft of the gap explanation contradicted the computed numbers. It was rewritten after comparing against the table.

## Commits

- 1ad97fd feat: Tonneau 2025 paper reproduction notebook (ticket 03)

## Time Spent

6m (ticket-start.js to ticket-complete.js)
