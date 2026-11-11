# Ticket 03: runtime-short-on-spans-docs

**Model:** heavy
**Effort:** high

**Acceptance Criteria:**
- [ ] An ON span shorter than `MIN_BEAM_ON_TIME` (last pulse cut by `max_dose`, micro-pulse cut by its train end, pulse clipped by `t_span`) in `prepare_chunked_run`, `prepare_single_pulse` and `run()` is not integrated: the radiolytic yield over the whole span (`compute_initial_cc_before(t0=span length)`, at the chunk dose rate) is added analytically, one result point is stored at the span end.
- [ ] A `logger.warning` names the span; `ODEResult.short_on_spans` (list of `(start, stop)`) records it, is concatenated by `__add__`, and is empty for results loaded from HDF5 (no file-format change).
- [ ] The following chunk starts from that point; the total delivered radiolytic production matches the expected dose (test with a pulsed beam whose max_dose leaves a short last pulse).
- [ ] `reconstruct_reactions` / `to_pandas` work on a result containing a short span.
- [ ] ON spans >= 1e-6 s behave exactly as before (golden tests unchanged).
- [ ] Docs: glossary term *Short ON span*; ADR 0008 (short spans approximated by a burst over their length instead of raising); CHANGELOG fragment (Added/Changed, Known issue replaced); README and CLAUDE.md (the "pulse shorter than T_ZERO raises ValueError, follow-up pending" text replaced).

**Files to Touch:**
- `radiopyo/simulation/unit_cell.py`
- `radiopyo/simulation/ode_result.py`
- `tests/simulation/`
- `.work/GLOSSARY.md`, `.work/adr/0008-*.md`, `.work/changelog/`, `README.md`, `CLAUDE.md`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
Only the two known golden time-grid failures; ruff and mypy report zero errors.

**Notes:**

Use `domain-doc.js glossary` / `domain-doc.js adr <slug> --title` for the glossary and ADR. A short span right at the end of a run (t_span clipping) still yields its end point. In `run()`, a t_span shorter than MIN_BEAM_ON_TIME with the beam ON is a short span too.
