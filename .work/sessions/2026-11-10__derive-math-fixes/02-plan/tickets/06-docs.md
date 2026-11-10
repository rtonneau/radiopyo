# Ticket 06: docs

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] README: force_positive gone; clipped concentrations and the depletion switch described where the negative products / enzymatic blocking are documented; catalase example consistent with the bundled value.
- [ ] CLAUDE.md solver and driver sections describe clipping, the switch and the smooth Jacobian; force_positive gone.
- [ ] CHANGELOG fragment for the session: Removed (force_positive), Changed (clipping, switch, catalase k), Added (poly-k warning).
- [ ] ADR 0003 gets a status line pointing to ADR 0005 for the piecewise part.
- [ ] `grep -rn force_positive radiopyo tests README.md CLAUDE.md` is empty.

**Files to Touch:**
- `README.md`
- `CLAUDE.md`
- `.work/changelog/2026-11-10__derive-math-fixes.md`
- `.work/adr/0003-analytic-jacobian-default-piecewise.md`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo && grep -rn force_positive radiopyo tests README.md CLAUDE.md
```

Expected:
All tests pass, ruff and mypy clean, and the grep prints nothing.

**Notes:**

Follow the existing CHANGELOG fragment format in .work/changelog/. Also update the "Current state" test count in CLAUDE.md.
