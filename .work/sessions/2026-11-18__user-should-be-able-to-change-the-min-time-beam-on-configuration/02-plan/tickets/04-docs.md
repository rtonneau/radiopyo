# Ticket 04: docs

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] ADR 0008 has an amendment: the minimum is a global setting (default 0.4e-6 s), spans <= T_ZERO always analytic, radkin takes it as input.
- [ ] Glossary "Short ON span" no longer states a fixed value; README beam timing section and CLAUDE.md describe `rp.set_min_beam_on_time`.

**Files to Touch:**
- `.work/adr/0008-short-on-spans-added-analytically.md`
- `.work/GLOSSARY.md`
- `README.md`
- `CLAUDE.md`

**Verification Step:**

Run:
```bash
grep -n "set_min_beam_on_time" README.md CLAUDE.md .work/adr/0008-short-on-spans-added-analytically.md
```

Expected:
At least one match in each file.

**Notes:**

`.work/GLOSSARY.md` already has unrelated uncommitted edits in the main checkout: commit only this session's hunk.
