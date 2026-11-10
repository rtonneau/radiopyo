# Ticket 05: docs

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] README documents `use_jac` (default, methods) and the enzymatic depletion block.
- [ ] CHANGELOG Unreleased lists the Jacobian, the default change and the F9 fix; CLAUDE.md solver/driver sections updated (no more "no jacobian").

**Files to Touch:**
- `README.md`
- `CHANGELOG.md`
- `CLAUDE.md`

**Verification Step:**

Run:
```bash
python -m pytest -q && ruff check radiopyo && mypy radiopyo
```

Expected:
All pass.

**Notes:**

ADR 0003 and glossary terms already written in the grill.
