# Ticket 04: docs-changelog

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] README and `docs/configuration.md` describe: disabling a reaction with a zero rate constant (also for included reactions; whole pair), H2O cannot be set in `[concentrations.*]`, the old acid/base syntax message, and `check_chemistry()` with its four checks and log-once behaviour.
- [ ] Changelog entry (Unreleased, or a `.work/changelog/` fragment per project convention) listing the new warnings, disabled reactions and the chemistry check.

**Files to Touch:**
- `README.md`
- `docs/configuration.md`
- `CHANGELOG.md` or `.work/changelog/<fragment>.md`

**Verification Step:**

Run:
```bash
git diff --stat && $PY -m pytest -q
```

Expected:
Only doc/changelog files changed in this ticket; all tests pass.

**Notes:**

Follow the existing style of `.work/changelog/2026-11-18__tuto-update-windows.md` for the fragment.
