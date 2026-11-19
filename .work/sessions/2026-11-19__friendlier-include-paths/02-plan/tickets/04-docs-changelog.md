# Ticket 04: docs-changelog

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] `docs/configuration.md` includes section: table and list forms, where files are looked for (incl. CWD warning), `rp.bundled_configs()`, example error messages, Windows-path tip (single quotes / forward slashes).
- [ ] CHANGELOG Unreleased: Added (list form, `bundled_configs`, `ConfigError` hierarchy, file origin in errors), Changed (unknown include key is an error; build errors are `ConfigError` instead of KeyError/TypeError/...; CWD fallback warns).
- [ ] CLAUDE.md parsing paragraph mentions `ConfigError`, origins and the list form.

**Files to Touch:**
- `docs/configuration.md`
- `CHANGELOG.md`
- `CLAUDE.md`

**Verification Step:**

Run:
```bash
C:/Users/rtonneau/AppData/Local/miniconda3/envs/radiobio_env/python.exe -m pytest -q
```

Expected:
Full suite passes.

**Notes:**

Match the existing doc tone and structure; link ADR 0012 in the CHANGELOG entry.
