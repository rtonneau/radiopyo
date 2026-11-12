# Ticket 01: docs-and-logo

**Model:** fast
**Effort:** low

**Acceptance Criteria:**
- [ ] docs/configuration.md contains, verbatim, the README sections "Configuring Simulation" (all TOML subsections), "Possible beam definitions" and "Ways of Running Simulation", with a short title/intro line.
- [ ] docs/logo.svg: valid SVG, an H2O/radical molecule with a radiation wave next to the "radiopyo" wordmark, readable on light and dark backgrounds, no external fonts/images required.
- [ ] README.md is not modified in this ticket.

**Files to Touch:**
- `docs/configuration.md`
- `docs/logo.svg`

**Verification Step:**

Run:
```bash
grep -c '^###' docs/configuration.md && python -c "import xml.etree.ElementTree as E; E.parse('docs/logo.svg'); print('svg ok')"
```

Expected:
A count >= 13, then "svg ok".

**Notes:**

Copy README.md lines from "## Configuring Simulation" to just before "## Acknowledgement". Fix relative links if any (e.g. `.work/adr/...` becomes `../.work/adr/...`).
