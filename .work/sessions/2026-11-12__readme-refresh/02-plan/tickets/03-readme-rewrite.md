# Ticket 03: readme-rewrite

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] README.md under 150 lines, sections: logo (docs/logo.svg) + badges (Tests workflow, version 0.6.0, Python >=3.11, LGPL-3.0, Ruff, mypy, DOI 10.1088/1361-6560/ae1fcc), short pitch, Installation (pip+git, editable), Update (`pip install --upgrade git+https://github.com/rtonneau/radiopyo.git`, `git pull` for editable), Usage (single run + chunked pulsed run, max_dose by keyword), Configuration (condensed table of TOML sections and beam types linking to docs/configuration.md), Development (short), Citation (Romain Tonneau et al 2025 Phys. Med. Biol. 70 235021, DOI 10.1088/1361-6560/ae1fcc, link https://iopscience.iop.org/article/10.1088/1361-6560/ae1fcc/meta, BibTeX), Acknowledgement (kept).
- [ ] Every relative link resolves.

**Files to Touch:**
- `README.md`

**Verification Step:**

Run:
```bash
wc -l README.md && grep -oE '\]\((docs|\.github)[^)]*\)' README.md
```

Expected:
Fewer than 150 lines; each listed path exists.

**Notes:**

Plain, direct prose. BibTeX title: check the paper page if possible; otherwise use the DOI-only fields (author, journal, volume, number 23?, pages 235021, year 2025, doi) and do not invent a title.
