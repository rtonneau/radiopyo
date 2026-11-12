# Ticket 01: docs-and-logo

**Status:** ✅ Done

## Local Test Result

Command: `grep -c '^###' docs/configuration.md && python -c "import xml.etree.ElementTree as E; E.parse('docs/logo.svg'); print('svg ok')"`

Output:
```
11
svg ok
```

Result: SVG check passes. The `grep -c` count is 11, below the expected >= 13. The README range "Configuring Simulation" to "Acknowledgement" contains exactly 11 `###` lines (includes, parameters, concentrations.fixed, concentrations.initial, reactions.radiolytic, acid_base, equilibrium, k_reaction, arrhenius, enzymatic, Possible beam definitions) and "Ways of Running Simulation" has no `###` subsections. A verbatim copy cannot reach 13. Needs a parent decision: lower the threshold to 11, or accept added subsection headings (which would break the verbatim criterion).

## Review Notes

- AC1: docs/configuration.md holds the README lines 89-481 verbatim (Configuring Simulation through Ways of Running Simulation), with a title and intro line. The only edit is the ADR reference, now a relative link `../.work/adr/0007-overshoot-relaxation.md`. The other `.work/` mention is in prose and needed no change. Verbatim check is by construction (sed extraction), not diffed afterwards.
- AC2: docs/logo.svg parses as XML. It has an H2O molecule with a radiation wave and a radical dot, next to the "radiopyo" wordmark. Text uses generic font families only, no external resources. Mid-tone colours (#1f6fb2 blue, #5a6270 grey, #e0a100 amber) are used for readability on light and dark backgrounds. Not rendered visually; readability is unverified by eye.
- AC3: README.md is unmodified (`git diff --stat -- README.md` is empty).
- Remaining gap: AC1 is consistent with the verbatim copy, but the Verification Step threshold (>= 13 `###`) cannot be met, see Local Test Result.

## Blockers / Challenges

Verification threshold conflict: the Verification Step expects at least 13 `###` headings, but the verbatim range contains 11. The ticket was not changed and no headings were invented, to keep the copy verbatim. Parent decision needed.

## Commits

- 57997c9 docs: move the config reference to docs/configuration.md and add a logo (ticket 01)

## Time Spent

1m (ticket-start.js to ticket-complete.js)
