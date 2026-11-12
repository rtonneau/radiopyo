# Ticket 03: readme-rewrite

**Status:** ✅ Done

## Local Test Result

`wc -l README.md && grep -oE '\]\((docs|\.github)[^)]*\)' README.md`

```
149 README.md
](docs/configuration.md)
```

Also checked: docs/configuration.md, docs/logo.svg and radiopyo/data/ exist.

## Review Notes

All sections present: logo, 7 badges, pitch, Installation, Update, Usage (single and chunked run, `max_dose=` by keyword), Configuration tables linking docs/configuration.md and radiopyo/data/, Development, Citation with BibTeX, Acknowledgement. 149 lines. Trimmed the first draft (162 lines) by condensing prose.

## Blockers / Challenges

None. The paper title and number 23 come from the ticket guidance and were not checked against the paper page.

## Commits

- 0dbaace docs: rewrite README: install, update, usage, citation (ticket 03)

## Time Spent

2m (ticket-start.js to ticket-complete.js)
