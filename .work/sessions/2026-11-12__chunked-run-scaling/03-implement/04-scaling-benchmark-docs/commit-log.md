# Ticket 04: scaling-benchmark-docs

**Status:** ✅ Done

## Local Test Result

`python examples/bench_chunked_scaling.py --pulses 100 400 1600` (config_Labarbe, peak 1e4 Gy/s, 10 us ON, 10 ms period):

| pulses | chunks | points | wall [s] | ms/chunk | first 10% | last 10% |
|---|---|---|---|---|---|---|
| 100 | 200 | 95969 | 10.5 | 52.6 | 56.2 | 47.9 |
| 400 | 800 | 365098 | 41.1 | 51.4 | 53.8 | 50.0 |
| 1600 | 3200 | 1425321 | 160.1 | 50.0 | 53.3 | 56.4 |

With `--points-on 20 --points-off 20`:

| pulses | chunks | points | wall [s] | ms/chunk | first 10% | last 10% |
|---|---|---|---|---|---|---|
| 1600 | 3200 | 64000 | 195.9 | 61.2 | 53.9 | 56.3 |
| 10000 | 20000 | 400000 | 986.4 | 49.3 | 54.8 | 38.5 |

Per-chunk cost flat (~50 ms) from 200 to 20000 chunks; grill baseline was 455 s for 1600 pulses (now 160 s). After merging origin/main (4f9991e): `python -m pytest -q` 786 passed; ruff "All checks passed!"; mypy "Success: no issues found in 42 source files".

## Review Notes

Mode inline.
- Benchmark in module banner layout, `--pulses`, `--points-on/--points-off`; prints chunks, stored points, wall, ms/chunk, first/last 10% per-chunk cost.
- Docs: origin/main (merged in 4f9991e, the branch had started from a stale main) moved the run section of the README into docs/configuration.md and released 0.7.0, so the thinning documentation went to docs/configuration.md and the CHANGELOG entries to a new [Unreleased] section above 0.7.0. CLAUDE.md: storage, thinning, benchmark notes.
- Thinning at 1600 pulses is not faster (the solver work dominates); it cuts stored points 22x.

## Blockers / Challenges

The branch started from a stale local main (23 commits behind origin/main); merged origin/main, only the uncommitted README/CHANGELOG edits conflicted and were relocated.

## Commits

- 680dd3a docs: chunked-run scaling benchmark and thinning docs (ticket 04)

## Time Spent

27m (ticket-start.js to ticket-complete.js)
