# Session Summary: Relax beam ON time constrain from 1e-6 to 0.4e-6s

**Session ID:** 2026-11-13__relax-beam-on-time-constrain-from-1e-6-to-0.4e-6s
**Created:** 2026-11-13T14:41:44.139Z
**Finished:** 2026-11-13T15:04:28.102Z
**Status:** Complete

## Grill

- Resume: [01-grill/resume.md](01-grill/resume.md)

## Plan

- Plan: [02-plan/plan.md](02-plan/plan.md)

## Tickets

- ✅ 01 radkin-min-on-time — [spec](02-plan/tickets/01-radkin-min-on-time.md) · [log](03-implement/01-radkin-min-on-time/commit-log.md)
- ✅ 02 radiopyo-min-on-time — [spec](02-plan/tickets/02-radiopyo-min-on-time.md) · [log](03-implement/02-radiopyo-min-on-time/commit-log.md)

## Changelog

- **Bump:** minor
- **Fragment:** [2026-11-13__relax-beam-on-time-constrain-from-1e-6-to-0.4e-6s.md](../../changelog/2026-11-13__relax-beam-on-time-constrain-from-1e-6-to-0.4e-6s.md) (merged into the CHANGELOG at release)

## Branch & PR

- **Branch:** `feat/min-beam-on-time-0.4us`
- **Base:** `main`
- **Pull request:** #16 (new): feat: Lower MIN_BEAM_ON_TIME to 0.4 µs in radiopyo and radkin (feat/min-beam-on-time-0.4us → main) — https://github.com/rtonneau/radiopyo/pull/16

## Timeline

| When | Phase | Event | Details | Files |
|---|---|---|---|---|
| 2026-11-13 15:41 | grill | session_started |  | [01-grill/resume.md](01-grill/resume.md) |
| 2026-11-13 15:43 | grill | auto_started | target: finish, steps: write:grill → plan → write:plan → ship → finish, shipMode: delegate |  |
| 2026-11-13 15:43 | plan-not-started | grill_written |  | [01-grill/resume.md](01-grill/resume.md) |
| 2026-11-13 15:43 | plan | plan_started |  | [02-plan/plan.md](02-plan/plan.md) |
| 2026-11-13 15:49 | plan | branch_created | branch: feat/min-beam-on-time-0.4us, base: main |  |
| 2026-11-13 15:49 | ship | plan_written | tickets: 2 | [02-plan/plan.md](02-plan/plan.md), [02-plan/tickets/01-radkin-min-on-time.md](02-plan/tickets/01-radkin-min-on-time.md), [02-plan/tickets/02-radiopyo-min-on-time.md](02-plan/tickets/02-radiopyo-min-on-time.md) |
| 2026-11-13 15:49 | ship | ticket_started | ticket: 01-radkin-min-on-time | [02-plan/tickets/01-radkin-min-on-time.md](02-plan/tickets/01-radkin-min-on-time.md), [03-implement/01-radkin-min-on-time/commit-log.md](03-implement/01-radkin-min-on-time/commit-log.md) |
| 2026-11-13 15:58 | ship | ticket_done | ticket: 01-radkin-min-on-time, commit: 5b3d9b9 | [03-implement/01-radkin-min-on-time/commit-log.md](03-implement/01-radkin-min-on-time/commit-log.md) |
| 2026-11-13 15:58 | ship | ticket_started | ticket: 02-radiopyo-min-on-time | [02-plan/tickets/02-radiopyo-min-on-time.md](02-plan/tickets/02-radiopyo-min-on-time.md), [03-implement/02-radiopyo-min-on-time/commit-log.md](03-implement/02-radiopyo-min-on-time/commit-log.md) |
| 2026-11-13 16:04 | finish-pending | ticket_done | ticket: 02-radiopyo-min-on-time, commit: cdc4f7e | [03-implement/02-radiopyo-min-on-time/commit-log.md](03-implement/02-radiopyo-min-on-time/commit-log.md) |
| 2026-11-13 16:04 | finish-pending | changelog_written | bump: minor |  |
| 2026-11-13 16:04 | finish-pending | pr_opened | url: https://github.com/rtonneau/radiopyo/pull/16 |  |
| 2026-11-13 16:04 | finished | session_finished |  | [INDEX.md](INDEX.md) |

## Next

Start a new feature with /gps start <next-feature>
