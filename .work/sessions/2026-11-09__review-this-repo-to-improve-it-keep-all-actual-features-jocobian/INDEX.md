# Session Summary: review this repo to improve it. Keep all actual features. Jocobian will be one of the first next job. Just improve flow and readability. Don't break how code is used by user. Other things can be adapted if needed

**Session ID:** 2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian
**Created:** 2026-11-09T20:05:46.317Z
**Finished:** 2026-11-10T09:31:19.755Z
**Status:** Complete

## Grill

- Resume: [01-grill/resume.md](01-grill/resume.md)

## Plan

- Plan: [02-plan/plan.md](02-plan/plan.md)

## Tickets

- ✅ 01 golden-outputs — [spec](02-plan/tickets/01-golden-outputs.md) · [log](03-implement/01-golden-outputs/commit-log.md)
- ✅ 02 port-tests — [spec](02-plan/tickets/02-port-tests.md) · [log](03-implement/02-port-tests/commit-log.md)
- ✅ 03 packaging-hygiene — [spec](02-plan/tickets/03-packaging-hygiene.md) · [log](03-implement/03-packaging-hygiene/commit-log.md)
- ✅ 04 dead-code-lint-zero — [spec](02-plan/tickets/04-dead-code-lint-zero.md) · [log](03-implement/04-dead-code-lint-zero/commit-log.md)
- ✅ 05 solver-core — [spec](02-plan/tickets/05-solver-core.md) · [log](03-implement/05-solver-core/commit-log.md)
- ✅ 06 run-flow — [spec](02-plan/tickets/06-run-flow.md) · [log](03-implement/06-run-flow/commit-log.md)
- ✅ 07 simenv-builders — [spec](02-plan/tickets/07-simenv-builders.md) · [log](03-implement/07-simenv-builders/commit-log.md)
- ✅ 08 split-sim-utils — [spec](02-plan/tickets/08-split-sim-utils.md) · [log](03-implement/08-split-sim-utils/commit-log.md)
- ✅ 09 chemistry-cleanup — [spec](02-plan/tickets/09-chemistry-cleanup.md) · [log](03-implement/09-chemistry-cleanup/commit-log.md)
- ✅ 10 beam-split — [spec](02-plan/tickets/10-beam-split.md) · [log](03-implement/10-beam-split/commit-log.md)
- ✅ 11 beam-generic-matrix — [spec](02-plan/tickets/11-beam-generic-matrix.md) · [log](03-implement/11-beam-generic-matrix/commit-log.md)
- ✅ 12 fix-run-burst — [spec](02-plan/tickets/12-fix-run-burst.md) · [log](03-implement/12-fix-run-burst/commit-log.md)
- ✅ 13 fix-temperature — [spec](02-plan/tickets/13-fix-temperature.md) · [log](03-implement/13-fix-temperature/commit-log.md)
- ✅ 14 fix-include-merge — [spec](02-plan/tickets/14-fix-include-merge.md) · [log](03-implement/14-fix-include-merge/commit-log.md)
- ✅ 15 fix-macro-beam — [spec](02-plan/tickets/15-fix-macro-beam.md) · [log](03-implement/15-fix-macro-beam/commit-log.md)
- ✅ 16 fix-species-order — [spec](02-plan/tickets/16-fix-species-order.md) · [log](03-implement/16-fix-species-order/commit-log.md)
- ✅ 17 fix-arrhenius-enzymes-add — [spec](02-plan/tickets/17-fix-arrhenius-enzymes-add.md) · [log](03-implement/17-fix-arrhenius-enzymes-add/commit-log.md)
- ✅ 18 docs — [spec](02-plan/tickets/18-docs.md) · [log](03-implement/18-docs/commit-log.md)

## Changelog

- **Bump:** minor
- **Reason:** breaking removals (RON format, run_multi_doses, e_LET, use_jac) in a 0.x version
- **Fragment:** [2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian.md](../../changelog/2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian.md) (merged into the CHANGELOG at release)

## Branch & PR

- **Branch:** `refactor/readability-and-jacobian-prep`
- **Base:** `main`
- **Pull request:** #1 (new): refactor: Refactor for readability and Jacobian prep; fix run burst, temperature, includes, beams and species order (refactor/readability-and-jacobian-prep → main) — https://github.com/rtonneau/radiopyo/pull/1

## Timeline

| When | Phase | Event | Details | Files |
|---|---|---|---|---|
| 2026-11-09 21:05 | grill | session_started |  | [01-grill/resume.md](01-grill/resume.md) |
| 2026-11-09 21:25 | plan-not-started | grill_written |  | [01-grill/resume.md](01-grill/resume.md) |
| 2026-11-09 21:25 | plan | plan_started |  | [02-plan/plan.md](02-plan/plan.md) |
| 2026-11-09 21:54 | plan | auto_started | target: finish, steps: write:plan → ship → finish, shipMode: delegate |  |
| 2026-11-09 21:54 | plan | branch_created | branch: refactor/readability-and-jacobian-prep, base: main |  |
| 2026-11-09 21:54 | ship | plan_written | tickets: 18 | [02-plan/plan.md](02-plan/plan.md), [02-plan/tickets/01-golden-outputs.md](02-plan/tickets/01-golden-outputs.md), [02-plan/tickets/02-port-tests.md](02-plan/tickets/02-port-tests.md), [02-plan/tickets/03-packaging-hygiene.md](02-plan/tickets/03-packaging-hygiene.md), [02-plan/tickets/04-dead-code-lint-zero.md](02-plan/tickets/04-dead-code-lint-zero.md), [02-plan/tickets/05-solver-core.md](02-plan/tickets/05-solver-core.md), [02-plan/tickets/06-run-flow.md](02-plan/tickets/06-run-flow.md), [02-plan/tickets/07-simenv-builders.md](02-plan/tickets/07-simenv-builders.md), [02-plan/tickets/08-split-sim-utils.md](02-plan/tickets/08-split-sim-utils.md), [02-plan/tickets/09-chemistry-cleanup.md](02-plan/tickets/09-chemistry-cleanup.md), [02-plan/tickets/10-beam-split.md](02-plan/tickets/10-beam-split.md), [02-plan/tickets/11-beam-generic-matrix.md](02-plan/tickets/11-beam-generic-matrix.md), [02-plan/tickets/12-fix-run-burst.md](02-plan/tickets/12-fix-run-burst.md), [02-plan/tickets/13-fix-temperature.md](02-plan/tickets/13-fix-temperature.md), [02-plan/tickets/14-fix-include-merge.md](02-plan/tickets/14-fix-include-merge.md), [02-plan/tickets/15-fix-macro-beam.md](02-plan/tickets/15-fix-macro-beam.md), [02-plan/tickets/16-fix-species-order.md](02-plan/tickets/16-fix-species-order.md), [02-plan/tickets/17-fix-arrhenius-enzymes-add.md](02-plan/tickets/17-fix-arrhenius-enzymes-add.md), [02-plan/tickets/18-docs.md](02-plan/tickets/18-docs.md) |
| 2026-11-09 21:54 | ship | ticket_started | ticket: 01-golden-outputs | [02-plan/tickets/01-golden-outputs.md](02-plan/tickets/01-golden-outputs.md), [03-implement/01-golden-outputs/commit-log.md](03-implement/01-golden-outputs/commit-log.md) |
| 2026-11-09 22:31 | ship | ticket_done | ticket: 01-golden-outputs, commit: dba163a | [03-implement/01-golden-outputs/commit-log.md](03-implement/01-golden-outputs/commit-log.md) |
| 2026-11-09 22:31 | ship | ticket_started | ticket: 02-port-tests | [02-plan/tickets/02-port-tests.md](02-plan/tickets/02-port-tests.md), [03-implement/02-port-tests/commit-log.md](03-implement/02-port-tests/commit-log.md) |
| 2026-11-09 22:38 | ship | ticket_done | ticket: 02-port-tests, commit: 9424c85 | [03-implement/02-port-tests/commit-log.md](03-implement/02-port-tests/commit-log.md) |
| 2026-11-09 22:39 | ship | ticket_started | ticket: 03-packaging-hygiene | [02-plan/tickets/03-packaging-hygiene.md](02-plan/tickets/03-packaging-hygiene.md), [03-implement/03-packaging-hygiene/commit-log.md](03-implement/03-packaging-hygiene/commit-log.md) |
| 2026-11-09 22:41 | ship | ticket_done | ticket: 03-packaging-hygiene, commit: 9a8e758 | [03-implement/03-packaging-hygiene/commit-log.md](03-implement/03-packaging-hygiene/commit-log.md) |
| 2026-11-09 22:41 | ship | ticket_started | ticket: 04-dead-code-lint-zero | [02-plan/tickets/04-dead-code-lint-zero.md](02-plan/tickets/04-dead-code-lint-zero.md), [03-implement/04-dead-code-lint-zero/commit-log.md](03-implement/04-dead-code-lint-zero/commit-log.md) |
| 2026-11-09 22:49 | ship | ticket_done | ticket: 04-dead-code-lint-zero, commit: ad88152 | [03-implement/04-dead-code-lint-zero/commit-log.md](03-implement/04-dead-code-lint-zero/commit-log.md) |
| 2026-11-09 22:49 | ship | ticket_started | ticket: 05-solver-core | [02-plan/tickets/05-solver-core.md](02-plan/tickets/05-solver-core.md), [03-implement/05-solver-core/commit-log.md](03-implement/05-solver-core/commit-log.md) |
| 2026-11-09 23:24 | ship | ticket_done | ticket: 05-solver-core, commit: f567a6f | [03-implement/05-solver-core/commit-log.md](03-implement/05-solver-core/commit-log.md) |
| 2026-11-09 23:24 | ship | ticket_started | ticket: 06-run-flow | [02-plan/tickets/06-run-flow.md](02-plan/tickets/06-run-flow.md), [03-implement/06-run-flow/commit-log.md](03-implement/06-run-flow/commit-log.md) |
| 2026-11-09 23:45 | ship | ticket_done | ticket: 06-run-flow, commit: d24bc92 | [03-implement/06-run-flow/commit-log.md](03-implement/06-run-flow/commit-log.md) |
| 2026-11-09 23:45 | ship | ticket_started | ticket: 07-simenv-builders | [02-plan/tickets/07-simenv-builders.md](02-plan/tickets/07-simenv-builders.md), [03-implement/07-simenv-builders/commit-log.md](03-implement/07-simenv-builders/commit-log.md) |
| 2026-11-09 23:56 | ship | ticket_done | ticket: 07-simenv-builders, commit: 3f8035f | [03-implement/07-simenv-builders/commit-log.md](03-implement/07-simenv-builders/commit-log.md) |
| 2026-11-09 23:56 | ship | ticket_started | ticket: 08-split-sim-utils | [02-plan/tickets/08-split-sim-utils.md](02-plan/tickets/08-split-sim-utils.md), [03-implement/08-split-sim-utils/commit-log.md](03-implement/08-split-sim-utils/commit-log.md) |
| 2026-11-10 00:02 | ship | ticket_done | ticket: 08-split-sim-utils, commit: aefda6e | [03-implement/08-split-sim-utils/commit-log.md](03-implement/08-split-sim-utils/commit-log.md) |
| 2026-11-10 00:02 | ship | ticket_started | ticket: 09-chemistry-cleanup | [02-plan/tickets/09-chemistry-cleanup.md](02-plan/tickets/09-chemistry-cleanup.md), [03-implement/09-chemistry-cleanup/commit-log.md](03-implement/09-chemistry-cleanup/commit-log.md) |
| 2026-11-10 00:22 | ship | ticket_done | ticket: 09-chemistry-cleanup, commit: a0bee93 | [03-implement/09-chemistry-cleanup/commit-log.md](03-implement/09-chemistry-cleanup/commit-log.md) |
| 2026-11-10 00:22 | ship | ticket_started | ticket: 10-beam-split | [02-plan/tickets/10-beam-split.md](02-plan/tickets/10-beam-split.md), [03-implement/10-beam-split/commit-log.md](03-implement/10-beam-split/commit-log.md) |
| 2026-11-10 00:31 | ship | ticket_done | ticket: 10-beam-split, commit: cd77b2b | [03-implement/10-beam-split/commit-log.md](03-implement/10-beam-split/commit-log.md) |
| 2026-11-10 00:31 | ship | ticket_started | ticket: 11-beam-generic-matrix | [02-plan/tickets/11-beam-generic-matrix.md](02-plan/tickets/11-beam-generic-matrix.md), [03-implement/11-beam-generic-matrix/commit-log.md](03-implement/11-beam-generic-matrix/commit-log.md) |
| 2026-11-10 00:40 | ship | ticket_blocked | ticket: 11-beam-generic-matrix, reason: merged into ticket 15: generic make_matrix cannot reproduce Macro/Constant quirks bit-identically; done with ticket 15 golden regeneration | [03-implement/11-beam-generic-matrix/commit-log.md](03-implement/11-beam-generic-matrix/commit-log.md) |
| 2026-11-10 00:40 | ship | ticket_started | ticket: 12-fix-run-burst | [02-plan/tickets/12-fix-run-burst.md](02-plan/tickets/12-fix-run-burst.md), [03-implement/12-fix-run-burst/commit-log.md](03-implement/12-fix-run-burst/commit-log.md) |
| 2026-11-10 07:05 | ship | ticket_done | ticket: 12-fix-run-burst, commit: af5b8ca | [03-implement/12-fix-run-burst/commit-log.md](03-implement/12-fix-run-burst/commit-log.md) |
| 2026-11-10 07:06 | ship | ticket_started | ticket: 13-fix-temperature | [02-plan/tickets/13-fix-temperature.md](02-plan/tickets/13-fix-temperature.md), [03-implement/13-fix-temperature/commit-log.md](03-implement/13-fix-temperature/commit-log.md) |
| 2026-11-10 07:16 | ship | ticket_done | ticket: 13-fix-temperature, commit: bafa16f | [03-implement/13-fix-temperature/commit-log.md](03-implement/13-fix-temperature/commit-log.md) |
| 2026-11-10 07:16 | ship | ticket_started | ticket: 14-fix-include-merge | [02-plan/tickets/14-fix-include-merge.md](02-plan/tickets/14-fix-include-merge.md), [03-implement/14-fix-include-merge/commit-log.md](03-implement/14-fix-include-merge/commit-log.md) |
| 2026-11-10 08:10 | ship | ticket_done | ticket: 14-fix-include-merge, commit: 662277a | [03-implement/14-fix-include-merge/commit-log.md](03-implement/14-fix-include-merge/commit-log.md) |
| 2026-11-10 08:10 | ship | ticket_started | ticket: 15-fix-macro-beam | [02-plan/tickets/15-fix-macro-beam.md](02-plan/tickets/15-fix-macro-beam.md), [03-implement/15-fix-macro-beam/commit-log.md](03-implement/15-fix-macro-beam/commit-log.md) |
| 2026-11-10 09:01 | ship | ticket_done | ticket: 15-fix-macro-beam, commit: 10978f3 | [03-implement/15-fix-macro-beam/commit-log.md](03-implement/15-fix-macro-beam/commit-log.md) |
| 2026-11-10 09:01 | ship | ticket_done | ticket: 11-beam-generic-matrix, commit: 10978f3 | [03-implement/11-beam-generic-matrix/commit-log.md](03-implement/11-beam-generic-matrix/commit-log.md) |
| 2026-11-10 09:01 | ship | ticket_started | ticket: 16-fix-species-order | [02-plan/tickets/16-fix-species-order.md](02-plan/tickets/16-fix-species-order.md), [03-implement/16-fix-species-order/commit-log.md](03-implement/16-fix-species-order/commit-log.md) |
| 2026-11-10 09:10 | ship | ticket_done | ticket: 16-fix-species-order, commit: 5b3b9e2 | [03-implement/16-fix-species-order/commit-log.md](03-implement/16-fix-species-order/commit-log.md) |
| 2026-11-10 09:10 | ship | ticket_started | ticket: 17-fix-arrhenius-enzymes-add | [02-plan/tickets/17-fix-arrhenius-enzymes-add.md](02-plan/tickets/17-fix-arrhenius-enzymes-add.md), [03-implement/17-fix-arrhenius-enzymes-add/commit-log.md](03-implement/17-fix-arrhenius-enzymes-add/commit-log.md) |
| 2026-11-10 10:10 | ship | ticket_done | ticket: 17-fix-arrhenius-enzymes-add, commit: 2ad881e | [03-implement/17-fix-arrhenius-enzymes-add/commit-log.md](03-implement/17-fix-arrhenius-enzymes-add/commit-log.md) |
| 2026-11-10 10:10 | ship | ticket_started | ticket: 18-docs | [02-plan/tickets/18-docs.md](02-plan/tickets/18-docs.md), [03-implement/18-docs/commit-log.md](03-implement/18-docs/commit-log.md) |
| 2026-11-10 10:30 | finish-pending | ticket_done | ticket: 18-docs, commit: 4ea057d | [03-implement/18-docs/commit-log.md](03-implement/18-docs/commit-log.md) |
| 2026-11-10 10:31 | finish-pending | changelog_written | bump: minor |  |
| 2026-11-10 10:31 | finish-pending | pr_opened | url: https://github.com/rtonneau/radiopyo/pull/1 |  |
| 2026-11-10 10:31 | finished | session_finished |  | [INDEX.md](INDEX.md) |

## Next

Start a new feature with /gps start <next-feature>
