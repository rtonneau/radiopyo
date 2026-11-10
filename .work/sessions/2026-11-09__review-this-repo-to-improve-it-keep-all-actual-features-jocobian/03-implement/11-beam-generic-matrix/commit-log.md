# Ticket 11: beam-generic-matrix

**Status:** ✅ Done

## Local Test Result

BLOCKED before implementation: no repository file changed, so the Verification Step was not run against a change. Baseline check on the unchanged tree: `python -m pytest -q tests/golden/test_golden.py::test_beam_matrices_match_golden tests/physics` gave 95 passed.

Feasibility check: `.scratch/tests/2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian/beam11_compare.py`, run with the base package (`git archive 8c89649`, base copy deleted afterwards) on `PYTHONPATH`. It prototypes one generic `make_matrix`: per-timer transition generators plus a single walk up to `min(stop, time_at_max_dose)`. Each beam's base `make_matrix`/`make_time_intervals` is compared bit for bit with the prototype on 257 (beam, window) cases. Result: 197 identical, 60 different. Full report in `evidence/11-generic-vs-base.txt`; macro extras (max_dose=inf, micro_start_time) in `evidence/11-macro-extra.txt`.
- PulsedBeam: 156/156 identical (matrix and intervals). The cases cover 8 configs: max_dose finite and inf, start_time 0 and 1e-3, n_pulse, frequency, high duty, single pulse. Windows include start/stop inside a pulse, exactly at pulse start and end, at T, and before start_time.
- ConstantBeam / SinglePulseBeam: 37/66 identical.
- MacroPulsedBeam: 4/35 identical.

## Review Notes

- Parent: base 8c8964990989bbd72ff06e5af4165b825bf2d0fc; pre-existing changes: gps session state only.

Implementer check against the Acceptance Criteria (no review rounds): AC1 (one walk up to `min(stop, time_at_max_dose)`, the four versions removed) and AC2 (beam goldens exact, quirks included) contradict each other for MacroPulsedBeam. They also clash for Constant/SinglePulse unless a type-specific branch is added. Per the ticket Notes nothing was implemented (see Blockers). AC3 (chunked-run golden) uses a PulsedBeam only. The prototype walk is bit-identical for PulsedBeam, so a generic version would leave that golden unchanged.

- Parent decision: ACs 1 and 2 contradict (generic walk bit-identical only for PulsedBeam 156/156; Macro 4/35, Constant/SinglePulse 37/66, evidence/11-generic-vs-base.txt). As anticipated by the ticket Notes, the generic make_matrix is merged into ticket 15, whose golden regeneration needs user approval anyway. Ticket 11 will be closed with ticket 15's commit.

- Parent acceptance: AC1 (one generic make_matrix walking timer transitions up to min(stop, time_at_max_dose)) delivered in ticket 15, commit 10978f3. AC2 (bit-identical quirks) dropped on purpose; the beam goldens were re-recorded there after user approval. AC3 (chunked-run golden) still passes.

## Blockers / Challenges

BLOCKED (ticket Notes case): reproducing the golden-pinned quirks needs beam-specific code paths. That is the hand-written versions under another name.

1. MacroPulsedBeam (decisive). Base `make_matrix` counts trains (`tot_trains = max_dose / dose_per_train`). Neither `stop` nor T bounds it.
   - Golden `macropulsed__3` (window 0..0.05): 53 rows running to 0.2105, three trains past stop. The walk gives 22 rows ending at 0.05.
   - Matrix end: the base partial train lasts `on_time * fraction` and ends at 0.2105, while T = 0.210402. The walk ends at T.
   - Duplicate `[0.01, 1e6]` row at train start (golden `__0`, `__3`).
   - Golden `__1` (start 0.01) starts `[0.01, 0], [0.01, 1e6]`.
   - Golden `__2` (start 0.05, in the OFF gap) starts `[0.05, 1e6]`, so the beam is wrongly ON. Its train times are `start + time_till_next_train` (0.11000000000000001 instead of 0.11).
   - With micro_start_time=1e-5 the micro start is added twice (rows 0.01001, 0.01002).
   - The first train uses the absolute micro timer; later trains are shifted copies of a relative `PulsedBeam.make_matrix`.
2. ConstantBeam / SinglePulseBeam. When stop < T and the beam is ON at stop, base closes with `[stop, dose_rate_peak]` (golden `constant__2`, `constant__5`, `singlepulse__4`). PulsedBeam closes the same situation with `[stop, 0]`. One walk can match both only with an `is_continuous` branch in the closing row. Choosing either convention changes the other beam. Boundary duplicates not in the golden are also lost: stop == T gives `[T,0],[T,0]`, and stop == start_time gives `[s0,dr],[s0,dr]`.

What a generic version gives instead: Pulsed unchanged. Constant/SinglePulse close with 0. Macro is bounded by min(stop, T), has no duplicate rows and is correct at an OFF start. It also terminates with max_dose=inf (62 rows for 0..0.25), which covers part of ticket 15 AC2. Suggestion: merge ticket 11 into 15 and regenerate the beam goldens there with user approval. Prototype design: transition generators on the timers (Timer: `t_on += period`, `t_off = t_on + on_time`; continuous: ON forever; MacroTimer: trains plus subtimer pulses cut at train end), plus one `ParticleBeam.make_matrix` walk. `make_time_intervals` is already a single shared implementation and was left as is.

**Blocked (2026-11-09T23:40:46.374Z):** merged into ticket 15: generic make_matrix cannot reproduce Macro/Constant quirks bit-identically; done with ticket 15 golden regeneration

## Commits

- 10978f3 fix: generic beam make_matrix, macro beam and timer fixes, average dose rate from start_time (tickets 15, 11)

## Time Spent

8h 30m (ticket-start.js to ticket-complete.js)
