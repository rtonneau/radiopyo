# Analytic Jacobian on by default, piecewise like derive

_ADR 0003, recorded 2026-11-10 (session 2026-11-10__analytic-jacobian)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

`use_jac` defaults to True: LSODA/BDF/Radau get the analytic Jacobian of `derive` (explicit methods get none), so default results move within solver tolerance, and golden tests check them against the recorded data with a stated tolerance instead of re-recording it. The Jacobian is the derivative on the current side of each on/off switch (a reaction blocked by a depleted negative product contributes 0, rows zeroed by `force_positive` are 0) rather than the smooth Jacobian of the unblocked rates, so that it agrees with `derive` everywhere except exactly at the threshold. The default only becomes True if the benchmark shows the analytic Jacobian is at least as fast and robust as finite differences on the bundled configs; otherwise it stays False.
