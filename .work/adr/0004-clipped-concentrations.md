# Rate laws read clipped concentrations; force_positive removed

_ADR 0004, recorded 2026-11-10 (session 2026-11-10__derive-math-fixes)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

A solver overshoot below 0 could drive a self-reacting species further negative (dA/dt = -2kA² < 0) or run a reaction backwards, and `force_positive` only patched the net derivative, which broke mass balance (dy/dt was no longer S·r). Every rate law (mass action and Michaelis-Menten) now reads the clipped concentration max(c, 0), always on, so dy/dt = S·r holds exactly and a species at 0 is no longer consumed; `force_positive` is removed outright (breaking: old callers get a TypeError). In the analytic Jacobian a species below 0 contributes 0, and at exactly 0 the derivative of the c > 0 side is used, so states starting at 0 keep an informative Jacobian.
