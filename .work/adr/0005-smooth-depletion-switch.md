# Smooth depletion switch replaces the depletion block

_ADR 0005, recorded 2026-11-10 (session 2026-11-10__derive-math-fixes)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

The depletion block (a hard step at DEPLETION_LIMIT) made `derive` discontinuous, which the stiff solvers' theory does not allow, and could still let a negative product overshoot within a step. The rate of a k or enzymatic reaction is now multiplied, for each of its negative products, by the smooth switch c⁺/(c⁺ + ε) with ε = DEPLETION_LIMIT = 1e-12 mol/L (half rate at c = ε, under 1 % slowdown above 1e-10 mol/L). This supersedes the piecewise part of ADR 0003: the analytic Jacobian now includes the switch derivative ε/(c + ε)² and is no longer piecewise except at the clipping kink c = 0.
