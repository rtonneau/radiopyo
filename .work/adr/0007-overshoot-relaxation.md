# Overshoot relaxation of species below 0; chunks follow the beam exactly

_ADR 0007, recorded 2026-11-11 (session 2026-11-10__lsoda-stalls-fenton-high-concentration)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

**Status:** amends ADR 0004 (clipped concentrations stay; dy/dt is no longer exactly S·r for a state with a species below 0).

With clipping alone, a stiff radical (e.g. OH_r consumed by RH at 1e9 s⁻¹) sitting near 0 has slope 0 below 0 and −k' above, and LSODA's Newton iterations stall at that kink (steps of ~1e-8 s for minutes; atol only moves the threshold). A dynamic species below 0 now gets an overshoot relaxation: its own row is continued linearly from c = 0⁺ (dc/dt += J_jj(0⁺)·c, the exact Jacobian includes it) while every other species still reads it as 0, so no reaction runs backwards and the mass imbalance is only the refill of an atol-sized hole. A beam discontinuity inside an integration interval also stalled the solver, so chunks are cut exactly at the pulse end (no GROWTH_TIME extension of the ON chunk) and a chunked run holds the dose rate of each chunk constant.

**Considered:** a looser default atol (moves the stall, does not remove it), a smooth softplus clip (drives a species with a tiny source deep below 0), unclipped first-order rates (runs reactions backwards, what ADR 0004 removed).
