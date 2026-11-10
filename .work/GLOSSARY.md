# radiopyo glossary

radiopyo simulates water radiolysis and radiobiology chemistry under UHDR and conventional irradiation by solving a system of chemical-kinetics ODEs, to compare how dose rate changes the chemistry.

<!-- Format adapted from mattpocock-skills' domain-modeling (GLOSSARY-FORMAT.md): https://github.com/mattpocock/skills.git -->

<!--
How to add a term (the moment it is resolved, not in a batch at the end):

**Term**:
One or two sentences: what it IS, not what it does.
_Avoid_: the other words people use for it

- Be opinionated: one word per concept; list the others under _Avoid_.
- Only terms specific to this project. General programming concepts and
  implementation details don't belong here: this is a glossary, not a spec.
- Group terms under "### <area>" subheadings once natural clusters emerge.
-->

## Language

**Dynamic species**:
A species whose concentration is solved by the ODE system; it is part of the state vector `y` (µmol/L).
_Avoid_: sim species, solved species

**Fixed species**:
A species held at a constant concentration (`[concentrations.fixed]`, e.g. H2O); appended to the state inside `derive` but never solved.
_Avoid_: constant species, cst species

**Initial burst**:
The radiolytic yield of the first `T_ZERO` (1e-7 s) of a beam-ON interval, added analytically to `y0` instead of being integrated.
_Avoid_: T0 yield, initial radiolytic injection

**n·k convention**:
For a scalar k on `n A -> B`, the configured value is n·k and is divided by n when the reaction is built (see ADR 0001).
_Avoid_: 2k rule

**Golden output**:
A reference result recorded from the code before a refactor, against which the refactored code must match (to a stated tolerance).
_Avoid_: snapshot, baseline

**Result-changing fix**:
A bug fix that changes numerical simulation output; delivered as its own ticket with before/after numbers for approval.
_Avoid_: behaviour fix

**Negative product**:
A species written `- species` on the product side of a reaction: the reaction consumes it although it is not in the rate law, so the reaction is slowed by a depletion switch as that species runs out.
_Avoid_: ope, negative species

**Depletion switch**:
The smooth factor c⁺/(c⁺ + ε), ε = `DEPLETION_LIMIT`, multiplying the rate of a reaction (k-reaction or enzymatic) for each of its negative products (ADR 0005); it replaces the former hard depletion block.
_Avoid_: depletion block, negative-product blocking, zeroing

**Clipped concentration**:
The concentration a rate law reads, max(c, 0): a species driven below 0 by a solver overshoot counts as absent, so it can neither be consumed further nor run a reaction backwards.
_Avoid_: positive part, clamped state

**Analytic Jacobian**:
The matrix d(dy/dt)/dy of the dynamic species computed in closed form from the reaction matrices, given to the stiff solvers in place of a finite-difference estimate.
_Avoid_: jac, exact Jacobian
