# Michaelis-Menten k is used as given (no n·k division)

_ADR 0006, recorded 2026-11-10 (session 2026-11-10__derive-math-fixes)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

For k-reactions a scalar configured k on `n A -> ...` is n·k and divided by n (ADR 0001); for Michaelis-Menten reactions the configured k_value is kept as written, even with a substrate of stoichiometry n, and the stoichiometry -n is applied to that rate. We keep the two conventions different rather than unifying them, so that enzymatic configs mean what they say; the bundled catalase value is corrected instead (1.324e8 -> 6.62e7, taking 6.62e7 as the per-reaction k of Labarbe's Rb10; that reading is not yet checked against Labarbe's MatLab model).
