# n·k convention applies to scalar k only

_ADR 0001, recorded 2026-11-09 (session 2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

For `n A -> B`, `KReaction` divides the configured k by n only when k is a single number; polynomial k(T), Arrhenius and equilibrium k are used as given. We keep this behaviour (and document it) rather than unifying it, because unifying would silently change results of user configs with polynomial k on nA reactions, while bundled configs only use scalar k.
