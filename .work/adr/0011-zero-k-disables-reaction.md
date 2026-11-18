# A zero rate constant disables the reaction

_ADR 0011, recorded 2026-11-18 (session 2026-11-18__config-file-warnings)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

An `[includes]` override cannot delete a reaction, so users wrote `k_value = 0` to switch one off; that silently built a reaction with log10(k) = -inf. A zero scalar rate constant (`k_value`, `A`, or either side of an acid/base or equilibrium pair, which then drops the whole pair) now removes the reaction from the merged config with a warning, rather than adding an explicit `enabled = false` field: no new syntax, and the existing configs keep their meaning, since a zero-k reaction already contributed nothing. Dropping happens at the config level so the Rust backend, which rebuilds from the config (ADR 0010), sees the same reactions; species only used by a disabled reaction disappear with it.
