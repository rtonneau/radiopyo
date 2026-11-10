# G(H3O+) is taken equal to G(e_aq)

_ADR 0002, recorded 2026-11-09 (session 2026-11-09__review-this-repo-to-improve-it-keep-all-actual-features-jocobian)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

When LET-dependent Boscolo 2020 G-values are applied, G(H3O_plus) is set to G(e_aq) rather than read from the CSV's own H3O_plus column. This is intended (charge balance of the primary yields) and is kept; the unused CSV column is documented as such rather than wired in.
