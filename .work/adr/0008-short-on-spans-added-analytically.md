# Short beam ON spans: error when configured, analytic yield at run time

_ADR 0008, recorded 2026-11-11 (session 2026-11-11__beam-on-time-min-and-single-pulse-off-step)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

A beam ON span shorter than T_ZERO could not be integrated (ValueError) and one shorter than about 1e-6 s is resolved poorly, so a beam whose *configured* ON time (pulse `on_time`, `micro_on_time`, macro train `on_time`, single pulse or `max_dose`-limited constant beam) is below `MIN_BEAM_ON_TIME` = 1e-6 s raises `ShortBeamOnError` when it is built. A span that only becomes short at run time (last pulse cut by `max_dose`, pulse clipped by `t_span`) is a legitimate input, so it is not integrated: its radiolytic yield over the whole span is added analytically at the span dose rate, a warning is logged and the span is recorded in `ODEResult.short_on_spans` (start, stop, dose rate); the stored step carries a dose rate of 0 so the production integral counts the yield once.

**Consequences:** the chemistry inside a short ON span is not resolved (the radicals are created all at once at its end); a result loaded from HDF5 has no short ON span record.

## Amendment 2026-11-13

`MIN_BEAM_ON_TIME` is now 0.4e-6 s (still 4x `T_ZERO` = 1e-7 s) instead of 1e-6 s, because real UHDR machines deliver 0.4 to 1 us pulses. Everything above holds with the new value; radkin carries the same threshold (rev `088f06b`).
