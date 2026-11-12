# Paper reproduction uses the code's GF(DR) fit, not the printed one

_ADR 0009, recorded 2026-11-12 (session 2026-11-12__paper-reproduction-notebook)._

<!-- Format adapted from mattpocock-skills' domain-modeling (ADR-FORMAT.md): https://github.com/mattpocock/skills.git -->

Fig. 3B of Tonneau et al. 2025 (PMB 70 235021) prints GF(DR) = 0.982·DR^-0.09, while `GValuesBoscolo._scale_factor` ships 0.972·DR^-0.1097 (commit a0821c3, "second fit with whole Blain data"). The paper reproduction notebook and its golden data run the shipped fit unchanged: the golden file records what radiopyo produces, not a re-tuned model, and the notebook states the gap. Making GF configurable to run the printed fit was rejected for this session.
