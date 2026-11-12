# Follow-up from gps session 2026-11-12__chunked-run-scaling

F1 (Worth exploring): with opt-in point thinning (points per chunk, ON/OFF) the production /
consumption integrals of ODEResult (trapezoid over the stored points: integrate_production_of,
integrate_consumption_of, detail_*, species_balance_evolution) become coarser. The session only
documents it and logs a warning. Proposed fix: integrate the reaction rates per chunk at the full
solver resolution during the run (e.g. accumulate the trapezoid of reaction_rates over the solver
steps, or integrate an augmented state) and store the per-reaction totals per chunk, so integrals
stay exact whatever the stored point density. Touches result math: needs golden checks / an ADR.
Files: radiopyo/simulation/ode_result.py, radiopyo/simulation/unit_cell.py, radiopyo/simulation/ode_functions_solver.py
