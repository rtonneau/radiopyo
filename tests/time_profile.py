# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from pathlib import Path

from scipy.integrate import solve_ivp  # type: ignore[import]

import radiopyo as rp
from radiopyo.ode_solver import derive, jac

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.physics.beam import ConstantBeam, PulsedBeam
from radiopyo.sim_env import SimEnv

file = Path(rf"{rp.__path__[0]}\data\reactions.ron")

if __name__ == "__main__":

    import sys

    env = SimEnv.from_ron_file(file)

    # -------------------------------------------------------------------------------- #

    beam = ConstantBeam("e", 5)

    t_span = [1e-9, 500e-6]
    yinit = env.initial_values()
    method = "LSODA"

    profiler = cProfile.Profile()
    profiler.enable()
    sol = solve_ivp(derive,
                    t_span,
                    yinit,
                    method=method,
                    rtol=1e-6,
                    args=(env, beam),
                    jac=jac,
                    max_step=1e-6,
                    )
    profiler.disable()
    stats = pstats.Stats(profiler)
    stats = stats.strip_dirs().sort_stats("ncalls")
    stats.print_stats()

