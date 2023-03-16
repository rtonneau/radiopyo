# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from functools import partial
from pathlib import Path
from typing import (
    TYPE_CHECKING,
)

import numpy as np
from scipy.integrate import solve_ivp  # type: ignore[import]

import radiopyo_vec as rp

# --------------------------------- Internal imports --------------------------------- #
from radiopyo_vec.physics.beam import ConstantBeam, ParticleBeam
from radiopyo_vec.sim_env import SimEnv

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from numpy.typing import NDArray


# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #

def jac(t: float,
        y: NDArray[np.float64],
        env: SimEnv,
        beam: ParticleBeam
        ) -> NDArray[np.float64]:
    pass


def derive(t: float,
           y: NDArray[np.float64],
           env: SimEnv,
           beam: ParticleBeam
           ) -> NDArray[np.float64]:
    pass


if __name__ == "__main__":
    file = Path(rf"{rp.__path__[0]}\data\reactions.ron")
    beam = ConstantBeam("e", 5)
    env = SimEnv.from_ron_file(file)

    t_span = [1e-9, 500e-6]
    yinit = env.initial_values()
    method = "Radau"
    sol = solve_ivp(derive, t_span, yinit, method,
                    rtol=1e-5, args=(env, beam), jac=partial(jac, env=env, beam=beam))
