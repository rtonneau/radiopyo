# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from functools import partial
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    List,
    Tuple,
    TypedDict,
)

import numpy as np
from scipy.integrate import solve_ivp  # type: ignore[import]

import radiopyo as rp

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.base_chemistry import ReactionRateEnum
from radiopyo.chemistry.k_reactions import KReaction
from radiopyo.physics.beam import ConstantBeam, ParticleBeam
from radiopyo.sim_env import SimEnv

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
    dfdy = np.zeros((env.number_tracked_species(),
                    env.number_tracked_species()))

    # Create a Dict[str, float] with the current cc
    sp_data = env.update_cc_species(y)  # Acid/Base eq. included
    jacs = env.compute_jacobians(sp_data)
    # Compute the Jacobian
    for sp in env.iter_TrackedSpecies():
        for rate in filter(lambda x: isinstance(x, KReaction), sp.kreactions):
            match rate.role:
                case ReactionRateEnum.Production:
                    dfdy[sp.index, :] += jacs[rate.reaction_index]
                case ReactionRateEnum.Consumption:
                    dfdy[sp.index, :] -= jacs[rate.reaction_index]
    return dfdy


def derive(t: float,
           y: NDArray[np.float64],
           env: SimEnv,
           beam: ParticleBeam
           ) -> NDArray[np.float64]:

    # Get the dose_rate for the time t:
    tm = beam.at(t)

    # Create a Dict[str, float] with the current cc
    sp_data = env.update_cc_species(y)  # Acid/Base eq. included

    # Compute production rate values from reaction list
    return env.compute_reactions(tm.dose_rate(), sp_data) * 1e6


class ArgsSolver(TypedDict, total=False):
    max_step: float
    min_step: float
    atol: float
    t_eval: NDArray


def simulate(t_span: List[float],
             y0: NDArray[np.float64],
             env: SimEnv,
             beam: ParticleBeam,
             **kwargs: Any,
             ) -> Any:

    return solve_ivp(derive,
                     t_span,
                     y0,
                     args=(env, beam),
                     jac=jac,
                     method="LSODA",
                     **kwargs,
                     )


def mono_pulse(env: SimEnv,
               beam: ParticleBeam,
               t_eval: NDArray,
               ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:

    t_span = [1e-9, t_eval[-1]]

    kwargs: ArgsSolver = {}
    kwargs["max_step"] = np.inf
    kwargs["min_step"] = 0.0
    kwargs["atol"] = 1e-6
    kwargs["t_eval"] = t_eval

    y0 = env.initial_values()
    sol = simulate(t_span, y0, env, beam, **kwargs)
    return sol.t, sol.y.T


def simulate_large_pulses(env: SimEnv,
                          beam: ParticleBeam,
                          end_time: float = 1,
                          ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    def simulate(t_span: List[float],
                 y0: NDArray[np.float64],
                 **kwargs: Any,
                 ) -> Any:
        return solve_ivp(derive,
                         t_span,
                         y0,
                         args=(env, beam),
                         jac=jac,
                         method="LSODA",
                         **kwargs,
                         )
    if beam.timer.duty_cycle() > 0.8:
        raise Exception("Too high duty cycle")
    MAX_STEP_ON = beam.timer.on_time / 100
    MAX_STEP_OFF = (beam.timer.period - beam.timer.on_time) / 10
    kwargs = {}
    kwargs["max_step"] = MAX_STEP_ON
    kwargs["first_step"] = 10e-9
    kwargs["atol"] = 1e-6

    t_span = [1e-7, beam.timer.on_time*1.1]
    y0 = env.initial_values()

    end_while = max(beam.time_at_max_dose, end_time)

    sol = simulate(t_span, y0, **kwargs)
    res = sol.y.T
    time = sol.t
    on_time = True

    while t_span[1] < end_while:
        if on_time:
            on_time = False
            t_span = [t_span[1], t_span[1] +
                      beam.timer.period-beam.timer.on_time*1.1]
            kwargs["max_step"] = MAX_STEP_OFF
        else:
            on_time = True
            t_span = [t_span[1], t_span[1]+beam.timer.on_time*1.1]
            kwargs["max_step"] = MAX_STEP_ON

        sol = simulate(t_span, res[-1, :], **kwargs)
        res = np.vstack([res, sol.y.T])
        time = np.append(time, sol.t)

    if end_while < end_time:
        kwargs["max_step"] = (end_time-t_span[1]) / 100
        sol = simulate([t_span[1], end_time], res[-1, :], **kwargs)
        res = np.vstack([res, sol.y.T])
        time = np.append(time, sol.t)
    return time, res


if __name__ == "__main__":
    file = Path(rf"{rp.__path__[0]}\data\reactions.ron")
    beam = ConstantBeam("e", 5)
    env = SimEnv.from_ron_file(file)

    t_span = [1e-9, 500e-6]
    yinit = env.initial_values()
    method = "Radau"
    sol = solve_ivp(derive, t_span, yinit, method,
                    rtol=1e-5, args=(env, beam), jac=partial(jac, env=env, beam=beam))
