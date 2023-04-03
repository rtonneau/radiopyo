# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import time
from functools import wraps
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Concatenate,
    Optional,
    ParamSpec,
    Tuple,
    TypedDict,
    TypeVar,
)

import numpy as np
import pandas as pd  # type: ignore[import]
from scipy.integrate import solve_ivp  # type: ignore[import]

# --------------------------------- Internal imports --------------------------------- #
from .chemistry.species import CstSpecies
from .exceptions import NoBeamDefineError
from .input.g_values import GValuesBoscolo
from .ode_solver import derive, jac
from .sim_env import SimEnv
from .sim_utils import reconstruct_KRi_from_res

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from .chemistry.species import SimSpecies
    from .physics.beam import ParticleBeam

# ----------------------------- Define useful decorators ----------------------------- #

P = ParamSpec("P")
R = TypeVar("R")


def with_beam(func: Callable[Concatenate[Simulation, P], R],
              ) -> Callable[Concatenate[Simulation, P], R]:
    @wraps(func)
    def wrapper(self: Simulation, *args: P.args, **kwargs: P.kwargs) -> R:
        if self.beam is None:
            raise NoBeamDefineError
        return func(self, *args, **kwargs)
    return wrapper

# ------------------------------------------------------------------------------------ #


class SimMatrices(TypedDict):
    ab_mat: Tuple[NDArray[np.int8], NDArray[np.int8], NDArray[np.float64]]
    ki_mat: Tuple[NDArray[np.float64], NDArray[np.int8],
                  NDArray[np.float32], NDArray[np.float32]]
    Gi_mat: NDArray[np.float64]
    Enzi: Optional[Tuple[NDArray[np.float32], NDArray[np.float32],
                         NDArray[np.float64], NDArray[np.int8],]]


class Simulation(object):
    file: Path
    env: SimEnv
    mat: SimMatrices
    beam: Optional[ParticleBeam]
    exec_time: float
    g_val: GValuesBoscolo

    def __init__(self,
                 file: str | Path,
                 beam: Optional[ParticleBeam] = None,
                 ) -> None:
        self.file = Path(file)
        self.env = SimEnv.from_ron_file(self.file)
        self.mat = SimMatrices(
            ab_mat=self.env.make_acidBase_matrix(),
            ki_mat=self.env.make_kreaction_matrix(),
            Gi_mat=self.env.make_radiolytic_vector(),
            Enzi=self.env.make_enzymatic_matrix(),
        )
        self.beam = beam
        self.exec_time = -1
        self.g_val = GValuesBoscolo()

    def set_beam(self, beam: ParticleBeam) -> None:
        self.beam = beam

    def change_LET(self, let: float) -> None:
        """Update the radiolytic yields from the simulation env.
        Args:
            let (float): in [kev/µm] and in range [0.14, 150]
        """
        available_sp = list(self.g_val.species())
        for reaction in self.env.reactions.radiolytic:
            sp = reaction.species()
            if sp.label not in available_sp or isinstance(sp, CstSpecies):
                continue

            reaction.update_from_ge(self.g_val.values_at(sp.label, let))
        self.mat["Gi_mat"] = self.env.make_radiolytic_vector()

    @with_beam
    def run(self,
            t_span: Tuple[float, float],
            t_sol: Optional[NDArray] = None,
            ) -> Any:

        t0 = time.time()
        sol = solve_ivp(derive,
                        t_span,
                        self.env.initial_values(),
                        args=(self.beam,
                              self.mat["ab_mat"],
                              self.mat["ki_mat"],
                              self.mat["Gi_mat"],
                              self.mat["Enzi"]),
                        t_eval=t_sol,
                        # jac=jac,
                        atol=1e-15,
                        method="LSODA",
                        first_step=1e-12,
                        )
        t1 = time.time()
        self.exec_time = t1 - t0
        return sol

    def reconstruct_KRi(self, sol: Any) -> pd.DataFrame:
        return reconstruct_KRi_from_res(sol.t,
                                        sol.y,
                                        self.env,
                                        self.mat["ab_mat"],
                                        self.mat["ki_mat"],
                                        )
