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

from scipy.integrate import solve_ivp  # type: ignore[import]

# --------------------------------- Internal imports --------------------------------- #
from .chemistry.species import CstSpecies
from .exceptions import NoBeamDefineError
from .input.g_values import GValuesBoscolo
from .ode_solver import derive, jacobian
from .sim_env import SimEnv

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray

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
    AB_mat: Tuple[NDArray[np.int8], NDArray[np.int8], NDArray[np.float64]]
    Ki_mat: Tuple[NDArray[np.float64], NDArray[np.int8],
                  NDArray[np.float32], NDArray[np.float32]]
    Gi_mat: NDArray[np.float64]
    Ez_mat: Optional[Tuple[NDArray[np.float32], NDArray[np.float32],
                           NDArray[np.float64], NDArray[np.int8],]]
    cst_sp: NDArray[np.float64]


class Simulation(object):
    file: Path
    env: SimEnv
    mat: SimMatrices
    beam: ParticleBeam
    exec_time: float
    g_val: GValuesBoscolo

    def __init__(self,
                 file: str | Path,
                 beam: ParticleBeam,
                 ) -> None:
        self.file = Path(file)
        self.env = SimEnv.from_ron_file(self.file)
        self.mat = SimMatrices(
            AB_mat=self.env.make_acidBase_matrix(),
            Ki_mat=self.env.make_kreaction_matrix(),
            Gi_mat=self.env.make_radiolytic_vector(),
            Ez_mat=self.env.make_enzymatic_matrix(),
            cst_sp=self.env.make_cst_species_vector(),
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
            t_eval: Optional[NDArray] = None,
            use_jac: bool = False,
            method: str = "LSODA",
            ) -> Any:

        t0 = time.time()
        sol = solve_ivp(derive,
                        t_span,
                        self.env.initial_values(),
                        args=(self.beam,
                              self.mat["cst_sp"],
                              self.mat["AB_mat"],
                              self.mat["Ki_mat"],
                              self.mat["Gi_mat"],
                              self.mat["Ez_mat"],
                              ),
                        t_eval=t_eval,
                        jac=jacobian if use_jac else None,
                        atol=1e-15,
                        method=method,
                        first_step=1e-12,
                        )
        t1 = time.time()
        self.exec_time = t1 - t0
        return sol
