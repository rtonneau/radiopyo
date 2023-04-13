# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import time
from contextlib import suppress
from itertools import repeat
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Iterable,
    List,
    NamedTuple,
    Optional,
    Tuple,
    TypedDict,
)

import numpy as np
from scipy.integrate import solve_ivp  # type: ignore[import]

# --------------------------------- Internal imports --------------------------------- #
from .chemistry.species import CstSpecies
from .input.g_values import GValuesBoscolo
from .ode_solver import ODEResult, ScipyODEResult, derive, jacobian
from .sim_env import SimEnv

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from .physics.beam import ParticleBeam

# ----------------------------- Define useful decorators ----------------------------- #

# ------------------------------------------------------------------------------------ #


class SimMatrices(TypedDict):
    AB_mat: Tuple[NDArray[np.int8], NDArray[np.int8], NDArray[np.float64]]
    Ki_mat: Tuple[NDArray[np.float64], NDArray[np.int8],
                  NDArray[np.float32], NDArray[np.float32]]
    Gi_mat: NDArray[np.float64]
    Ez_mat: Optional[Tuple[NDArray[np.float32], NDArray[np.float32],
                           NDArray[np.float64], NDArray[np.int8],]]
    cst_sp: NDArray[np.float64]


class UnitCell(object):
    FIRST_STEP: float = 1e-12

    file: Path
    env: SimEnv
    mat: SimMatrices
    beam: ParticleBeam
    last_exec_time: float
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
        self.last_exec_time = -1
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

    def run(self,
            t_span: Tuple[float, float],
            t_eval: Optional[NDArray[np.float64]] = None,
            use_jac: bool = False,
            method: str = "LSODA",
            atol: float = 1e-8,
            max_step: float = np.inf,
            y0: Optional[NDArray[np.float64]] = None,
            ) -> ScipyODEResult:

        if y0 is None:
            y0 = self.env.initial_values()

        first_step = self.FIRST_STEP
        if first_step >= max_step:
            first_step = max_step / 10

        t0 = time.time()
        sol: ScipyODEResult = solve_ivp(derive,
                                        t_span,
                                        y0,
                                        args=(self.beam,
                                              self.mat["cst_sp"],
                                              self.mat["AB_mat"],
                                              self.mat["Ki_mat"],
                                              self.mat["Gi_mat"],
                                              self.mat["Ez_mat"],
                                              ),
                                        t_eval=t_eval,
                                        jac=jacobian if use_jac else None,
                                        atol=atol,
                                        method=method,
                                        max_step=max_step,
                                        first_step=first_step,
                                        )
        t1 = time.time()
        self.last_exec_time = t1 - t0
        return sol

    def run_chunked(self,
                    t_chunks: NDArray[np.float64],
                    t_eval: Optional[NDArray[np.float64]] = None,
                    use_jac: bool = False,
                    method: str = "LSODA",
                    atol: float = 1e-10,
                    max_step: float | Iterable[float] = np.inf,
                    ) -> ODEResult:
        sol = ODEResult()

        # Check up max_step
        if isinstance(max_step, float):
            max_step = repeat(max_step, len(t_chunks)-1)
        elif t_chunks.shape[0] != len(list(max_step)):
            raise ValueError(f"t_chunks (len={len(t_chunks)}) and max_step "
                             f"(len={len(list(max_step))}) must have the same "
                             f"length if max_step is a list or and array")

        # Checkup t_eval
        iter_t_eval: Iterable | List
        if t_eval is None:
            iter_t_eval = repeat(None, len(t_chunks)-1)
        else:
            t_eval = np.array(t_eval)
            iter_t_eval = []
            for idx, stop in enumerate(t_chunks[1:], 1):
                to_push = []
                with suppress(IndexError):
                    while t_eval[0] < stop:
                        to_push.append(t_eval[0])
                        # remove first element
                        t_eval = np.delete(t_eval, (0))

                # if nothing to add in this interval, at least at startup point
                if len(to_push) == 0:
                    to_push.append(t_chunks[idx-1])
                    to_push.append(t_chunks[idx-1]+(stop-t_chunks[idx-1])/2)
                iter_t_eval.append(np.array(to_push))

        # Loop through
        for idx, (stop, step_size, _t_eval) in enumerate(zip(t_chunks[1:],
                                                             max_step,
                                                             iter_t_eval,),
                                                         start=0):
            res = self.run(t_span=(t_chunks[idx], stop),
                           use_jac=use_jac,
                           method=method,
                           atol=atol,
                           max_step=step_size,
                           t_eval=_t_eval,
                           y0=sol.y[:, -1] if len(sol) > 0 else None,
                           )
            sol.append_sol(res, exec_time=self.last_exec_time)
        return sol

    def prepare_chunked_run(self,
                            t_chunks: NDArray[np.float64],
                            use_jac: bool = False,
                            method: str = "LSODA",
                            atol: float = 1e-10,
                            max_step: float | Iterable[float] = np.inf,
                            ) -> SimChunksIter:

        # Check up max_step
        if isinstance(max_step, float):
            max_step = repeat(max_step, len(t_chunks)-1)
        elif t_chunks.shape[0] != len(list(max_step)):
            raise ValueError(f"t_chunks (len={len(t_chunks)}) and max_step "
                             f"(len={len(list(max_step))}) must have the same "
                             f"length if max_step is a list or and array")

        chunk_list: List[SimChunk] = []
        for idx, (stop, step_size) in enumerate(zip(t_chunks[1:],
                                                    max_step,),
                                                start=0):
            chunk = SimChunk(start=t_chunks[idx],
                             stop=stop,
                             max_step_size=step_size)
            chunk_list.append(chunk)
        return SimChunksIter(sim=self,
                             chunks=chunk_list,
                             atol=atol,
                             use_jac=use_jac,
                             method=method,)


class SimChunk(NamedTuple):
    start: float
    stop: float
    max_step_size: float


class SimChunksIter():
    sim: UnitCell
    chunks: List[SimChunk]
    sol: ODEResult
    current: int
    atol: float
    use_jac: bool
    method: str

    def __init__(self, sim: UnitCell,
                 chunks: List[SimChunk],
                 atol: float = 1e-10,
                 use_jac: bool = True,
                 method: str = "LSODA",
                 ):
        self.sim = sim
        self.chunks = chunks
        self.sol = ODEResult()
        self.atol = atol
        self.use_jac = use_jac
        self.method = method
        self.current = -1

    def __iter__(self):
        return self

    def __next__(self):
        self.current += 1
        try:
            chunk = self.chunks[self.current]
        except IndexError:
            raise StopIteration from None
        res = self.sim.run(t_span=(chunk.start, chunk.stop),
                           use_jac=self.use_jac,
                           method=self.method,
                           atol=self.atol,
                           max_step=chunk.max_step_size,
                           y0=self.sol.y[:, -1] if self.current > 0 else None,
                           )
        self.sol.append_sol(res, exec_time=self.sim.last_exec_time)
        return self.sol

    def __len__(self):
        return len(self.chunks)
