# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Optional,
    Tuple,
)

import numpy as np
import pandas as pd  # type: ignore[import]
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


def matprint(mat: NDArray, fmt: str = "g") -> None:
    col_maxes = [max([len(("{:"+fmt+"}").format(x))
                     for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
        print("")


def jac(t: float,
        y: NDArray[np.float64],
        env: SimEnv,
        beam: ParticleBeam
        ) -> NDArray[np.float64]:
    pass


def derive(t: float,
           y: NDArray[np.float64],
           ab_mat: Tuple[NDArray[np.int8],
                         NDArray[np.int8],
                         NDArray[np.float64]],
           ki_mat: Tuple[int,
                         NDArray[np.float64],
                         NDArray[np.int8]],
           Gi_mat: NDArray[np.float64],
           beam: ParticleBeam,
           Enzi: Optional[Tuple[NDArray[np.int8],
                                NDArray[np.int8],
                                NDArray[np.float64],
                                ]] = None,
           ) -> NDArray[np.float64]:

    #print(f"\nTime is {t}")
    # Force negative cc values to 0
    np.set_printoptions(precision=1)
    print(f"\ny(t={t:.2g}) = {y*1e-6} (before clipping)")
    y = y.clip(min=0.0,) * 1e-6
    #print(f"===> Before A/B: {y}")
    # First step ==> Handling acid/base equilibrium
    # ABy, ABi, ABo => ab_mat
    y = ab_mat[0]*y + np.matmul(ab_mat[2], np.matmul(y, ab_mat[1]).T)
    print(f"===> After A/B: {y}")

    # Build u vector => [1, y(t), y(t),...] n times y(t) with n = max
    # stoichiometry of reactants specie s (max_stoi).
    # max_stoi, Ri, Riy => ki_mat
    #u = np.append([1], np.tile(y, ki_mat[0]))
    #print(f"y is: {y}")
    Rki = np.nan_to_num(np.einsum("j, jik -> ik ", y, ki_mat[2]), nan=1)
    # matprint(Rki)

    # Compute dydt, first only with k reactions term:
    dydt = np.matmul(ki_mat[1], (np.prod(Rki, axis=1) * ki_mat[0]).T)
    print(f"===> dydt From kreactions : {dydt[3:13]}")

    # Add Radiolytic term ==> Need deposition rate!
    # => Only term function of time... probably need refactorization => TODO
    dydt_radio = Gi_mat * beam.at(t).dose_rate()

    print(f"===> dydt From Radiolytics : {dydt_radio}")

    dydt += dydt_radio
    # Add Enzymatic term
    if Enzi is not None:
        # En_enzy, En_subs, En_kmki => Enzi[0], Enzi[1], Enzi[2]
        # Extract substrate cc of each enzymatic reaction.
        S = np.matmul(y, Enzi[1])
        beta = np.matmul(Enzi[2][0, :] * Enzi[2][1, :], Enzi[0])
        dydt += ((Enzi[2][0, :] + S) / beta)
    #print(f"===> dydt From radiolytics: {dydt[3:]}")
    return dydt * 1e6


def reconstruct_Ri_from_res(time: NDArray,
                            res: NDArray,
                            env: SimEnv,
                            ab_mat: Tuple[NDArray[np.int8],
                                          NDArray[np.int8],
                                          NDArray[np.float64]],
                            ki_mat: Tuple[int,
                                          NDArray[np.float64],
                                          NDArray[np.int8]],
                            Gi_mat: NDArray[np.float64],
                            beam: ParticleBeam,
                            Enzi: Optional[Tuple[NDArray[np.int8],
                                                 NDArray[np.int8],
                                                 NDArray[np.float64],
                                                 ]] = None,
                            ) -> pd.DataFrame:
    out = pd.DataFrame(columns=[f"R_{idx:d}" for idx in
                                range(len(env.reactions.k_reaction))],
                       index=time,
                       )
    for ((idx, t), y) in zip(enumerate(time), res.T):
        y = y.clip(min=0.0,) * 1e-6

        # Acid/Base
        y = ab_mat[0]*y + np.matmul(ab_mat[2], np.matmul(y, ab_mat[1]).T)

        # K Reactions
        Rki = np.nan_to_num(np.einsum("j, jik -> ik ", y, ki_mat[2]), nan=1)
        Rki = np.prod(Rki, axis=1) * ki_mat[0]
        out.iloc[idx, :] = Rki
    return out


if __name__ == "__main__":
    file = Path(rf"{rp.__path__[0]}\data\reactions.ron")
    beam = ConstantBeam("e", 5)
    env = SimEnv.from_ron_file(file)

    t_span = [1e-9, 500e-6]
    yinit = env.initial_values()
    method = "Radau"
    sol = solve_ivp(derive, t_span, yinit, method,
                    rtol=1e-5, args=(env, beam))
