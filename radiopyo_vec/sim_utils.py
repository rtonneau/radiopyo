# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from itertools import chain
from typing import (
    TYPE_CHECKING,
    Optional,
    Tuple,
)

import numpy as np
import pandas as pd  # type: ignore[import]

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo_vec.physics.beam import ParticleBeam

    from .sim_env import SimEnv


__all__ = [
    "reconstruct_reactions_from_res",
]


def reconstruct_reactions_from_res(t: NDArray,
                                   res: NDArray,
                                   env: SimEnv,
                                   beam: ParticleBeam,
                                   cst_sp: NDArray[np.float64],
                                   ab_mat: Tuple[NDArray[np.int8],
                                                 NDArray[np.int8],
                                                 NDArray[np.float64]],
                                   ki_mat: Tuple[NDArray[np.float64],
                                                 NDArray[np.int8],
                                                 NDArray[np.float32],
                                                 NDArray[np.float32]],
                                   Gi_mat: NDArray[np.float64],
                                   Enzi: Optional[Tuple[NDArray[np.float32],
                                                        NDArray[np.float32],
                                                        NDArray[np.float64],
                                                        NDArray[np.int8],
                                                        ]] = None,
                                   ) -> pd.DataFrame:
    # Compute column names and Pre compute slices to ease access in the dataframe
    # columns
    cols_kr = [kr.as_label() for kr in env.reactions.k_reaction]
    idx_kr = slice(0, len(cols_kr))

    cols_gi = [gr.as_label() for gr in env.reactions.radiolytic]
    idx_gi = slice(idx_kr.stop, idx_kr.stop+len(cols_gi))

    if Enzi is not None:
        cols_ez = [er.as_label() for er in env.reactions.enzymatic]
        idx_ez = slice(idx_gi.stop, idx_gi.stop+len(cols_ez))

    out = pd.DataFrame(columns=list(chain.from_iterable([cols_kr,
                                                         cols_gi,
                                                         cols_ez])),
                       index=t,)

    # Prefill values for radiolytic reactions
    for label, reaction in zip(cols_gi, env.reactions.radiolytic):
        out.loc[:, label] = reaction.kr()

    for ((idx, t_val), y) in zip(enumerate(t), res.T):
        y = y.clip(min=0.0,) * 1e-6
        y = np.append(y, cst_sp)

        # K Reactions
        Rki = np.nan_to_num(np.einsum("j, jik -> ik ", y, ki_mat[2]), nan=1)
        Rki = np.prod(Rki, axis=1) * ki_mat[0]
        out.iloc[idx, idx_kr] = Rki

        # Radiolytic Reactions
        out.iloc[idx, idx_gi] *= beam.at(t_val).dose_rate()

        # Enzimatic reactions
        if Enzi is not None:
            E = np.matmul(y, Enzi[1])
            S = np.matmul(y, Enzi[0])
            MM = Enzi[2][1, :] * E * S / (Enzi[2][0, :] + S)
            out.iloc[idx, idx_ez] = MM

    return out
    
