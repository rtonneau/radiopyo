# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

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
    "reconstruct_KRi_from_res",
    "reconstruct_ERi_from_res",
]


def reconstruct_KRi_from_res(time: NDArray,
                             res: NDArray,
                             env: SimEnv,
                             ab_mat: Tuple[NDArray[np.int8],
                                           NDArray[np.int8],
                                           NDArray[np.float64]],
                             ki_mat: Tuple[NDArray[np.float64],
                                           NDArray[np.int8],
                                           NDArray[np.float32],
                                           NDArray[np.float32]],
                             ) -> pd.DataFrame:
    out = pd.DataFrame(columns=[f"KR_{idx:d}" for idx in
                                range(len(env.reactions.k_reaction))],
                       index=time,
                       )
    for ((idx, _), y) in zip(enumerate(time), res.T):
        y = y.clip(min=0.0,) * 1e-6

        # Acid/Base
        #y = ab_mat[0]*y + np.matmul(ab_mat[2], np.matmul(y, ab_mat[1]).T)

        # K Reactions
        Rki = np.nan_to_num(np.einsum("j, jik -> ik ", y, ki_mat[2]), nan=1)
        Rki = np.prod(Rki, axis=1) * ki_mat[0]
        out.iloc[idx, :] = Rki
    return out
