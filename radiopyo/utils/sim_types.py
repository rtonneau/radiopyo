#!/usr/bin/env python3
"""
The sim_types submodule contains class types hint (TypedDict) in order to help users to
interface to the module.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 08/05/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp

# --------------------------------- Internal imports --------------------------------- #


# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = [
    "KwargsRun",
    "ConfigDict",
    "ReactionsDict",
    "ConcentrationDict",
    "default_ConfigDict",
    "SimMatrices",
]
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class KwargsRun(tp.TypedDict):
    t_span: tp.Tuple[float, float]
    t_eval: tp.Optional[NDArray[np.float64]]
    use_jac: bool
    method: str
    atol: tp.Optional[float]
    rtol: tp.Optional[float]
    max_step: tp.Optional[float]
    y0: tp.Optional[NDArray[np.float64]]


class ConcentrationDict(tp.TypedDict):
    fixed: tp.Dict[str, float]
    initial: tp.Dict[str, float]


class ReactionsDict(tp.TypedDict):
    radiolytic: tp.Dict
    acid_base: tp.List[tp.Dict]
    k_reaction: tp.List[tp.Dict]
    enzymatic: tp.List[tp.Dict]


class ConfigDict(tp.TypedDict):
    includes: tp.Dict[str, tp.List[str]]
    bio_param: tp.Dict[str, float]
    concentrations: ConcentrationDict
    reactions: ReactionsDict
    beam: tp.Dict[str, tp.List[tp.Dict]]


class SimMatrices(tp.TypedDict):
    """Convenient TypedDict use to leverage type checking and store matrices used in ODE
    solver."""
    AB_mat: tp.Tuple[NDArray[np.int8], NDArray[np.int8], NDArray[np.float64]]
    Ki_mat: tp.Tuple[NDArray[np.float64], NDArray[np.int8],
                     NDArray[np.float32], NDArray[np.float32]]
    Gi_mat: NDArray[np.float64]
    Ez_mat: tp.Optional[tp.Tuple[NDArray[np.float32], NDArray[np.float32],
                                 NDArray[np.float64], NDArray[np.int8],]]
    cst_sp: NDArray[np.float64]
    O2_mat: tp.Optional[tp.Tuple[NDArray[np.float64], NDArray[np.float64]]]


def default_ConfigDict() -> ConfigDict:
    """ """
    return {"includes": {},
            "bio_param": {},
            "concentrations": {
                "fixed": {},
                "initial": {},
    },
        "beam": {},
        "reactions": {
            "acid_base": [],
            "k_reaction": [],
            "enzymatic": [],
            "radiolytic": {},
    },
    }
