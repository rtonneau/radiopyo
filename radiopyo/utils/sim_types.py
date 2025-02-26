#!/usr/bin/env python3
"""
The sim_types submodule contains class types hint (TypedDict) in order
to help users to interface to the module.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 08/05/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.chemistry.utils.units import Temperature

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "ConcentrationDict",
    "ConfigDict",
    "KwargsRun",
    "ReactionsDict",
    "SimMatrices",
    "deep_update",
    "default_ConfigDict",
]
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
#                            FUNCTION DEFINITIONS                             #
# --------------------------------------------------------------------------- #


def deep_update(
    original: tp.Dict[str, tp.Any],
    updates: tp.Dict[str, tp.Any],
) -> tp.Dict[str, tp.Any]:
    """Recursively updates a nested dictionary with another dictionary."""
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(original.get(key), dict):
            deep_update(original[key], value)  # Recursive update
        else:
            original[key] = value  # Replace value
    return original


# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


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
    equilibrium: tp.List[tp.Dict]
    arrhenius: tp.List[tp.Dict]
    constant: tp.List[tp.Dict]


class SimParameters(tp.TypedDict, total=False):
    pH: float
    force_cst_pH: bool
    species_order: tp.Dict[str, tp.List[str]]
    temperature: float
    timestamp: float
    t_zero_pulse: float
    force_positive: bool
    variable_GValues: bool


class ConfigDict(tp.TypedDict, total=False):
    includes: tp.Dict[str, tp.List[str]]
    parameters: SimParameters
    concentrations: ConcentrationDict
    reactions: ReactionsDict
    beam: tp.Dict[str, tp.List[tp.Dict]]


class SimMatrices(tp.TypedDict):
    """
    Convenient TypedDict use to leverage type checking and store
    matrices used in ODE solver.
    """

    Ki_mat: tp.Tuple[
        NDArray[np.float64],
        NDArray[np.int8],
        NDArray[np.float32],
        NDArray[np.float32],
    ]
    Gi_mat: NDArray[np.float64]
    Ez_mat: tp.Optional[
        tp.Tuple[
            NDArray[np.float32],
            NDArray[np.float32],
            NDArray[np.float64],
            NDArray[np.int8],
        ]
    ]
    cst_sp: NDArray[np.float64]


def default_SimParameters() -> SimParameters:
    return SimParameters(
        pH=7.0,
        force_cst_pH=False,
        species_order={},
        temperature=Temperature._DEFAULT,
        timestamp=0.0,
        t_zero_pulse=1e-7,
        force_positive=True,
        variable_GValues=False,
    )


def default_ConfigDict() -> ConfigDict:
    """ """
    return {
        "includes": {},
        "parameters": default_SimParameters(),
        "concentrations": {
            "fixed": {},
            "initial": {},
        },
        "beam": {},
        "reactions": {
            "acid_base": [],
            "k_reaction": [],
            "equilibrium": [],
            "enzymatic": [],
            "radiolytic": {},
            "arrhenius": [],
            "constant": [],
        },
    }
