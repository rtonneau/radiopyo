#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 27/08/2024
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

import numpy as np

# ----------------------------- Internal imports ---------------------------- #

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    pass

# ------------------------- Expose importable stuffs ------------------------ #
__all__: tp.List = [
    "polynomial_factory",
]

# ---------------------------------- Logging -------------------------------- #

# --------------------------------------------------------------------------- #
#                            DECORATOR DEFINITIONS                            #
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #
def polynomial_factory(
    value: float | tp.List[float] | np.polynomial.Polynomial,
    log10: bool = False,
) -> np.polynomial.Polynomial:
    """
    A factory function to create a numpy polynomial from a given value.

    Args:
        value (float | tp.List[float] | np.polynomial.Polynomial):
            The value to be converted into a numpy polynomial. It can be
            a float, a list of floats or a numpy polynomial.

    Returns:
        np.polynomial.Polynomial:
            A numpy polynomial created from the given value.

    Raises:
        ValueError:
            If the given value is not a valid input for creating a
            numpy polynomial.
    """
    if isinstance(value, np.polynomial.Polynomial):
        return value
    if isinstance(value, (list, tuple, np.ndarray)):
        if len(value) >= 1:
            _coef = [np.log10(elt) if log10 else elt for elt in value]
        else:
            raise ValueError("Cannot make a polynomial with no coeff!")
    else:
        if value <= 0:
            _coef = [-np.inf]
        else:
            _coef = [np.log10(value)] if log10 else [value]

    return np.polynomial.Polynomial(_coef)


# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #
