#!/usr/bin/env python3
"""
Test radiopyo.chemistry.utils.factories module.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 27/08/2024
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations
import pytest
from numpy.polynomial import Polynomial
import numpy as np

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.chemistry.utils.factories import polynomial_factory


# ----------------------------------- Tests---------------------------------- #
@pytest.mark.parametrize(
    "value, log10, result",
    [
        (1.0, False, Polynomial([1.0])),
        (1.5, True, Polynomial([np.log10(1.5)])),
        ([1.0, 2.0, 2.5], False, Polynomial([1.0, 2.0, 2.5])),
        (
            [1.0, 2.0, 2.5],
            True,
            Polynomial([1.0, 2.0, 2.5]),
        ),
        (Polynomial([1.5]), False, Polynomial([1.5])),
        (Polynomial([1.5]), True, Polynomial([1.5])),
    ],
)
def test_polynomial_factory(value, log10, result):
    assert polynomial_factory(value, log10=log10) == result


def test_polynomial_factory_error():
    with pytest.raises(ValueError):
        polynomial_factory([])
