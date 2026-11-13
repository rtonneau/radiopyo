#!/usr/bin/env python3
"""
Fixtures of the Rust backend tests (radiopyo.rust, ADR 0010).

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy
import typing as tp

import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.simulation.unit_cell import UnitCell

CONFIGS = {
    "pure_water": rp.PURE_WATER_CONFIG_FILE,
    "labarbe": rp.LABARBE_CONFIG_FILE,
}

_TEMPLATES: tp.Dict[str, UnitCell] = {}


def template(name: str) -> UnitCell:
    """A fresh copy of the cell built from CONFIGS[name] (built once)."""
    if name not in _TEMPLATES:
        _TEMPLATES[name] = UnitCell.from_toml(CONFIGS[name])
    return copy.deepcopy(_TEMPLATES[name])


@pytest.fixture
def labarbe() -> UnitCell:
    return template("labarbe")


@pytest.fixture
def pure_water() -> UnitCell:
    return template("pure_water")
